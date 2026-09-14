#!/usr/bin/env python3
"""SCION-X burst detector and 1 Hz GFSK carrier search.

The GNU Radio flowgraph feeds normalized complex64 IQ samples into
``StreamingBurstScanner``.  For each detected burst the scanner:

1. computes an FM/GFSK discriminator in Hz;
2. applies a 6.5 kHz post-demodulation low-pass filter;
3. searches every requested carrier hypothesis and several fractional symbol
   phases for complete HDLC candidates and checks their CRC;
4. retains signed soft symbol values and, after a hard CRC failure, runs bounded
   stuffing-aware beam plus CRC-syndrome repair;
5. writes each selected frame's *raw serial bits* (flags and stuffing retained),
   packed LSB-first into a separate ``Frame_NNN/GNU_GFSK_Raw_HEX.txt``.

Frequency translation followed by a zero-threshold slicer is mathematically
equivalent to slicing the discriminator at the candidate carrier frequency.
This lets a 40001-point, 1 Hz scan run without constructing 40001 DSP chains.
"""

from __future__ import annotations

import argparse
import csv
import importlib
import importlib.util
import json
import sys
from collections import deque
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any, Iterable

import numpy as np
from scipy.signal import fftconvolve, firwin


FLAG_BITS = np.asarray([0, 1, 1, 1, 1, 1, 1, 0], dtype=np.uint8)
PREAMBLE_BITS = np.tile(FLAG_BITS, 4)
SCIONX_FRAME_BYTES = 274
SCIONX_FRAME_BITS = SCIONX_FRAME_BYTES * 8
SCIONX_HEADER_BYTES = bytes.fromhex(
    "84 9C 60 86 AA 40 60 84 9C 60 A6 86 B0 E1 03 F0"
)


def _default_soft_repair_library() -> Path:
    """Return a bundled or repository-local SCION-X decoder package root."""
    module_dir = Path(__file__).resolve().parent
    candidates = [module_dir / "soft_repair"]
    try:
        candidates.append(
            module_dir.parents[1]
            / "Reference"
            / "beacon analysis"
            / "scionx-decode-writeup-main"
        )
    except IndexError:
        pass
    for root in candidates:
        if (root / "scionx" / "advanced_decoder.py").is_file():
            return root
    # Preserve the portable location in the eventual error message.
    return candidates[0]


@lru_cache(maxsize=4)
def _load_advanced_decoder(library_root: str | None = None) -> Any:
    """Load the repair package under a private name to avoid import collisions."""
    root = Path(library_root).resolve() if library_root else _default_soft_repair_library()
    package_dir = root / "scionx"
    package_init = package_dir / "__init__.py"
    advanced_file = package_dir / "advanced_decoder.py"
    if not package_init.is_file() or not advanced_file.is_file():
        raise FileNotFoundError(
            f"SCION-X advanced decoder package was not found below {root}"
        )

    package_name = "_scionx_burst_soft_repair_" + str(abs(hash(str(root))))
    if package_name not in sys.modules:
        spec = importlib.util.spec_from_file_location(
            package_name,
            package_init,
            submodule_search_locations=[str(package_dir)],
        )
        if spec is None or spec.loader is None:
            raise ImportError(f"Cannot load SCION-X decoder package from {package_init}")
        package = importlib.util.module_from_spec(spec)
        sys.modules[package_name] = package
        spec.loader.exec_module(package)
    return importlib.import_module(f"{package_name}.advanced_decoder")


def crc16_x25(data: bytes) -> int:
    """Return the AX.25/HDLC CRC-16/X-25 value for *data*."""
    crc = 0xFFFF
    for value in data:
        crc ^= value
        for _ in range(8):
            crc = (crc >> 1) ^ 0x8408 if crc & 1 else crc >> 1
    return (~crc) & 0xFFFF


def fcs_ok(frame: bytes) -> bool:
    """Check the little-endian two-byte AX.25 FCS at the end of a frame."""
    if len(frame) <= 2:
        return False
    expected = crc16_x25(frame[:-2])
    received = frame[-2] | (frame[-1] << 8)
    return expected == received


def pack_lsb_first(bits: np.ndarray) -> bytes:
    """Pack chronological serial bits into bytes, first bit in bit position 0."""
    bits = np.asarray(bits, dtype=np.uint8).ravel()
    if bits.size % 8:
        bits = np.pad(bits, (0, 8 - bits.size % 8))
    if bits.size == 0:
        return b""
    octets = bits.reshape(-1, 8)
    weights = (1 << np.arange(8, dtype=np.uint16))[None, :]
    return bytes(np.sum(octets.astype(np.uint16) * weights, axis=1).astype(np.uint8))


def find_flag_positions(bits: np.ndarray) -> np.ndarray:
    """Return every bit index at which an HDLC 0x7E flag begins."""
    bits = np.asarray(bits, dtype=np.uint8).ravel()
    if bits.size < 8:
        return np.empty(0, dtype=np.int64)
    windows = np.lib.stride_tricks.sliding_window_view(bits, 8)
    return np.flatnonzero(np.all(windows == FLAG_BITS, axis=1))


def destuff_bits(bits: np.ndarray) -> tuple[np.ndarray, bool]:
    """Remove HDLC stuffed zeroes from bits between two flags."""
    output: list[int] = []
    ones = 0
    index = 0
    source = np.asarray(bits, dtype=np.uint8).ravel()
    while index < source.size:
        bit = int(source[index])
        output.append(bit)
        if bit:
            ones += 1
            if ones > 6:
                return np.asarray(output, dtype=np.uint8), False
        else:
            ones = 0
        index += 1
        if ones == 5:
            if index >= source.size or int(source[index]) != 0:
                return np.asarray(output, dtype=np.uint8), False
            index += 1
            ones = 0
    return np.asarray(output, dtype=np.uint8), True


@dataclass
class FrameCandidate:
    opening_flag_index: int
    closing_flag_index: int
    frame: bytes
    crc_ok: bool
    boundary_mode: str = "STRICT"
    header_bit_errors: int = 0
    closing_flag_bit_errors: int = 0
    destuff_violations: int = 0
    soft_repair_summary: dict[str, object] | None = None


def _bytes_to_lsb_bits(data: bytes) -> np.ndarray:
    """Expand bytes into chronological AX.25/HDLC LSB-first serial bits."""
    if not data:
        return np.empty(0, dtype=np.uint8)
    values = np.frombuffer(data, dtype=np.uint8)
    return ((values[:, None] >> np.arange(8, dtype=np.uint8)) & 1).astype(
        np.uint8
    ).reshape(-1)


SCIONX_HEADER_BITS = _bytes_to_lsb_bits(SCIONX_HEADER_BYTES)


def _header_bit_errors(frame: bytes) -> int:
    """Return Hamming distance from the fixed 16-byte SCION-X AX.25 header."""
    available = min(len(frame), len(SCIONX_HEADER_BYTES))
    if available == 0:
        return len(SCIONX_HEADER_BYTES) * 8
    actual = _bytes_to_lsb_bits(frame[:available])
    expected = SCIONX_HEADER_BITS[: actual.size]
    missing = (len(SCIONX_HEADER_BYTES) - available) * 8
    return int(np.count_nonzero(actual != expected)) + missing


def _preamble_opening_flags(bits: np.ndarray, minimum_flags: int = 4) -> list[int]:
    """Return the last exact flag of each contiguous preamble flag run."""
    flags = find_flag_positions(bits)
    if flags.size < minimum_flags:
        return []
    openings: list[int] = []
    run_start = 0
    for index in range(1, flags.size + 1):
        run_ended = index == flags.size or int(flags[index] - flags[index - 1]) != 8
        if not run_ended:
            continue
        if index - run_start >= minimum_flags:
            openings.append(int(flags[index - 1]))
        run_start = index
    return openings


def _destuff_fixed_length(
    bits: np.ndarray,
    first_data_bit: int,
    output_bit_count: int,
) -> tuple[np.ndarray, int, int, bool]:
    """Tolerantly destuff exactly *output_bit_count* bits without padding.

    A bit immediately following five data ones is always treated as the
    stuffing slot.  A zero is normal; a one is discarded as a hypothesized
    corrupted stuffed zero and increments the violation count.  The returned
    raw cursor therefore points at the expected closing flag position.
    """
    source = np.asarray(bits, dtype=np.uint8).ravel()
    output: list[int] = []
    raw_index = int(first_data_bit)
    ones = 0
    violations = 0
    while len(output) < output_bit_count and raw_index < source.size:
        bit = int(source[raw_index])
        raw_index += 1
        output.append(bit)
        if bit:
            ones += 1
        else:
            ones = 0
        if ones == 5:
            if raw_index >= source.size:
                return np.asarray(output, dtype=np.uint8), raw_index, violations, False
            if int(source[raw_index]) != 0:
                violations += 1
            raw_index += 1
            ones = 0
    complete = len(output) == output_bit_count
    return np.asarray(output, dtype=np.uint8), raw_index, violations, complete


def recover_frame_boundaries(
    bits: np.ndarray,
    max_header_bit_errors: int = 16,
    max_closing_flag_bit_errors: int = 2,
    max_destuff_violations: int = 8,
) -> list[FrameCandidate]:
    """Recover fixed-length SCION-X frames after an exact flag preamble.

    This is deliberately conservative: no bits are padded and no frame bytes
    are accepted unless 2192 destuffed bits exist, the expected closing flag
    is within the configured Hamming distance, and the fixed AX.25 header is
    sufficiently close.  CRC is evaluated last and may fail; such a frame is
    still exported with an explicit untrusted status.
    """
    source = np.asarray(bits, dtype=np.uint8).ravel()
    recovered: list[FrameCandidate] = []
    for opening in _preamble_opening_flags(source):
        destuffed, closing, violations, complete = _destuff_fixed_length(
            source, opening + 8, SCIONX_FRAME_BITS
        )
        if not complete or closing + 8 > source.size:
            continue
        closing_errors = int(
            np.count_nonzero(source[closing:closing + 8] != FLAG_BITS)
        )
        if closing_errors > max_closing_flag_bit_errors:
            continue
        frame = pack_lsb_first(destuffed)
        header_errors = _header_bit_errors(frame)
        if header_errors > max_header_bit_errors or violations > max_destuff_violations:
            continue
        recovered.append(
            FrameCandidate(
                opening_flag_index=opening,
                closing_flag_index=closing,
                frame=frame,
                crc_ok=fcs_ok(frame),
                boundary_mode="RECOVERED",
                header_bit_errors=header_errors,
                closing_flag_bit_errors=closing_errors,
                destuff_violations=violations,
            )
        )
    return recovered


def diagnostic_frame_boundaries(bits: np.ndarray) -> list[FrameCandidate]:
    """Create untrusted 274-byte candidates after a reliable flag preamble.

    This fallback deliberately relaxes the closing-flag, header and stuffing
    quality gates, but it never pads missing data bits.  It exists only to let
    an operator inspect best-effort telemetry after a FLAGS_ONLY result.
    """
    source = np.asarray(bits, dtype=np.uint8).ravel()
    diagnostics: list[FrameCandidate] = []
    for opening in _preamble_opening_flags(source):
        destuffed, closing, violations, complete = _destuff_fixed_length(
            source, opening + 8, SCIONX_FRAME_BITS
        )
        if not complete:
            continue
        available_closing = source[closing:min(closing + 8, source.size)]
        closing_errors = 8 - int(available_closing.size)
        if available_closing.size:
            closing_errors += int(
                np.count_nonzero(
                    available_closing != FLAG_BITS[: available_closing.size]
                )
            )
        frame = pack_lsb_first(destuffed)
        diagnostics.append(
            FrameCandidate(
                opening_flag_index=opening,
                closing_flag_index=closing,
                frame=frame,
                crc_ok=fcs_ok(frame),
                boundary_mode="FLAGS_ONLY_DIAGNOSTIC",
                header_bit_errors=_header_bit_errors(frame),
                closing_flag_bit_errors=closing_errors,
                destuff_violations=violations,
            )
        )
    return diagnostics


def frame_candidates(bits: np.ndarray) -> list[FrameCandidate]:
    """Decode candidates bounded by consecutive HDLC flags."""
    flags = find_flag_positions(bits)
    frames: list[FrameCandidate] = []
    for left, right in zip(flags[:-1], flags[1:]):
        # Adjacent flags are an idle/preamble run, not an empty frame.
        if right <= left + 8:
            continue
        raw = bits[left + 8:right]
        destuffed, valid = destuff_bits(raw)
        if not valid or destuffed.size % 8:
            continue
        frame = pack_lsb_first(destuffed)
        if len(frame) < 3:
            continue
        frames.append(
            FrameCandidate(
                int(left),
                int(right),
                frame,
                fcs_ok(frame),
                header_bit_errors=_header_bit_errors(frame),
            )
        )
    return frames


def _pattern_count_by_threshold(
    values: np.ndarray,
    pattern: np.ndarray,
    scan_min_hz: int,
    scan_max_hz: int,
    scan_step_hz: int,
) -> np.ndarray:
    """Count pattern windows for every integer carrier threshold efficiently."""
    count = (scan_max_hz - scan_min_hz) // scan_step_hz + 1
    pattern = np.asarray(pattern, dtype=np.uint8).ravel()
    if values.size < pattern.size:
        return np.zeros(count, dtype=np.int32)
    windows = np.lib.stride_tricks.sliding_window_view(values, pattern.size)
    zero_ceiling = np.max(windows[:, pattern == 0], axis=1)
    one_floor = np.min(windows[:, pattern == 1], axis=1)

    # bit = discriminator > carrier.  A flag exists when both zero locations
    # are <= carrier and all six one locations are > carrier.
    low = np.ceil((zero_ceiling - scan_min_hz) / scan_step_hz).astype(np.int64)
    high = np.ceil((one_floor - scan_min_hz) / scan_step_hz).astype(np.int64) - 1
    low = np.maximum(low, 0)
    high = np.minimum(high, count - 1)
    valid = high >= low
    low = low[valid]
    high = high[valid]

    difference = np.zeros(count + 1, dtype=np.int32)
    np.add.at(difference, low, 1)
    np.add.at(difference, high + 1, -1)
    return np.cumsum(difference[:-1])


def _exact_pattern_count(bits: np.ndarray, pattern: np.ndarray) -> int:
    bits = np.asarray(bits, dtype=np.uint8).ravel()
    pattern = np.asarray(pattern, dtype=np.uint8).ravel()
    if bits.size < pattern.size:
        return 0
    windows = np.lib.stride_tricks.sliding_window_view(bits, pattern.size)
    return int(np.count_nonzero(np.all(windows == pattern, axis=1)))


@dataclass
class TimingHypothesis:
    preamble_score: int
    flag_score: int
    phase_fraction: float
    phase_index: int
    carrier_hz: int
    symbol_rate: float
    polarity: int
    values: np.ndarray
    timing_mode: str = "FIXED"
    reported_carrier_hz: int | None = None
    filter_method: str = "RAW"


@dataclass
class HypothesisResult:
    hypothesis: TimingHypothesis
    frames: list[FrameCandidate]
    recovered_frames: list[FrameCandidate]
    diagnostic_frames: list[FrameCandidate]
    bits: np.ndarray
    key: tuple[object, ...]

    @property
    def trusted_frames(self) -> list[FrameCandidate]:
        strict = [
            frame for frame in self.frames if len(frame.frame) == SCIONX_FRAME_BYTES
        ]
        return strict if strict else self.recovered_frames


@dataclass
class SoftRepairSelection:
    """One uniquely repaired frame and its original soft observation."""

    hypothesis: TimingHypothesis
    frame: FrameCandidate
    original_bits: np.ndarray
    corrected_bits: np.ndarray
    decision_values: np.ndarray
    signed_soft_values: np.ndarray
    soft_distances: np.ndarray
    status: str
    reason: str


def _hard_and_soft_decisions(
    hypothesis: TimingHypothesis,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return hard bits plus raw, signed-normalized, and magnitude soft values."""
    decision_values = np.asarray(hypothesis.values, dtype=np.float64).reshape(-1)
    centered = (
        decision_values - float(hypothesis.carrier_hz)
    ) * float(hypothesis.polarity)
    residual = np.abs(centered)
    if residual.size:
        scale = max(
            float(np.percentile(residual, 75.0)),
            float(np.std(centered)) * 0.25,
            1e-9,
        )
    else:
        scale = 1.0
    signed_soft = np.clip(centered / scale, -6.0, 6.0)
    signed_soft = np.nan_to_num(signed_soft, nan=0.0, posinf=6.0, neginf=-6.0)
    distances = np.abs(signed_soft)
    bits = (signed_soft > 0.0).astype(np.uint8)
    return bits, decision_values.copy(), signed_soft, distances


def _soft_repair_audit_metrics(
    repair_summary: dict[str, object] | None,
) -> dict[str, object]:
    """Return stable per-frame repair counts for text/CSV audit outputs."""
    if repair_summary is None:
        return {
            "applied": False,
            "bit_flip_count": 0,
            "structural_flip_count": 0,
            "syndrome_flip_count": 0,
            "global_flip_indices": [],
            "method": "",
        }
    raw_indices = repair_summary.get(
        "global_flip_indices", repair_summary.get("flip_indices", [])
    )
    if not isinstance(raw_indices, (list, tuple)):
        raw_indices = []
    indices = sorted({int(index) for index in raw_indices})
    structural = repair_summary.get("beam_flip_indices", [])
    if not isinstance(structural, (list, tuple)):
        structural = []
    return {
        "applied": bool(indices),
        "bit_flip_count": len(indices),
        "structural_flip_count": len({int(index) for index in structural}),
        "syndrome_flip_count": int(
            repair_summary.get("syndrome_flip_count", 0) or 0
        ),
        "global_flip_indices": indices,
        "method": str(repair_summary.get("repair_method", "")),
    }


def _attempt_advanced_soft_repair(
    hypotheses: Iterable[TimingHypothesis],
    *,
    library_root: str | None,
    maximum_header_bit_errors: int,
    maximum_observations: int,
    repair_effort: str,
) -> tuple[SoftRepairSelection | None, str, str]:
    """Run bounded stuffing-aware soft/CRC repair on promising HDLC routes."""
    try:
        advanced = _load_advanced_decoder(library_root)
    except Exception as error:
        return None, "ADVANCED_DECODER_UNAVAILABLE", str(error)

    # First reduce the potentially large 1 Hz/baud acquisition bank using its
    # exact preamble score.  Header distance then ranks the remaining openings.
    unique_hypotheses: dict[tuple[object, ...], TimingHypothesis] = {}
    for hypothesis in hypotheses:
        if hypothesis.preamble_score < 1:
            continue
        key = (
            hypothesis.filter_method,
            hypothesis.timing_mode,
            hypothesis.phase_index,
            hypothesis.phase_fraction,
            hypothesis.carrier_hz,
            hypothesis.reported_carrier_hz,
            hypothesis.symbol_rate,
            hypothesis.polarity,
        )
        previous = unique_hypotheses.get(key)
        if previous is None or hypothesis.flag_score > previous.flag_score:
            unique_hypotheses[key] = hypothesis
    preselected = sorted(
        unique_hypotheses.values(),
        key=lambda item: (item.preamble_score, item.flag_score),
        reverse=True,
    )[:64]

    records: list[dict[str, object]] = []
    for hypothesis in preselected:
        bits, decision_values, signed_soft, distances = _hard_and_soft_decisions(
            hypothesis
        )
        for opening in _preamble_opening_flags(bits):
            first_data_bit = opening + len(FLAG_BITS)
            minimum_raw_bits = SCIONX_FRAME_BITS + len(FLAG_BITS)
            if bits.size - first_data_bit < minimum_raw_bits:
                continue
            header_bits, _, _, complete = _destuff_fixed_length(
                bits, first_data_bit, SCIONX_HEADER_BITS.size
            )
            if not complete:
                continue
            header_errors = int(np.count_nonzero(header_bits != SCIONX_HEADER_BITS))
            if header_errors > maximum_header_bit_errors:
                continue
            records.append({
                "hypothesis": hypothesis,
                "opening": opening,
                "first_data_bit": first_data_bit,
                "bits": bits,
                "decision_values": decision_values,
                "signed_soft": signed_soft,
                "distances": distances,
                "header_errors": header_errors,
            })

    records.sort(key=lambda item: (
        int(item["header_errors"]),
        -int(item["hypothesis"].preamble_score),
        -int(item["hypothesis"].flag_score),
    ))
    records = records[:max(1, maximum_observations)]
    if not records:
        return (
            None,
            "SOFT_REPAIR_NOT_ELIGIBLE",
            "No exact preamble route had enough symbols and a sufficiently close fixed header.",
        )

    observations: list[Any] = []
    phase_records: dict[int, dict[str, object]] = {}
    for phase, record in enumerate(records):
        first_data_bit = int(record["first_data_bit"])
        observation = advanced.PhaseObservation(
            phase=phase,
            offset=first_data_bit,
            raw_bits=record["bits"][first_data_bit:].copy(),
            soft_distance=record["distances"][first_data_bit:].copy(),
            header_errors=int(record["header_errors"]),
        )
        observations.append(observation)
        phase_records[phase] = record

    effort = str(repair_effort).strip().upper()
    if effort == "MAXIMUM":
        repair_limits = {
            "beam_width": 1024,
            "beam_branch_pool_size": 160,
            "maximum_structural_flips": 6,
            "maximum_beam_candidates_per_phase": 12,
            "maximum_syndrome_bits": 4,
            "syndrome_pool_size": 96,
        }
    elif effort == "STANDARD":
        repair_limits = {
            "beam_width": 512,
            "beam_branch_pool_size": 96,
            "maximum_structural_flips": 4,
            "maximum_beam_candidates_per_phase": 8,
            "maximum_syndrome_bits": 4,
            "syndrome_pool_size": 64,
        }
    else:
        raise ValueError("soft_repair_effort must be STANDARD or MAXIMUM")

    try:
        repair = advanced.advanced_repair_frame(
            observations,
            **repair_limits,
            minimum_arbitration_margin=0.05,
            # Existing GNU/MATLAB captures retain one closing flag, rather than
            # the three flags required by the decoder's conservative default.
            minimum_closing_flags=1,
        )
    except Exception as error:
        return None, "SOFT_REPAIR_FAILED", str(error)

    summary = repair.to_summary()
    status = str(summary.get("status", "SOFT_REPAIR_FAILED"))
    reason = str(summary.get("reason", ""))
    if (
        repair.frame is None
        or repair.selected_phase is None
        or status not in ("ORIGINAL_CRC_PASS", "AUTO_REPAIRED_UNIQUE")
    ):
        return None, status, reason

    record = phase_records[int(repair.selected_phase)]
    hypothesis = record["hypothesis"]
    opening = int(record["opening"])
    first_data_bit = int(record["first_data_bit"])
    original_bits = record["bits"].copy()
    corrected_bits = original_bits.copy()
    global_flips = tuple(first_data_bit + int(index) for index in repair.flip_indices)
    if any(index < 0 or index >= corrected_bits.size for index in global_flips):
        return None, "SOFT_REPAIR_FAILED", "The decoder returned an out-of-range bit index."
    if global_flips:
        corrected_bits[list(global_flips)] ^= 1
    closing = first_data_bit + int(repair.selected_flag_start)
    available_closing = original_bits[closing:min(closing + 8, original_bits.size)]
    closing_errors = 8 - int(available_closing.size)
    if available_closing.size:
        closing_errors += int(np.count_nonzero(
            available_closing != FLAG_BITS[:available_closing.size]
        ))

    summary.update({
        "opening_flag_bit_index": opening,
        "closing_flag_bit_index": closing,
        "global_flip_indices": list(global_flips),
        "observation_header_bit_errors": int(record["header_errors"]),
        "carrier_hz": (
            hypothesis.reported_carrier_hz
            if hypothesis.reported_carrier_hz is not None
            else hypothesis.carrier_hz
        ),
        "symbol_rate": hypothesis.symbol_rate,
        "polarity": "NORMAL" if hypothesis.polarity == 1 else "INVERTED",
        "timing_recovery": hypothesis.timing_mode,
        "filter_method": hypothesis.filter_method,
        "soft_repair_effort": effort,
    })
    frame = FrameCandidate(
        opening_flag_index=opening,
        closing_flag_index=closing,
        frame=repair.frame,
        crc_ok=fcs_ok(repair.frame),
        boundary_mode=(
            "SOFT_REPAIRED" if status == "AUTO_REPAIRED_UNIQUE" else "SOFT_VALIDATED"
        ),
        header_bit_errors=_header_bit_errors(repair.frame),
        closing_flag_bit_errors=closing_errors,
        destuff_violations=0,
        soft_repair_summary=summary,
    )
    if not frame.crc_ok:
        return None, "SOFT_REPAIR_FAILED", "Repaired frame failed independent CRC validation."
    return (
        SoftRepairSelection(
            hypothesis=hypothesis,
            frame=frame,
            original_bits=original_bits,
            corrected_bits=corrected_bits,
            decision_values=record["decision_values"].copy(),
            signed_soft_values=record["signed_soft"].copy(),
            soft_distances=record["distances"].copy(),
            status=status,
            reason=reason,
        ),
        status,
        reason,
    )


def _build_timing_shortlist(
    filtered: np.ndarray,
    sample_rate: float,
    symbol_rates: Iterable[float],
    phase_count: int,
    scan_min_hz: int,
    scan_max_hz: int,
    scan_step_hz: int,
    top_thresholds_per_phase: int,
    test_inverted_polarity: bool,
    filter_method: str = "RAW",
) -> list[TimingHypothesis]:
    source_index = np.arange(filtered.size, dtype=np.float64)
    centers = np.arange(scan_min_hz, scan_max_hz + 1, scan_step_hz, dtype=np.int32)
    shortlist: list[TimingHypothesis] = []
    polarities = (1, -1) if test_inverted_polarity else (1,)
    for candidate_symbol_rate in symbol_rates:
        samples_per_symbol = sample_rate / float(candidate_symbol_rate)
        for phase_index in range(phase_count):
            phase_fraction = phase_index / phase_count
            start = phase_fraction * samples_per_symbol
            positions = start + np.arange(
                int(max(0, (filtered.size - 1 - start) / samples_per_symbol)) + 1,
                dtype=np.float64,
            ) * samples_per_symbol
            positions = positions[positions < filtered.size - 1]
            if positions.size < 16:
                continue
            values = np.interp(positions, source_index, filtered)
            for polarity in polarities:
                flag_pattern = FLAG_BITS if polarity == 1 else 1 - FLAG_BITS
                preamble_pattern = (
                    PREAMBLE_BITS if polarity == 1 else 1 - PREAMBLE_BITS
                )
                flag_counts = _pattern_count_by_threshold(
                    values,
                    flag_pattern,
                    scan_min_hz,
                    scan_max_hz,
                    scan_step_hz,
                )
                preamble_counts = _pattern_count_by_threshold(
                    values,
                    preamble_pattern,
                    scan_min_hz,
                    scan_max_hz,
                    scan_step_hz,
                )
                rank = preamble_counts.astype(np.int64) * 1_000_000 + flag_counts
                take = min(max(1, top_thresholds_per_phase), rank.size)
                candidate_indices = np.argpartition(rank, -take)[-take:]
                for candidate_index in candidate_indices:
                    shortlist.append(
                        TimingHypothesis(
                            preamble_score=int(preamble_counts[candidate_index]),
                            flag_score=int(flag_counts[candidate_index]),
                            phase_fraction=phase_fraction,
                            phase_index=phase_index,
                            carrier_hz=int(centers[candidate_index]),
                            symbol_rate=float(candidate_symbol_rate),
                            polarity=polarity,
                            values=values,
                            filter_method=filter_method,
                        )
                    )
    return shortlist


def _evaluate_hypothesis(
    hypothesis: TimingHypothesis,
    enable_boundary_recovery: bool,
    recovery_max_header_bit_errors: int,
    recovery_max_closing_flag_bit_errors: int,
    recovery_max_destuff_violations: int,
    allow_diagnostic: bool,
) -> HypothesisResult:
    bits, _, _, _ = _hard_and_soft_decisions(hypothesis)
    normal_bits = (hypothesis.values > hypothesis.carrier_hz).astype(np.uint8)
    frames = frame_candidates(bits)
    complete_frames = [
        frame for frame in frames if len(frame.frame) == SCIONX_FRAME_BYTES
    ]
    valid_frames = [frame for frame in complete_frames if frame.crc_ok]
    recovered_frames: list[FrameCandidate] = []
    if (
        enable_boundary_recovery
        and not complete_frames
        and hypothesis.preamble_score >= 1
    ):
        recovered_frames = recover_frame_boundaries(
            bits,
            max_header_bit_errors=recovery_max_header_bit_errors,
            max_closing_flag_bit_errors=recovery_max_closing_flag_bit_errors,
            max_destuff_violations=recovery_max_destuff_violations,
        )
    recovered_valid = [frame for frame in recovered_frames if frame.crc_ok]
    best_recovery = min(
        recovered_frames,
        key=lambda frame: (
            frame.header_bit_errors,
            frame.closing_flag_bit_errors,
            frame.destuff_violations,
        ),
        default=None,
    )
    diagnostic_frames: list[FrameCandidate] = []
    if (
        allow_diagnostic
        and not complete_frames
        and not recovered_frames
        and hypothesis.preamble_score >= 1
    ):
        candidates = diagnostic_frame_boundaries(bits)
        if candidates:
            diagnostic_frames = [
                min(
                    candidates,
                    key=lambda frame: (
                        frame.header_bit_errors,
                        frame.closing_flag_bit_errors,
                        frame.destuff_violations,
                    ),
                )
            ]
    best_diagnostic = diagnostic_frames[0] if diagnostic_frames else None
    longest = max((len(frame.frame) for frame in valid_frames), default=0)
    longest_candidate = max((len(frame.frame) for frame in frames), default=0)
    closest_size_error = min(
        (abs(len(frame.frame) - SCIONX_FRAME_BYTES) for frame in frames),
        default=SCIONX_FRAME_BYTES,
    )
    low_values = hypothesis.values[normal_bits == 0]
    high_values = hypothesis.values[normal_bits == 1]
    if low_values.size and high_values.size:
        low_tone = float(np.median(low_values))
        high_tone = float(np.median(high_values))
        tone_midpoint = 0.5 * (low_tone + high_tone)
        midpoint_error = abs(hypothesis.carrier_hz - tone_midpoint)
        eye_margin = min(
            hypothesis.carrier_hz - low_tone,
            high_tone - hypothesis.carrier_hz,
        )
    else:
        midpoint_error = float("inf")
        eye_margin = float("-inf")
    key = (
        len(valid_frames),
        len(recovered_valid),
        longest,
        len(complete_frames),
        len(recovered_frames),
        -(best_recovery.header_bit_errors if best_recovery else 10_000),
        -(best_recovery.closing_flag_bit_errors if best_recovery else 10_000),
        -(best_recovery.destuff_violations if best_recovery else 10_000),
        len(diagnostic_frames),
        -(best_diagnostic.header_bit_errors if best_diagnostic else 10_000),
        -(
            best_diagnostic.closing_flag_bit_errors
            if best_diagnostic
            else 10_000
        ),
        -(best_diagnostic.destuff_violations if best_diagnostic else 10_000),
        -closest_size_error,
        longest_candidate,
        eye_margin,
        -midpoint_error,
        hypothesis.preamble_score,
        hypothesis.flag_score,
    )
    return HypothesisResult(
        hypothesis=hypothesis,
        frames=frames,
        recovered_frames=recovered_frames,
        diagnostic_frames=diagnostic_frames,
        bits=bits,
        key=key,
    )


def _best_hypothesis(
    shortlist: Iterable[TimingHypothesis],
    enable_boundary_recovery: bool,
    recovery_max_header_bit_errors: int,
    recovery_max_closing_flag_bit_errors: int,
    recovery_max_destuff_violations: int,
    allow_diagnostic: bool,
) -> HypothesisResult | None:
    best: HypothesisResult | None = None
    for hypothesis in shortlist:
        result = _evaluate_hypothesis(
            hypothesis,
            enable_boundary_recovery,
            recovery_max_header_bit_errors,
            recovery_max_closing_flag_bit_errors,
            recovery_max_destuff_violations,
            allow_diagnostic,
        )
        if best is None or result.key > best.key:
            best = result
    return best


def _gardner_symbol_values(
    soft_signal: np.ndarray,
    start_sample: float,
    nominal_sps: float,
    loop_bandwidth: float,
    feedback_sign: float,
) -> np.ndarray:
    """Track symbol phase and clock rate with a bounded Gardner timing loop."""
    damping = 1.0
    theta = loop_bandwidth / (damping + 1.0 / (4.0 * damping))
    denominator = 1.0 + 2.0 * damping * theta + theta * theta
    alpha = 4.0 * damping * theta / denominator
    beta = 4.0 * theta * theta / denominator
    minimum_sps = nominal_sps * 0.98
    maximum_sps = nominal_sps * 1.02
    source_index = np.arange(soft_signal.size, dtype=np.float64)
    time = float(start_sample)
    omega = float(nominal_sps)
    values: list[float] = []
    previous = float(np.interp(max(0.0, time - omega), source_index, soft_signal))
    while time + 0.55 * omega < soft_signal.size - 1.0:
        current = float(np.interp(time, source_index, soft_signal))
        midpoint = float(np.interp(time - 0.5 * omega, source_index, soft_signal))
        error = float(np.clip(midpoint * (previous - current), -2.0, 2.0))
        omega = float(
            np.clip(
                omega + feedback_sign * beta * error,
                minimum_sps,
                maximum_sps,
            )
        )
        values.append(current)
        time += omega + feedback_sign * alpha * error
        previous = current
    return np.asarray(values, dtype=np.float64)


def _mm_symbol_values(
    soft_signal: np.ndarray,
    start_sample: float,
    nominal_sps: float,
    gain_mu: float,
    gain_omega: float,
    feedback_sign: float,
) -> np.ndarray:
    """Decision-directed Mueller and Muller symbol clock tracking."""
    source_index = np.arange(soft_signal.size, dtype=np.float64)
    time = float(start_sample)
    omega = float(nominal_sps)
    minimum_sps = nominal_sps * 0.98
    maximum_sps = nominal_sps * 1.02
    values: list[float] = []
    previous = 0.0
    previous_decision = 0.0
    first = True
    while time < soft_signal.size - 1.0:
        current = float(np.interp(time, source_index, soft_signal))
        decision = 1.0 if current >= 0.0 else -1.0
        error = 0.0
        if not first:
            error = previous_decision * current - decision * previous
            error = float(np.clip(error, -2.0, 2.0))
            omega = float(
                np.clip(
                    omega + feedback_sign * gain_omega * error,
                    minimum_sps,
                    maximum_sps,
                )
            )
        values.append(current)
        time += omega + feedback_sign * gain_mu * error
        previous = current
        previous_decision = decision
        first = False
    return np.asarray(values, dtype=np.float64)


def _clock_recovery_shortlist(
    filtered: np.ndarray,
    acquisition: HypothesisResult,
    sample_rate: float,
) -> list[TimingHypothesis]:
    hypothesis = acquisition.hypothesis
    centered = filtered - hypothesis.carrier_hz
    scale = float(np.percentile(np.abs(centered), 70.0))
    soft_signal = np.clip(centered / max(scale, 1.0), -3.0, 3.0)
    samples_per_symbol = sample_rate / hypothesis.symbol_rate
    base_start = hypothesis.phase_fraction * samples_per_symbol
    variants = [
        (fraction, bandwidth, feedback_sign)
        for bandwidth in (0.002, 0.008, 0.02)
        for feedback_sign in (-1.0, 1.0)
        for fraction in (-0.45, -0.20, 0.0, 0.20, 0.45)
    ]
    recovered: list[TimingHypothesis] = []
    for phase_index, (fraction, bandwidth, feedback_sign) in enumerate(variants):
        values = _gardner_symbol_values(
            soft_signal,
            max(0.0, base_start + fraction * samples_per_symbol),
            samples_per_symbol,
            bandwidth,
            feedback_sign,
        )
        if values.size < 16:
            continue
        normal_bits = (values > 0.0).astype(np.uint8)
        bits = normal_bits if hypothesis.polarity == 1 else 1 - normal_bits
        recovered.append(
            TimingHypothesis(
                preamble_score=_exact_pattern_count(bits, PREAMBLE_BITS),
                flag_score=_exact_pattern_count(bits, FLAG_BITS),
                phase_fraction=fraction,
                phase_index=phase_index,
                carrier_hz=0,
                symbol_rate=hypothesis.symbol_rate,
                polarity=hypothesis.polarity,
                values=values,
                timing_mode="GARDNER",
                reported_carrier_hz=hypothesis.carrier_hz,
                filter_method=hypothesis.filter_method,
            )
        )
    next_phase_index = len(recovered)
    mm_variants = [
        (fraction, gain_mu, gain_omega, feedback_sign)
        for gain_mu, gain_omega in ((0.02, 0.0002), (0.05, 0.0008), (0.10, 0.002))
        for feedback_sign in (-1.0, 1.0)
        for fraction in (-0.45, -0.20, 0.0, 0.20, 0.45)
    ]
    for variant_index, (fraction, gain_mu, gain_omega, feedback_sign) in enumerate(
        mm_variants
    ):
        values = _mm_symbol_values(
            soft_signal,
            max(0.0, base_start + fraction * samples_per_symbol),
            samples_per_symbol,
            gain_mu,
            gain_omega,
            feedback_sign,
        )
        if values.size < 16:
            continue
        normal_bits = (values > 0.0).astype(np.uint8)
        bits = normal_bits if hypothesis.polarity == 1 else 1 - normal_bits
        recovered.append(
            TimingHypothesis(
                preamble_score=_exact_pattern_count(bits, PREAMBLE_BITS),
                flag_score=_exact_pattern_count(bits, FLAG_BITS),
                phase_fraction=fraction,
                phase_index=next_phase_index + variant_index,
                carrier_hz=0,
                symbol_rate=hypothesis.symbol_rate,
                polarity=hypothesis.polarity,
                values=values,
                timing_mode="M&M",
                reported_carrier_hz=hypothesis.carrier_hz,
                filter_method=hypothesis.filter_method,
            )
        )
    return recovered


@dataclass
class ScanResult:
    carrier_hz: int
    symbol_rate: float
    polarity: str
    timing_recovery: str
    filter_method: str
    phase_index: int
    phase_fraction: float
    preamble_count: int
    flag_count: int
    frame_count: int
    complete_frame_count: int
    recovered_frame_count: int
    diagnostic_frame_count: int
    valid_frames: list[FrameCandidate]
    output_frames: list[FrameCandidate]
    raw_bits: np.ndarray
    raw_frame_bits: list[np.ndarray]
    raw_frame_original_bits: list[np.ndarray]
    raw_frame_decision_values: list[np.ndarray]
    raw_frame_signed_soft_values: list[np.ndarray]
    raw_frame_soft_distances: list[np.ndarray]
    raw_frame_symbol_starts: list[int]
    soft_repair_status: str
    soft_repair_reason: str


def _raw_span_for_result(bits: np.ndarray, frames: list[FrameCandidate]) -> np.ndarray:
    """Keep the preamble flags, selected frame and its closing flag."""
    flags = find_flag_positions(bits)
    valid = [frame for frame in frames if frame.crc_ok]
    if valid:
        chosen = max(valid, key=lambda frame: len(frame.frame))
        opening = chosen.opening_flag_index
        run_start = opening
        prior = flags[flags < opening]
        for position in prior[::-1]:
            if run_start - int(position) == 8:
                run_start = int(position)
            else:
                break
        return bits[run_start:chosen.closing_flag_index + 8].copy()
    if flags.size >= 2:
        return bits[int(flags[0]):int(flags[-1]) + 8].copy()
    return bits.copy()


def _raw_bounds_for_frame(bits: np.ndarray, frame: FrameCandidate) -> tuple[int, int]:
    """Return bounds for one frame with contiguous preamble and closing flag."""
    flags = find_flag_positions(bits)
    opening = frame.opening_flag_index
    run_start = opening
    prior = flags[flags < opening]
    for position in prior[::-1]:
        if run_start - int(position) == 8:
            run_start = int(position)
        else:
            break
    return run_start, frame.closing_flag_index + 8


def _raw_span_for_frame(bits: np.ndarray, frame: FrameCandidate) -> np.ndarray:
    """Return one complete frame with its contiguous preamble and closing flag."""
    first, last = _raw_bounds_for_frame(bits, frame)
    return bits[first:last].copy()


def scan_burst(
    iq: np.ndarray,
    sample_rate: float = 100000.0,
    symbol_rate: float = 9600.0,
    deviation_hz: float = 2500.0,
    scan_min_hz: int = -20000,
    scan_max_hz: int = 20000,
    scan_step_hz: int = 1,
    phase_count: int = 32,
    lowpass_hz: float = 6500.0,
    top_thresholds_per_phase: int = 100,
    enable_boundary_recovery: bool = True,
    recovery_max_header_bit_errors: int = 16,
    recovery_max_closing_flag_bit_errors: int = 2,
    recovery_max_destuff_violations: int = 8,
    enable_flags_only_diagnostic: bool = True,
    test_inverted_polarity: bool = True,
    enable_symbol_rate_search: bool = True,
    symbol_rate_min: float = 9550.0,
    symbol_rate_max: float = 9650.0,
    symbol_rate_step: float = 1.0,
    baud_search_phase_count: int = 4,
    baud_top_thresholds_per_phase: int = 1,
    enable_clock_recovery: bool = True,
    enable_fir_path: bool = True,
    enable_soft_repair: bool = True,
    soft_repair_library: str | None = None,
    soft_repair_max_header_bit_errors: int = 24,
    soft_repair_max_observations: int = 32,
    soft_repair_effort: str = "MAXIMUM",
) -> ScanResult:
    """Select carrier, polarity, baud and tracked symbol timing for one burst."""
    del deviation_hz  # Kept as an explicit interface/documentation parameter.
    iq = np.asarray(iq, dtype=np.complex64).ravel()
    if iq.size < 32:
        raise ValueError("Burst is too short for GFSK/HDLC acquisition")
    discriminator_hz = np.angle(iq[1:] * np.conj(iq[:-1])).astype(np.float64)
    discriminator_hz *= sample_rate / (2.0 * np.pi)

    taps = firwin(81, lowpass_hz, fs=sample_rate, window=("kaiser", 6.76))
    raw_filtered = fftconvolve(discriminator_hz, taps, mode="same")
    filter_routes: list[tuple[str, np.ndarray]] = [("RAW", raw_filtered)]
    if enable_fir_path:
        smoothing_taps = np.asarray([0.25, 0.50, 0.25], dtype=np.float64)
        filter_routes.append((
            "FIR_025_050_025",
            np.convolve(raw_filtered, smoothing_taps, mode="same"),
        ))

    nominal_by_filter: dict[str, list[TimingHypothesis]] = {}
    nominal_shortlist: list[TimingHypothesis] = []
    for filter_method, filtered in filter_routes:
        route_shortlist = _build_timing_shortlist(
            filtered,
            sample_rate,
            [symbol_rate],
            phase_count,
            scan_min_hz,
            scan_max_hz,
            scan_step_hz,
            top_thresholds_per_phase,
            test_inverted_polarity,
            filter_method,
        )
        nominal_by_filter[filter_method] = route_shortlist
        nominal_shortlist.extend(route_shortlist)
    if not nominal_shortlist:
        raise ValueError("No carrier/timing candidates were generated")
    repair_hypotheses = list(nominal_shortlist)
    best = _best_hypothesis(
        nominal_shortlist,
        enable_boundary_recovery,
        recovery_max_header_bit_errors,
        recovery_max_closing_flag_bit_errors,
        recovery_max_destuff_violations,
        enable_flags_only_diagnostic,
    )
    assert best is not None

    trusted_crc_valid = any(frame.crc_ok for frame in best.trusted_frames)
    if enable_symbol_rate_search and not trusted_crc_valid:
        if symbol_rate_step <= 0 or symbol_rate_min > symbol_rate_max:
            raise ValueError("Invalid symbol-rate search range")
        baud_values = np.arange(
            symbol_rate_min,
            symbol_rate_max + 0.5 * symbol_rate_step,
            symbol_rate_step,
            dtype=np.float64,
        )
        expanded_shortlist: list[TimingHypothesis] = []
        for filter_method, filtered in filter_routes:
            route_shortlist = _build_timing_shortlist(
                filtered,
                sample_rate,
                baud_values,
                baud_search_phase_count,
                scan_min_hz,
                scan_max_hz,
                scan_step_hz,
                baud_top_thresholds_per_phase,
                test_inverted_polarity,
                filter_method,
            )
            nominal_by_filter[filter_method].extend(route_shortlist)
            expanded_shortlist.extend(route_shortlist)
        repair_hypotheses.extend(expanded_shortlist)
        expanded_best = _best_hypothesis(
            expanded_shortlist,
            enable_boundary_recovery,
            recovery_max_header_bit_errors,
            recovery_max_closing_flag_bit_errors,
            recovery_max_destuff_violations,
            enable_flags_only_diagnostic,
        )
        if expanded_best is not None and expanded_best.key > best.key:
            best = expanded_best

    trusted_crc_valid = any(frame.crc_ok for frame in best.trusted_frames)
    if enable_clock_recovery and not trusted_crc_valid:
        clock_shortlist: list[TimingHypothesis] = []
        filter_signal = dict(filter_routes)
        for filter_method, route_hypotheses in nominal_by_filter.items():
            route_best = _best_hypothesis(
                route_hypotheses,
                enable_boundary_recovery,
                recovery_max_header_bit_errors,
                recovery_max_closing_flag_bit_errors,
                recovery_max_destuff_violations,
                enable_flags_only_diagnostic,
            )
            if route_best is not None:
                clock_shortlist.extend(_clock_recovery_shortlist(
                    filter_signal[filter_method], route_best, sample_rate
                ))
        repair_hypotheses.extend(clock_shortlist)
        clock_best = _best_hypothesis(
            clock_shortlist,
            enable_boundary_recovery,
            recovery_max_header_bit_errors,
            recovery_max_closing_flag_bit_errors,
            recovery_max_destuff_violations,
            enable_flags_only_diagnostic,
        )
        if clock_best is not None and clock_best.key > best.key:
            best = clock_best

    hypothesis = best.hypothesis
    frames = best.frames
    recovered_frames = best.recovered_frames
    diagnostic_frames = best.diagnostic_frames
    bits = best.bits.copy()
    original_bits = bits.copy()
    (
        _,
        decision_values,
        signed_soft_values,
        soft_distances,
    ) = _hard_and_soft_decisions(hypothesis)
    strict_complete = [frame for frame in frames if len(frame.frame) == SCIONX_FRAME_BYTES]
    if strict_complete:
        output_frames = strict_complete
    elif recovered_frames:
        output_frames = recovered_frames
    else:
        output_frames = diagnostic_frames
    valid_frames = [frame for frame in output_frames if frame.crc_ok]
    soft_repair_status = "NOT_NEEDED" if valid_frames else "NOT_RUN"
    soft_repair_reason = "An original/recovered frame already passed CRC." if valid_frames else ""
    if enable_soft_repair and not valid_frames:
        repair_hypotheses.append(best.hypothesis)
        selection, soft_repair_status, soft_repair_reason = _attempt_advanced_soft_repair(
            repair_hypotheses,
            library_root=soft_repair_library,
            maximum_header_bit_errors=soft_repair_max_header_bit_errors,
            maximum_observations=soft_repair_max_observations,
            repair_effort=soft_repair_effort,
        )
        if selection is not None:
            hypothesis = selection.hypothesis
            frames = [selection.frame]
            recovered_frames = []
            diagnostic_frames = []
            output_frames = [selection.frame]
            valid_frames = [selection.frame]
            original_bits = selection.original_bits
            bits = selection.corrected_bits
            decision_values = selection.decision_values
            signed_soft_values = selection.signed_soft_values
            soft_distances = selection.soft_distances

    frame_bounds = [_raw_bounds_for_frame(bits, frame) for frame in output_frames]
    return ScanResult(
        carrier_hz=(
            hypothesis.reported_carrier_hz
            if hypothesis.reported_carrier_hz is not None
            else hypothesis.carrier_hz
        ),
        symbol_rate=hypothesis.symbol_rate,
        polarity="NORMAL" if hypothesis.polarity == 1 else "INVERTED",
        timing_recovery=hypothesis.timing_mode,
        filter_method=hypothesis.filter_method,
        phase_index=hypothesis.phase_index,
        phase_fraction=hypothesis.phase_fraction,
        preamble_count=hypothesis.preamble_score,
        flag_count=hypothesis.flag_score,
        frame_count=len(frames),
        complete_frame_count=len(output_frames),
        recovered_frame_count=sum(
            frame.boundary_mode == "RECOVERED" for frame in output_frames
        ),
        diagnostic_frame_count=sum(
            frame.boundary_mode == "FLAGS_ONLY_DIAGNOSTIC"
            for frame in output_frames
        ),
        valid_frames=valid_frames,
        output_frames=output_frames,
        raw_bits=_raw_span_for_result(bits, frames),
        raw_frame_bits=[bits[first:last].copy() for first, last in frame_bounds],
        raw_frame_original_bits=[
            original_bits[first:last].copy() for first, last in frame_bounds
        ],
        raw_frame_decision_values=[
            decision_values[first:last].copy() for first, last in frame_bounds
        ],
        raw_frame_signed_soft_values=[
            signed_soft_values[first:last].copy() for first, last in frame_bounds
        ],
        raw_frame_soft_distances=[
            soft_distances[first:last].copy() for first, last in frame_bounds
        ],
        raw_frame_symbol_starts=[first for first, _ in frame_bounds],
        soft_repair_status=soft_repair_status,
        soft_repair_reason=soft_repair_reason,
    )


class StreamingBurstScanner:
    """Streaming power detector for complex IQ or FM-discriminator audio."""

    REPORT_FIELDS = [
        "burst_index",
        "start_seconds",
        "duration_seconds",
        "power_db",
        "noise_floor_db",
        "detector_method",
        "detector_score",
        "detector_energy_z",
        "best_carrier_hz",
        "best_symbol_rate",
        "polarity",
        "timing_recovery",
        "filter_method",
        "extended_search_used",
        "phase_index",
        "phase_fraction",
        "preamble_count",
        "flag_count",
        "frame_candidate_count",
        "complete_frame_count",
        "recovered_frame_count",
        "diagnostic_frame_count",
        "soft_repaired_count",
        "soft_repair_bit_flip_count",
        "crc_valid_count",
        "crc_failed_count",
        "raw_bit_count",
        "raw_hex_byte_count",
        "output_frame_count",
        "boundary_modes",
        "output_crc_statuses",
        "gnu_hex_files",
        "soft_symbol_files",
        "soft_repair_files",
        "soft_repair_status",
        "soft_repair_reason",
        "status",
    ]

    def __init__(
        self,
        output_txt: str,
        report_csv: str,
        sample_rate: float = 100000.0,
        symbol_rate: float = 9600.0,
        deviation_hz: float = 2500.0,
        scan_min_hz: int = -20000,
        scan_max_hz: int = 20000,
        scan_step_hz: int = 1,
        phase_count: int = 32,
        detection_margin_db: float = 4.0,
        block_size: int = 1000,
        pretrigger_ms: float = 100.0,
        end_hold_ms: float = 350.0,
        minimum_burst_ms: float = 80.0,
        maximum_burst_ms: float = 5000.0,
        enable_boundary_recovery: bool = True,
        recovery_max_header_bit_errors: int = 16,
        recovery_max_closing_flag_bit_errors: int = 2,
        recovery_max_destuff_violations: int = 8,
        enable_flags_only_diagnostic: bool = True,
        test_inverted_polarity: bool = True,
        enable_symbol_rate_search: bool = True,
        symbol_rate_min: float = 9550.0,
        symbol_rate_max: float = 9650.0,
        symbol_rate_step: float = 1.0,
        baud_search_phase_count: int = 4,
        enable_clock_recovery: bool = True,
        enable_fir_path: bool = True,
        extended_search_min_power_excess_db: float = 2.0,
        no_preamble_full_search_top_n: int = 3,
        enable_soft_repair: bool = True,
        soft_repair_library: str | None = None,
        soft_repair_max_header_bit_errors: int = 24,
        soft_repair_max_observations: int = 32,
        soft_repair_effort: str = "MAXIMUM",
        input_mode: str = "IQ",
        audio_hz_per_unit: float = 0.0,
        emit_decoded_candidate_files: bool = False,
    ) -> None:
        self.output_txt = Path(output_txt)
        self.report_csv = Path(report_csv)
        self.sample_rate = float(sample_rate)
        self.symbol_rate = float(symbol_rate)
        self.deviation_hz = float(deviation_hz)
        self.scan_min_hz = int(scan_min_hz)
        self.scan_max_hz = int(scan_max_hz)
        self.scan_step_hz = int(scan_step_hz)
        self.phase_count = int(phase_count)
        self.enable_boundary_recovery = bool(enable_boundary_recovery)
        self.recovery_max_header_bit_errors = int(recovery_max_header_bit_errors)
        self.recovery_max_closing_flag_bit_errors = int(
            recovery_max_closing_flag_bit_errors
        )
        self.recovery_max_destuff_violations = int(
            recovery_max_destuff_violations
        )
        self.enable_flags_only_diagnostic = bool(enable_flags_only_diagnostic)
        self.test_inverted_polarity = bool(test_inverted_polarity)
        self.enable_symbol_rate_search = bool(enable_symbol_rate_search)
        self.symbol_rate_min = float(symbol_rate_min)
        self.symbol_rate_max = float(symbol_rate_max)
        self.symbol_rate_step = float(symbol_rate_step)
        self.baud_search_phase_count = int(baud_search_phase_count)
        self.enable_clock_recovery = bool(enable_clock_recovery)
        self.enable_fir_path = bool(enable_fir_path)
        self.no_preamble_full_search_top_n = int(no_preamble_full_search_top_n)
        self.enable_soft_repair = bool(enable_soft_repair)
        self.soft_repair_library = soft_repair_library
        self.soft_repair_max_header_bit_errors = int(
            soft_repair_max_header_bit_errors
        )
        self.soft_repair_max_observations = int(soft_repair_max_observations)
        self.soft_repair_effort = str(soft_repair_effort).strip().upper()
        self.input_mode = str(input_mode).strip().upper()
        self.audio_hz_per_unit = float(audio_hz_per_unit)
        self.emit_decoded_candidate_files = bool(emit_decoded_candidate_files)
        self.extended_search_min_power_excess_db = float(
            extended_search_min_power_excess_db
        )
        self.detection_margin_db = float(detection_margin_db)
        self.block_size = int(block_size)
        self.pretrigger_blocks = max(1, int(round(pretrigger_ms * 1e-3 * self.sample_rate / self.block_size)))
        self.end_hold_blocks = max(1, int(round(end_hold_ms * 1e-3 * self.sample_rate / self.block_size)))
        self.minimum_samples = int(round(minimum_burst_ms * 1e-3 * self.sample_rate))
        self.maximum_samples = int(round(maximum_burst_ms * 1e-3 * self.sample_rate))

        if self.input_mode not in ("IQ", "AUDIO_DISCRIMINATOR"):
            raise ValueError("input_mode must be IQ or AUDIO_DISCRIMINATOR")
        if self.audio_hz_per_unit < 0:
            raise ValueError("audio_hz_per_unit must be zero (auto) or positive")
        if self.scan_step_hz <= 0 or self.scan_min_hz > self.scan_max_hz:
            raise ValueError("Invalid carrier scan range")
        if (
            self.symbol_rate_step <= 0
            or self.symbol_rate_min > self.symbol_rate_max
            or self.baud_search_phase_count <= 0
            or self.soft_repair_max_header_bit_errors < 0
            or self.soft_repair_max_observations <= 0
            or self.no_preamble_full_search_top_n < 0
        ):
            raise ValueError("Invalid symbol-rate search configuration")
        if self.soft_repair_effort not in ("STANDARD", "MAXIMUM"):
            raise ValueError("soft_repair_effort must be STANDARD or MAXIMUM")
        # GRC instantiates Embedded Python blocks while inspecting their ports.
        # Do not truncate user output merely because a flowgraph is opened or
        # compiled; initialize files only when real samples first arrive.
        self._started = False
        self._remainder = np.empty(0, dtype=np.complex64)
        self._pretrigger: deque[np.ndarray] = deque(maxlen=self.pretrigger_blocks)
        self._noise_history: deque[float] = deque(maxlen=1000)
        self._active_blocks: list[np.ndarray] = []
        self._active = False
        self._below_count = 0
        self._block_index = 0
        self._burst_index = 0
        self._output_frame_index = 0
        self._burst_start_sample = 0
        self._explicit_noise_floor_db = float("nan")
        self._force_extended_search_once = False
        self._detector_method = "POWER_RMS"
        self._detector_score = float("nan")
        self._detector_energy_z = float("nan")
        # Quick-scan failures retain only the strongest N raw bursts.  They are
        # reprocessed at EOF with polarity, baud and clock-recovery searches.
        self._no_preamble_candidates: list[dict[str, Any]] = []
        self._force_full_scan_once = False

    def feed(self, samples: np.ndarray) -> None:
        data = np.asarray(samples, dtype=np.complex64).ravel()
        if data.size and not self._started:
            self._start_outputs()
        if self._remainder.size:
            data = np.concatenate((self._remainder, data))
        count = data.size // self.block_size
        for index in range(count):
            block = data[index * self.block_size:(index + 1) * self.block_size].copy()
            self._process_power_block(block)
        self._remainder = data[count * self.block_size:].copy()

    def finish(self) -> None:
        if not self._started:
            return
        if self._remainder.size:
            padded = np.pad(self._remainder, (0, self.block_size - self._remainder.size))
            self._process_power_block(padded.astype(np.complex64))
            self._remainder = np.empty(0, dtype=np.complex64)
        if self._active and self._active_blocks:
            self._finalize_burst()
        self._run_no_preamble_full_searches()

    def process_explicit_burst(
        self,
        samples: np.ndarray,
        start_sample: int,
        *,
        noise_floor_db: float = float("nan"),
        force_extended_search: bool = True,
        force_full_search: bool = False,
        detector_method: str = "EXPLICIT",
        detector_score: float = float("nan"),
        detector_energy_z: float = float("nan"),
    ) -> None:
        """Decode one externally detected burst through the normal scan/output path.

        This entry point lets specialised front ends (notably the OGG matched-
        correlation detector) select burst windows without duplicating GFSK,
        clock-recovery, soft-repair, or raw-HDLC output handling.
        """
        burst = np.asarray(samples, dtype=np.complex64).ravel()
        if burst.size == 0:
            return
        if not self._started:
            self._start_outputs()
        self._burst_start_sample = max(0, int(start_sample))
        self._explicit_noise_floor_db = float(noise_floor_db)
        self._force_extended_search_once = bool(force_extended_search)
        self._force_full_scan_once = bool(force_full_search)
        self._detector_method = str(detector_method)
        self._detector_score = float(detector_score)
        self._detector_energy_z = float(detector_energy_z)
        self._active_blocks = [burst.copy()]
        self._active = True
        self._below_count = 0
        self._finalize_burst()

    def _start_outputs(self) -> None:
        self.output_txt.parent.mkdir(parents=True, exist_ok=True)
        self.report_csv.parent.mkdir(parents=True, exist_ok=True)
        # Remove only files generated by this demo, so a new run cannot pick up
        # stale Frame_NNN results from an older recording.
        for old_dir in self.output_txt.parent.glob("Frame_[0-9][0-9][0-9]"):
            if not old_dir.is_dir():
                continue
            for generated_name in (
                self.output_txt.name,
                "CRC_Status.txt",
                "Boundary_Recovery_Status.txt",
                "Recovered_Frame_HEX.txt",
                "Soft_Repaired_Frame_HEX.txt",
                "Diagnostic_Frame_HEX.txt",
                "Soft_Symbols.csv",
                "Soft_Repair_Status.json",
                "Raw_bit.txt",
                "TLM.txt",
            ):
                generated = old_dir / generated_name
                if generated.is_file():
                    generated.unlink()
            try:
                old_dir.rmdir()
            except OSError:
                pass
        with self.report_csv.open("w", newline="", encoding="utf-8-sig") as handle:
            csv.DictWriter(handle, fieldnames=self.REPORT_FIELDS).writeheader()
        self._started = True

    def _process_power_block(self, block: np.ndarray) -> None:
        power_db = 10.0 * np.log10(float(np.mean(np.abs(block) ** 2)) + 1e-20)
        noise_floor = float(np.median(self._noise_history)) if self._noise_history else power_db
        is_high = len(self._noise_history) >= 30 and power_db > noise_floor + self.detection_margin_db

        if not self._active:
            if is_high:
                self._active = True
                self._below_count = 0
                self._burst_start_sample = max(
                    0,
                    (self._block_index - len(self._pretrigger)) * self.block_size,
                )
                self._active_blocks = list(self._pretrigger)
                self._active_blocks.append(block)
            else:
                self._noise_history.append(power_db)
                self._pretrigger.append(block)
        else:
            self._active_blocks.append(block)
            if power_db <= noise_floor + self.detection_margin_db:
                self._below_count += 1
            else:
                self._below_count = 0
            active_samples = len(self._active_blocks) * self.block_size
            if self._below_count >= self.end_hold_blocks or active_samples >= self.maximum_samples:
                self._finalize_burst()
                self._pretrigger.clear()

        self._block_index += 1

    def _finalize_burst(self) -> None:
        burst = np.concatenate(self._active_blocks) if self._active_blocks else np.empty(0, np.complex64)
        self._active = False
        self._active_blocks = []
        self._below_count = 0
        explicit_noise_floor_db = self._explicit_noise_floor_db
        force_extended_search = self._force_extended_search_once
        detector_method = self._detector_method
        detector_score = self._detector_score
        detector_energy_z = self._detector_energy_z
        force_full_scan = self._force_full_scan_once
        self._explicit_noise_floor_db = float("nan")
        self._force_extended_search_once = False
        self._detector_method = "POWER_RMS"
        self._detector_score = float("nan")
        self._detector_energy_z = float("nan")
        self._force_full_scan_once = False
        if burst.size < self.minimum_samples:
            return
        self._burst_index += 1
        burst_power_db = 10.0 * np.log10(float(np.mean(np.abs(burst) ** 2)) + 1e-20)
        noise_floor_db = (
            explicit_noise_floor_db
            if np.isfinite(explicit_noise_floor_db)
            else (
                float(np.median(self._noise_history))
                if self._noise_history
                else float("nan")
            )
        )
        row = {
            "burst_index": self._burst_index,
            "start_seconds": self._burst_start_sample / self.sample_rate,
            "duration_seconds": burst.size / self.sample_rate,
            "power_db": burst_power_db,
            "noise_floor_db": noise_floor_db,
            "detector_method": detector_method,
            "detector_score": detector_score if np.isfinite(detector_score) else "",
            "detector_energy_z": (
                detector_energy_z if np.isfinite(detector_energy_z) else ""
            ),
            "best_carrier_hz": "",
            "best_symbol_rate": "",
            "polarity": "",
            "timing_recovery": "",
            "filter_method": "",
            "extended_search_used": "",
            "phase_index": "",
            "phase_fraction": "",
            "preamble_count": 0,
            "flag_count": 0,
            "frame_candidate_count": 0,
            "complete_frame_count": 0,
            "recovered_frame_count": 0,
            "diagnostic_frame_count": 0,
            "soft_repaired_count": 0,
            "soft_repair_bit_flip_count": 0,
            "crc_valid_count": 0,
            "crc_failed_count": 0,
            "raw_bit_count": 0,
            "raw_hex_byte_count": 0,
            "output_frame_count": 0,
            "boundary_modes": "",
            "output_crc_statuses": "",
            "gnu_hex_files": "",
            "soft_symbol_files": "",
            "soft_repair_files": "",
            "soft_repair_status": "NOT_RUN",
            "soft_repair_reason": "",
            "status": "SCAN_FAILED",
        }
        try:
            power_excess_db = burst_power_db - noise_floor_db
            extended_search_used = bool(
                force_extended_search
                or (
                    np.isfinite(power_excess_db)
                    and power_excess_db >= self.extended_search_min_power_excess_db
                )
            )
            analysis_iq = burst
            if self.input_mode == "AUDIO_DISCRIMINATOR":
                audio = np.real(burst).astype(np.float64)
                audio -= float(np.median(audio))
                if self.audio_hz_per_unit > 0:
                    audio_scale = self.audio_hz_per_unit
                else:
                    robust_level = float(np.percentile(np.abs(audio), 90.0))
                    audio_scale = self.deviation_hz / max(robust_level, 1e-9)
                discriminator_hz = np.clip(
                    audio * audio_scale,
                    -0.45 * self.sample_rate,
                    0.45 * self.sample_rate,
                )
                phase = np.cumsum(
                    2.0 * np.pi * discriminator_hz / self.sample_rate,
                    dtype=np.float64,
                )
                analysis_iq = np.exp(1j * phase).astype(np.complex64)
            scan_options = {
                "sample_rate": self.sample_rate,
                "symbol_rate": self.symbol_rate,
                "deviation_hz": self.deviation_hz,
                "scan_min_hz": self.scan_min_hz,
                "scan_max_hz": self.scan_max_hz,
                "scan_step_hz": self.scan_step_hz,
                "phase_count": self.phase_count,
                "enable_boundary_recovery": self.enable_boundary_recovery,
                "recovery_max_header_bit_errors": self.recovery_max_header_bit_errors,
                "recovery_max_closing_flag_bit_errors": (
                    self.recovery_max_closing_flag_bit_errors
                ),
                "recovery_max_destuff_violations": (
                    self.recovery_max_destuff_violations
                ),
                "enable_flags_only_diagnostic": self.enable_flags_only_diagnostic,
                "test_inverted_polarity": self.test_inverted_polarity,
                "enable_symbol_rate_search": (
                    self.enable_symbol_rate_search and extended_search_used
                ),
                "symbol_rate_min": self.symbol_rate_min,
                "symbol_rate_max": self.symbol_rate_max,
                "symbol_rate_step": self.symbol_rate_step,
                "baud_search_phase_count": self.baud_search_phase_count,
                "enable_clock_recovery": (
                    self.enable_clock_recovery and extended_search_used
                ),
                "enable_fir_path": self.enable_fir_path,
                "enable_soft_repair": self.enable_soft_repair,
                "soft_repair_library": self.soft_repair_library,
                "soft_repair_max_header_bit_errors": (
                    self.soft_repair_max_header_bit_errors
                ),
                "soft_repair_max_observations": self.soft_repair_max_observations,
                "soft_repair_effort": self.soft_repair_effort,
            }
            quick_options = dict(scan_options)
            quick_options.update({
                "enable_boundary_recovery": False,
                "enable_flags_only_diagnostic": False,
                "enable_symbol_rate_search": False,
                "enable_clock_recovery": False,
                "enable_soft_repair": False,
            })
            result = scan_burst(analysis_iq, **quick_options)
            full_search_used = bool(force_full_scan or result.preamble_count >= 1)
            if full_search_used:
                full_options = dict(scan_options)
                # A structural preamble or strongest-N fallback is stronger
                # evidence than average-power gating.  Always enable the
                # configured baud and clock searches for these selected bursts.
                full_options["enable_symbol_rate_search"] = self.enable_symbol_rate_search
                full_options["enable_clock_recovery"] = self.enable_clock_recovery
                result = scan_burst(analysis_iq, **full_options)
            elif self.no_preamble_full_search_top_n:
                rank_score = (
                    power_excess_db if np.isfinite(power_excess_db) else burst_power_db
                )
                self._remember_no_preamble_candidate({
                    "rank_score": float(rank_score),
                    "burst": burst.copy(),
                    "burst_index": self._burst_index,
                    "start_sample": self._burst_start_sample,
                    "noise_floor_db": noise_floor_db,
                    "detector_method": detector_method,
                    "detector_score": detector_score,
                    "detector_energy_z": detector_energy_z,
                })
            if result.complete_frame_count:
                if any(
                    frame.boundary_mode == "SOFT_REPAIRED"
                    for frame in result.output_frames
                ):
                    status = "SOFT_REPAIRED_CRC_VALID"
                elif result.diagnostic_frame_count:
                    if any(frame.crc_ok for frame in result.output_frames):
                        status = "FLAGS_ONLY_DIAGNOSTIC_CRC_MATCH_UNTRUSTED"
                    else:
                        status = "FLAGS_ONLY_DIAGNOSTIC_CRC_FAIL"
                elif result.recovered_frame_count:
                    if len(result.valid_frames) == result.complete_frame_count:
                        status = "BOUNDARY_RECOVERED_CRC_VALID"
                    elif result.valid_frames:
                        status = "BOUNDARY_RECOVERED_CRC_MIXED"
                    else:
                        status = "BOUNDARY_RECOVERED_CRC_FAILED_CONTINUED"
                else:
                    if len(result.valid_frames) == result.complete_frame_count:
                        status = "CRC_VALID"
                    elif result.valid_frames:
                        status = "CRC_MIXED"
                    else:
                        status = "CRC_FAILED_CONTINUED"
            else:
                status = "FLAGS_ONLY" if result.preamble_count >= 1 else "NO_HDLC_PREAMBLE"
            written_files: list[str] = []
            soft_symbol_files: list[str] = []
            soft_repair_files: list[str] = []
            output_crc_statuses: list[str] = []
            boundary_modes: list[str] = []
            raw_bit_count = 0
            raw_hex_byte_count = 0
            soft_repair_bit_flip_count = 0
            frame_streams = zip(
                result.output_frames,
                result.raw_frame_bits,
                result.raw_frame_original_bits,
                result.raw_frame_decision_values,
                result.raw_frame_signed_soft_values,
                result.raw_frame_soft_distances,
                result.raw_frame_symbol_starts,
            )
            for (
                frame,
                raw_frame_bits,
                raw_original_bits,
                raw_decision_values,
                raw_signed_soft,
                raw_soft_distances,
                raw_symbol_start,
            ) in frame_streams:
                self._output_frame_index += 1
                frame_dir = self.output_txt.parent / f"Frame_{self._output_frame_index:03d}"
                frame_dir.mkdir(parents=True, exist_ok=True)
                frame_file = frame_dir / self.output_txt.name
                packed = pack_lsb_first(raw_frame_bits)
                line = " ".join(f"{value:02X}" for value in packed)
                frame_file.write_text(line + "\n", encoding="ascii", newline="\n")
                if (
                    self.emit_decoded_candidate_files
                    and frame.boundary_mode in ("RECOVERED", "SOFT_REPAIRED", "SOFT_VALIDATED")
                ):
                    recovered_file = frame_dir / "Recovered_Frame_HEX.txt"
                    recovered_file.write_text(
                        " ".join(f"{value:02X}" for value in frame.frame) + "\n",
                        encoding="ascii",
                        newline="\n",
                    )
                    if frame.boundary_mode == "SOFT_REPAIRED":
                        soft_repaired_file = frame_dir / "Soft_Repaired_Frame_HEX.txt"
                        soft_repaired_file.write_text(
                            " ".join(f"{value:02X}" for value in frame.frame) + "\n",
                            encoding="ascii",
                            newline="\n",
                        )
                elif (
                    self.emit_decoded_candidate_files
                    and frame.boundary_mode == "FLAGS_ONLY_DIAGNOSTIC"
                ):
                    diagnostic_file = frame_dir / "Diagnostic_Frame_HEX.txt"
                    diagnostic_file.write_text(
                        " ".join(f"{value:02X}" for value in frame.frame) + "\n",
                        encoding="ascii",
                        newline="\n",
                    )
                soft_symbol_file = frame_dir / "Soft_Symbols.csv"
                with soft_symbol_file.open("w", newline="", encoding="utf-8-sig") as handle:
                    fieldnames = [
                        "relative_symbol_index",
                        "burst_symbol_index",
                        "hard_bit_original",
                        "hard_bit_output",
                        "corrected",
                        "decision_value",
                        "signed_soft_value",
                        "soft_distance",
                    ]
                    writer = csv.DictWriter(handle, fieldnames=fieldnames)
                    writer.writeheader()
                    for relative_index in range(raw_frame_bits.size):
                        original_bit = int(raw_original_bits[relative_index])
                        output_bit = int(raw_frame_bits[relative_index])
                        writer.writerow({
                            "relative_symbol_index": relative_index,
                            "burst_symbol_index": int(raw_symbol_start) + relative_index,
                            "hard_bit_original": original_bit,
                            "hard_bit_output": output_bit,
                            "corrected": int(original_bit != output_bit),
                            "decision_value": "{:.12g}".format(
                                float(raw_decision_values[relative_index])
                            ),
                            "signed_soft_value": "{:.12g}".format(
                                float(raw_signed_soft[relative_index])
                            ),
                            "soft_distance": "{:.12g}".format(
                                float(raw_soft_distances[relative_index])
                            ),
                        })
                soft_symbol_files.append(str(soft_symbol_file))

                repair_summary = frame.soft_repair_summary
                repair_audit = _soft_repair_audit_metrics(repair_summary)
                soft_repair_bit_flip_count += int(repair_audit["bit_flip_count"])
                repair_status_document = {
                    "soft_repair_status": result.soft_repair_status,
                    "soft_repair_reason": result.soft_repair_reason,
                    "repair_evaluated": repair_summary is not None,
                    "repair_applied": repair_audit["applied"],
                    "repair_bit_flip_count": repair_audit["bit_flip_count"],
                    "repair_structural_flip_count": repair_audit[
                        "structural_flip_count"
                    ],
                    "repair_syndrome_flip_count": repair_audit[
                        "syndrome_flip_count"
                    ],
                    "burst_index": self._burst_index,
                    "frame_index": self._output_frame_index,
                    "filter_method": result.filter_method,
                    "boundary_mode": frame.boundary_mode,
                    "header_bit_errors": frame.header_bit_errors,
                    "closing_flag_bit_errors": frame.closing_flag_bit_errors,
                    "destuff_violations": frame.destuff_violations,
                    "crc_ok": frame.crc_ok,
                    "algorithm_result": repair_summary,
                }
                # Always emit an audit record.  Previously this file existed only
                # after a successful repair, which made NOT_ELIGIBLE/FAILED look
                # indistinguishable from a disabled soft-repair stage.
                repair_file = frame_dir / "Soft_Repair_Status.json"
                repair_file.write_text(
                    json.dumps(
                        repair_status_document,
                        ensure_ascii=False,
                        indent=2,
                    ) + "\n",
                    encoding="utf-8",
                    newline="\n",
                )
                soft_repair_files.append(str(repair_file))
                crc_status = "PASS" if frame.crc_ok else "FAIL"
                calculated_fcs = crc16_x25(frame.frame[:-2])
                received_fcs = frame.frame[-2] | (frame.frame[-1] << 8)
                (frame_dir / "CRC_Status.txt").write_text(
                    "CRC Check: {status}\n"
                    "Calculated CRC-16/X-25: 0x{calculated:04X}\n"
                    "Received FCS: 0x{received:04X}\n"
                    "Burst index: {burst}\n"
                    "Carrier offset: {carrier:+d} Hz\n"
                    "Symbol rate: {baud:.3f} sym/s\n"
                    "FSK polarity: {polarity}\n"
                    "Timing recovery: {timing}\n"
                    "Filter method: {filter_method}\n"
                    "Boundary mode: {mode}\n"
                    "Soft repair status: {soft_status}\n"
                    "Soft repair method: {soft_method}\n"
                    "Soft repair applied: {soft_applied}\n"
                    "Soft repair bit flips: {soft_flips}\n"
                    "Soft repair structural flips: {structural_flips}\n"
                    "Soft repair CRC syndrome flips: {syndrome_flips}\n"
                    "CRC PASS after soft repair: {pass_after_repair}\n".format(
                        status=crc_status,
                        calculated=calculated_fcs,
                        received=received_fcs,
                        burst=self._burst_index,
                        carrier=result.carrier_hz,
                        baud=result.symbol_rate,
                        polarity=result.polarity,
                        timing=result.timing_recovery,
                        filter_method=result.filter_method,
                        mode=frame.boundary_mode,
                        soft_status=result.soft_repair_status,
                        soft_method=repair_audit["method"],
                        soft_applied="YES" if repair_audit["applied"] else "NO",
                        soft_flips=repair_audit["bit_flip_count"],
                        structural_flips=repair_audit["structural_flip_count"],
                        syndrome_flips=repair_audit["syndrome_flip_count"],
                        pass_after_repair=(
                            "YES"
                            if repair_audit["applied"] and frame.crc_ok
                            else "NO"
                        ),
                    ),
                    encoding="ascii",
                    newline="\n",
                )
                validation = (
                    "FAIL - DIAGNOSTIC ONLY"
                    if frame.boundary_mode == "FLAGS_ONLY_DIAGNOSTIC"
                    else "PASS"
                )
                (frame_dir / "Boundary_Recovery_Status.txt").write_text(
                    "Boundary Mode: {mode}\n"
                    "Opening Flag Bit Index: {opening}\n"
                    "Expected Closing Flag Bit Index: {closing}\n"
                    "Header Bit Errors: {header}\n"
                    "Closing Flag Bit Errors: {closing_errors}\n"
                    "Destuff Violations: {violations}\n"
                    "Frame Validation: {validation}\n"
                    "Telemetry Use: {telemetry_use}\n"
                    "CRC Check: {crc}\n"
                    "Filter method: {filter_method}\n"
                    "Soft repair applied: {soft_applied}\n"
                    "Soft repair bit flips: {soft_flips}\n"
                    "Soft repair structural flips: {structural_flips}\n"
                    "Soft repair CRC syndrome flips: {syndrome_flips}\n"
                    "CRC PASS after soft repair: {pass_after_repair}\n".format(
                        mode=frame.boundary_mode,
                        opening=frame.opening_flag_index,
                        closing=frame.closing_flag_index,
                        header=frame.header_bit_errors,
                        closing_errors=frame.closing_flag_bit_errors,
                        violations=frame.destuff_violations,
                        validation=validation,
                        telemetry_use=(
                            "REFERENCE ONLY - NOT VALIDATED"
                            if frame.boundary_mode == "FLAGS_ONLY_DIAGNOSTIC"
                            else (
                                "CRC-VALIDATED SOFT REPAIR"
                                if frame.boundary_mode == "SOFT_REPAIRED"
                                else "VALIDATED BOUNDARY"
                            )
                        ),
                        crc=crc_status,
                        filter_method=result.filter_method,
                        soft_applied="YES" if repair_audit["applied"] else "NO",
                        soft_flips=repair_audit["bit_flip_count"],
                        structural_flips=repair_audit["structural_flip_count"],
                        syndrome_flips=repair_audit["syndrome_flip_count"],
                        pass_after_repair=(
                            "YES"
                            if repair_audit["applied"] and frame.crc_ok
                            else "NO"
                        ),
                    ),
                    encoding="ascii",
                    newline="\n",
                )
                written_files.append(str(frame_file))
                output_crc_statuses.append(crc_status)
                boundary_modes.append(frame.boundary_mode)
                raw_bit_count += int(raw_frame_bits.size)
                raw_hex_byte_count += len(packed)
            row.update({
                "best_carrier_hz": result.carrier_hz,
                "best_symbol_rate": result.symbol_rate,
                "polarity": result.polarity,
                "timing_recovery": result.timing_recovery,
                "filter_method": result.filter_method,
                "extended_search_used": full_search_used,
                "phase_index": result.phase_index,
                "phase_fraction": result.phase_fraction,
                "preamble_count": result.preamble_count,
                "flag_count": result.flag_count,
                "frame_candidate_count": result.frame_count,
                "complete_frame_count": result.complete_frame_count,
                "recovered_frame_count": result.recovered_frame_count,
                "diagnostic_frame_count": result.diagnostic_frame_count,
                "soft_repaired_count": sum(
                    frame.boundary_mode == "SOFT_REPAIRED"
                    for frame in result.output_frames
                ),
                "soft_repair_bit_flip_count": soft_repair_bit_flip_count,
                "crc_valid_count": len(result.valid_frames),
                "crc_failed_count": result.complete_frame_count - len(result.valid_frames),
                "raw_bit_count": raw_bit_count,
                "raw_hex_byte_count": raw_hex_byte_count,
                "output_frame_count": len(written_files),
                "boundary_modes": ";".join(boundary_modes),
                "output_crc_statuses": ";".join(output_crc_statuses),
                "gnu_hex_files": ";".join(written_files),
                "soft_symbol_files": ";".join(soft_symbol_files),
                "soft_repair_files": ";".join(soft_repair_files),
                "soft_repair_status": result.soft_repair_status,
                "soft_repair_reason": result.soft_repair_reason,
                "status": status,
            })
            print(
                "SCION-X burst {idx}: t={time:.3f}s carrier={carrier:+d}Hz "
                "baud={baud:.1f} polarity={polarity} timing={timing} "
                "filter={filter_method} "
                "phase={phase} preambles={preambles} flags={flags} "
                "complete={complete} recovered={recovered} diagnostic={diagnostic} "
                "crc_valid={valid} "
                "crc_failed={failed} "
                "soft_repair={soft_repair} status={status}".format(
                    idx=self._burst_index,
                    time=row["start_seconds"],
                    carrier=result.carrier_hz,
                    baud=result.symbol_rate,
                    polarity=result.polarity,
                    timing=result.timing_recovery,
                    filter_method=result.filter_method,
                    phase=result.phase_index,
                    preambles=result.preamble_count,
                    flags=result.flag_count,
                    complete=result.complete_frame_count,
                    recovered=result.recovered_frame_count,
                    diagnostic=result.diagnostic_frame_count,
                    valid=len(result.valid_frames),
                    failed=result.complete_frame_count - len(result.valid_frames),
                    soft_repair=result.soft_repair_status,
                    status=status,
                )
            )
        except Exception as error:  # Keep later bursts decodable if one scan fails.
            row["status"] = "SCAN_FAILED: {}".format(error)
            print("SCION-X burst {} scan failed: {}".format(self._burst_index, error))
        with self.report_csv.open("a", newline="", encoding="utf-8-sig") as handle:
            csv.DictWriter(handle, fieldnames=self.REPORT_FIELDS).writerow(row)

    def _remember_no_preamble_candidate(self, candidate: dict[str, Any]) -> None:
        """Keep only the strongest quick-scan failures for the EOF retry."""
        self._no_preamble_candidates.append(candidate)
        self._no_preamble_candidates.sort(
            key=lambda item: (float(item["rank_score"]), -int(item["burst_index"])),
            reverse=True,
        )
        del self._no_preamble_candidates[self.no_preamble_full_search_top_n:]

    def _run_no_preamble_full_searches(self) -> None:
        """Run the expensive search on the strongest preamble-free bursts."""
        if not self._no_preamble_candidates:
            return
        candidates = list(self._no_preamble_candidates)
        self._no_preamble_candidates.clear()
        final_burst_index = self._burst_index
        print(
            "SCION-X fallback: full polarity/9550..9650 sym/s/clock search "
            "for burst(s) {}".format(
                ", ".join(str(item["burst_index"]) for item in candidates)
            )
        )
        for candidate in candidates:
            self._burst_index = int(candidate["burst_index"]) - 1
            self._burst_start_sample = int(candidate["start_sample"])
            self._explicit_noise_floor_db = float(candidate["noise_floor_db"])
            self._force_extended_search_once = True
            self._force_full_scan_once = True
            self._detector_method = str(candidate["detector_method"])
            self._detector_score = float(candidate["detector_score"])
            self._detector_energy_z = float(candidate["detector_energy_z"])
            self._active_blocks = [candidate["burst"]]
            self._active = True
            self._below_count = 0
            self._finalize_burst()
        self._burst_index = final_burst_index
        self._deduplicate_report_rows()

    def _deduplicate_report_rows(self) -> None:
        """Replace quick rows with their later full-search rows and restore order."""
        if not self.report_csv.is_file():
            return
        with self.report_csv.open("r", newline="", encoding="utf-8-sig") as handle:
            rows = list(csv.DictReader(handle))
        latest: dict[int, dict[str, str]] = {}
        for row in rows:
            latest[int(row["burst_index"])] = row
        with self.report_csv.open("w", newline="", encoding="utf-8-sig") as handle:
            writer = csv.DictWriter(handle, fieldnames=self.REPORT_FIELDS)
            writer.writeheader()
            for burst_index in sorted(latest):
                writer.writerow(latest[burst_index])


def decode_cs16_file(
    input_file: str,
    output_txt: str,
    report_csv: str,
    **scanner_options: object,
) -> None:
    """Standalone path using exactly the same scanner as the GNU Radio block."""
    scanner = StreamingBurstScanner(output_txt, report_csv, **scanner_options)
    raw = np.memmap(input_file, dtype="<i2", mode="r")
    usable = raw.size - raw.size % 2
    complex_count = usable // 2
    chunk_samples = 1_000_000
    for first in range(0, complex_count, chunk_samples):
        part = np.asarray(raw[2 * first:2 * min(complex_count, first + chunk_samples)], dtype=np.float32)
        iq = (part[0::2] + 1j * part[1::2]).astype(np.complex64) / 32767.0
        scanner.feed(iq)
    scanner.finish()


def _argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_file")
    parser.add_argument("--output", default="GNU_GFSK_Raw_HEX.txt")
    parser.add_argument("--report", default="GNU_Burst_Scan_Report.csv")
    parser.add_argument("--sample-rate", type=float, default=100000.0)
    parser.add_argument("--symbol-rate", type=float, default=9600.0)
    parser.add_argument("--deviation", type=float, default=2500.0)
    parser.add_argument("--scan-min", type=int, default=-20000)
    parser.add_argument("--scan-max", type=int, default=20000)
    parser.add_argument("--scan-step", type=int, default=1)
    parser.add_argument("--phase-count", type=int, default=32)
    parser.add_argument("--detection-margin-db", type=float, default=4.0)
    parser.add_argument("--disable-boundary-recovery", action="store_true")
    parser.add_argument("--disable-flags-only-diagnostic", action="store_true")
    parser.add_argument("--normal-polarity-only", action="store_true")
    parser.add_argument("--disable-symbol-rate-search", action="store_true")
    parser.add_argument("--symbol-rate-min", type=float, default=9550.0)
    parser.add_argument("--symbol-rate-max", type=float, default=9650.0)
    parser.add_argument("--symbol-rate-step", type=float, default=1.0)
    parser.add_argument("--baud-search-phase-count", type=int, default=4)
    parser.add_argument("--disable-clock-recovery", action="store_true")
    parser.add_argument("--disable-fir-path", action="store_true")
    parser.add_argument("--no-preamble-full-search-top-n", type=int, default=3)
    parser.add_argument("--disable-soft-repair", action="store_true")
    parser.add_argument(
        "--soft-repair-library",
        default=str(_default_soft_repair_library()),
        help="Root containing the repository-local scionx advanced decoder package",
    )
    parser.add_argument("--soft-repair-max-header-bit-errors", type=int, default=24)
    parser.add_argument("--soft-repair-max-observations", type=int, default=32)
    parser.add_argument(
        "--soft-repair-effort", choices=("STANDARD", "MAXIMUM"), default="MAXIMUM"
    )
    parser.add_argument("--extended-search-min-power-excess-db", type=float, default=2.0)
    parser.add_argument("--recovery-max-header-bit-errors", type=int, default=16)
    parser.add_argument("--recovery-max-closing-flag-bit-errors", type=int, default=2)
    parser.add_argument("--recovery-max-destuff-violations", type=int, default=8)
    return parser


def main(argv: Iterable[str] | None = None) -> int:
    args = _argument_parser().parse_args(argv)
    decode_cs16_file(
        args.input_file,
        args.output,
        args.report,
        sample_rate=args.sample_rate,
        symbol_rate=args.symbol_rate,
        deviation_hz=args.deviation,
        scan_min_hz=args.scan_min,
        scan_max_hz=args.scan_max,
        scan_step_hz=args.scan_step,
        phase_count=args.phase_count,
        detection_margin_db=args.detection_margin_db,
        enable_boundary_recovery=not args.disable_boundary_recovery,
        recovery_max_header_bit_errors=args.recovery_max_header_bit_errors,
        recovery_max_closing_flag_bit_errors=(
            args.recovery_max_closing_flag_bit_errors
        ),
        recovery_max_destuff_violations=args.recovery_max_destuff_violations,
        enable_flags_only_diagnostic=not args.disable_flags_only_diagnostic,
        test_inverted_polarity=not args.normal_polarity_only,
        enable_symbol_rate_search=not args.disable_symbol_rate_search,
        symbol_rate_min=args.symbol_rate_min,
        symbol_rate_max=args.symbol_rate_max,
        symbol_rate_step=args.symbol_rate_step,
        baud_search_phase_count=args.baud_search_phase_count,
        enable_clock_recovery=not args.disable_clock_recovery,
        enable_fir_path=not args.disable_fir_path,
        no_preamble_full_search_top_n=args.no_preamble_full_search_top_n,
        enable_soft_repair=not args.disable_soft_repair,
        soft_repair_library=args.soft_repair_library,
        soft_repair_max_header_bit_errors=(
            args.soft_repair_max_header_bit_errors
        ),
        soft_repair_max_observations=args.soft_repair_max_observations,
        soft_repair_effort=args.soft_repair_effort,
        extended_search_min_power_excess_db=(
            args.extended_search_min_power_excess_db
        ),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
