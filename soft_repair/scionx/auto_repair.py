"""Conservative HDLC repair for the SCION-X fixed-length beacon frame.

The repair search is deliberately bounded.  It first finds an exact run of
closing HDLC flags, then applies strict bit de-stuffing to the raw decisions.
Only the least-confident raw decisions are eligible for at most two flips.
A repair is accepted only when exactly one distinct 274-byte frame satisfies
all of these constraints:

* 2192 de-stuffed bits (272-byte AX.25 payload plus 2-byte FCS),
* the fixed 16-byte SCION-X address/control/PID anchor,
* an exact run of at least three closing 0x7E flags, and
* CRC-16/X.25.

CRC is used as a validator, not as permission for an unbounded brute-force
search.  When no unique solution exists, the caller retains the original
decisions for manual review.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from typing import Iterable, Sequence

import numpy as np

from .hdlc import crc16_x25


FLAG_BITS = (0, 1, 1, 1, 1, 1, 1, 0)
EXPECTED_FRAME_BITS = 274 * 8
PAYLOAD_BYTES = 272
HARD_ANCHOR = bytes.fromhex("849c6086aa4060849c60a686b0e103f0")


@dataclass(frozen=True)
class DestuffResult:
    valid: bool
    bits: tuple[int, ...]
    source_indices: tuple[int, ...]
    excluded_indices: tuple[int, ...]
    violation_index: int | None = None


@dataclass(frozen=True)
class RepairCandidate:
    frame: bytes
    flag_start: int
    flip_indices: tuple[int, ...]
    excluded_indices: tuple[int, ...]
    source_indices: tuple[int, ...]
    cost: float


@dataclass(frozen=True)
class AutoRepairResult:
    status: str
    reason: str
    original_stuffing_aligned: bool
    original_crc_pass: bool
    hard_anchor_match: bool
    closing_flag_found: bool
    selected_flag_start: int | None
    destuffed_bit_count: int
    excluded_indices: tuple[int, ...]
    flip_indices: tuple[int, ...]
    source_indices: tuple[int, ...]
    candidate_pool: tuple[int, ...]
    candidate_evaluations: int
    crc_solution_count: int
    repair_cost: float | None
    frame: bytes | None

    @property
    def auto_repaired(self) -> bool:
        return self.status == "AUTO_REPAIRED_UNIQUE"

    def to_summary(self) -> dict:
        return {
            "status": self.status,
            "reason": self.reason,
            "original_stuffing_aligned": self.original_stuffing_aligned,
            "original_crc_pass": self.original_crc_pass,
            "hard_anchor_match": self.hard_anchor_match,
            "closing_flag_found": self.closing_flag_found,
            "selected_flag_start": self.selected_flag_start,
            "destuffed_bit_count": self.destuffed_bit_count,
            "excluded_indices": list(self.excluded_indices),
            "flip_indices": list(self.flip_indices),
            "candidate_pool": list(self.candidate_pool),
            "candidate_evaluations": self.candidate_evaluations,
            "crc_solution_count": self.crc_solution_count,
            "repair_cost": self.repair_cost,
            "frame_hex": self.frame.hex() if self.frame is not None else "",
        }


def bits_lsb_first(data: bytes) -> tuple[int, ...]:
    return tuple((byte >> bit) & 1 for byte in data for bit in range(8))


def bits_to_bytes(bits: Sequence[int]) -> bytes:
    if len(bits) % 8:
        raise ValueError("The bit count must be divisible by 8.")
    output = bytearray(len(bits) // 8)
    for byte_index in range(len(output)):
        value = 0
        for bit_index in range(8):
            value |= int(bits[byte_index * 8 + bit_index]) << bit_index
        output[byte_index] = value
    return bytes(output)


def stuff_bits(bits: Sequence[int]) -> tuple[int, ...]:
    output: list[int] = []
    run = 0
    for value in bits:
        bit = int(value)
        output.append(bit)
        if bit:
            run += 1
            if run == 5:
                output.append(0)
                run = 0
        else:
            run = 0
    return tuple(output)


def strict_destuff(raw_bits: Sequence[int]) -> DestuffResult:
    output: list[int] = []
    sources: list[int] = []
    excluded: list[int] = []
    run = 0
    for raw_index, value in enumerate(raw_bits):
        bit = int(value)
        if run == 5:
            if bit != 0:
                return DestuffResult(
                    False, tuple(output), tuple(sources), tuple(excluded), raw_index
                )
            excluded.append(raw_index)
            run = 0
            continue
        output.append(bit)
        sources.append(raw_index)
        run = run + 1 if bit else 0
    return DestuffResult(True, tuple(output), tuple(sources), tuple(excluded))


def find_closing_flag_starts(
    raw_bits: Sequence[int],
    expected_destuffed_bits: int = EXPECTED_FRAME_BITS,
    minimum_flags: int = 3,
    search_extra_bits: int = 128,
) -> tuple[int, ...]:
    bits = tuple(int(value) for value in raw_bits)
    pattern = FLAG_BITS * minimum_flags
    low = expected_destuffed_bits
    high = min(len(bits) - len(pattern), expected_destuffed_bits + search_extra_bits)
    if high < low:
        return tuple()

    starts: list[int] = []
    for index in range(low, high + 1):
        if bits[index:index + len(pattern)] != pattern:
            continue
        # Record only the first flag in a contiguous run, rather than the
        # second/third flag that also happens to begin three valid flags.
        if index >= 8 and bits[index - 8:index] == FLAG_BITS:
            continue
        starts.append(index)
    return tuple(starts)


def _evaluate(
    corrected_raw: np.ndarray,
    flag_start: int,
    flip_indices: tuple[int, ...],
    soft_distance: np.ndarray,
) -> RepairCandidate | None:
    destuffed = strict_destuff(corrected_raw[:flag_start])
    if not destuffed.valid or len(destuffed.bits) != EXPECTED_FRAME_BITS:
        return None
    frame = bits_to_bytes(destuffed.bits)
    if frame[:len(HARD_ANCHOR)] != HARD_ANCHOR:
        return None
    received_fcs = int.from_bytes(frame[PAYLOAD_BYTES:PAYLOAD_BYTES + 2], "little")
    if crc16_x25(frame[:PAYLOAD_BYTES]) != received_fcs:
        return None
    cost = float(sum(float(soft_distance[index]) for index in flip_indices))
    return RepairCandidate(
        frame=frame,
        flag_start=flag_start,
        flip_indices=flip_indices,
        excluded_indices=destuffed.excluded_indices,
        source_indices=destuffed.source_indices,
        cost=cost,
    )


def _header_raw_mismatches(raw_bits: np.ndarray) -> tuple[int, ...]:
    expected = stuff_bits(bits_lsb_first(HARD_ANCHOR))
    count = min(len(expected), len(raw_bits))
    return tuple(
        index for index in range(count) if int(raw_bits[index]) != expected[index]
    )


def _protocol_priority_indices(raw_bits: np.ndarray, limit: int) -> set[int]:
    """Return indices implicated in one-bit-near HDLC stuffing patterns."""
    priority: set[int] = set()
    upper = min(limit, len(raw_bits))
    for index in range(5, upper):
        previous = raw_bits[index - 5:index]
        ones = int(np.count_nonzero(previous))
        current = int(raw_bits[index])
        if ones == 5 and current == 1:
            priority.update(range(index - 5, index + 1))
        elif ones == 4 and current == 0:
            priority.add(index)
            priority.update(
                position
                for position in range(index - 5, index)
                if int(raw_bits[position]) == 0
            )
    return priority


def _build_candidate_pool(
    raw_bits: np.ndarray,
    soft_distance: np.ndarray,
    flag_starts: Sequence[int],
    maximum_pool_size: int,
) -> tuple[int, ...]:
    limit = max(flag_starts) if flag_starts else len(raw_bits)
    limit = min(limit, len(raw_bits))
    mandatory = set(_header_raw_mismatches(raw_bits[:limit]))
    protocol = _protocol_priority_indices(raw_bits, limit)

    ranked = sorted(
        range(limit),
        key=lambda index: (
            0 if index in mandatory else 1,
            float(soft_distance[index]),
            0 if index in protocol else 1,
            index,
        ),
    )
    selected = ranked[:maximum_pool_size]
    return tuple(sorted(selected))


def _original_alignment(
    raw_bits: np.ndarray, flag_starts: Sequence[int]
) -> tuple[DestuffResult | None, int | None, bool, bool]:
    best: tuple[tuple[int, int], DestuffResult, int, bool] | None = None
    for flag_start in flag_starts:
        destuffed = strict_destuff(raw_bits[:flag_start])
        aligned = destuffed.valid and len(destuffed.bits) == EXPECTED_FRAME_BITS
        hard_anchor = False
        if aligned:
            frame = bits_to_bytes(destuffed.bits)
            hard_anchor = frame[:len(HARD_ANCHOR)] == HARD_ANCHOR
        score = (0 if aligned else 1, abs(flag_start - (EXPECTED_FRAME_BITS + 8)))
        if best is None or score < best[0]:
            best = (score, destuffed, flag_start, hard_anchor)
    if best is None:
        return None, None, False, False
    return best[1], best[2], best[1].valid and len(best[1].bits) == EXPECTED_FRAME_BITS, best[3]


def auto_repair_frame(
    raw_bits: Sequence[int],
    soft_distance: Sequence[float],
    *,
    maximum_flips: int = 2,
    maximum_pool_size: int = 32,
    maximum_candidate_evaluations: int = 2000,
) -> AutoRepairResult:
    """Attempt a bounded, unique CRC-valid repair of one raw HDLC frame."""
    bits = np.asarray(raw_bits, dtype=np.uint8).reshape(-1)
    distances = np.asarray(soft_distance, dtype=np.float64).reshape(-1)
    if bits.size != distances.size:
        raise ValueError("raw_bits and soft_distance must have equal lengths.")
    if maximum_flips < 0 or maximum_flips > 2:
        raise ValueError("maximum_flips must be between 0 and 2.")
    if maximum_pool_size <= 0 or maximum_candidate_evaluations <= 0:
        raise ValueError("Search limits must be positive.")
    distances = np.nan_to_num(distances, nan=np.inf, posinf=np.inf, neginf=np.inf)

    flag_starts = find_closing_flag_starts(bits)
    original, original_flag, original_aligned, original_anchor = _original_alignment(
        bits, flag_starts
    )
    original_crc = False
    original_frame: bytes | None = None
    if original_aligned and original is not None:
        original_frame = bits_to_bytes(original.bits)
        received_fcs = int.from_bytes(
            original_frame[PAYLOAD_BYTES:PAYLOAD_BYTES + 2], "little"
        )
        original_crc = (
            crc16_x25(original_frame[:PAYLOAD_BYTES]) == received_fcs
            and original_anchor
        )

    if not flag_starts:
        return AutoRepairResult(
            status="UNRECOVERABLE",
            reason="No exact run of three closing 0x7E flags was found near the frame end.",
            original_stuffing_aligned=False,
            original_crc_pass=False,
            hard_anchor_match=False,
            closing_flag_found=False,
            selected_flag_start=None,
            destuffed_bit_count=len(original.bits) if original is not None else 0,
            excluded_indices=tuple(),
            flip_indices=tuple(),
            source_indices=tuple(),
            candidate_pool=tuple(),
            candidate_evaluations=0,
            crc_solution_count=0,
            repair_cost=None,
            frame=None,
        )

    pool = _build_candidate_pool(bits, distances, flag_starts, maximum_pool_size)
    solutions: dict[bytes, RepairCandidate] = {}
    evaluations = 0

    for flip_count in range(maximum_flips + 1):
        for flip_indices in combinations(pool, flip_count):
            if evaluations >= maximum_candidate_evaluations:
                break
            evaluations += 1
            corrected = bits.copy()
            if flip_indices:
                corrected[list(flip_indices)] ^= 1
            for flag_start in flag_starts:
                candidate = _evaluate(
                    corrected, flag_start, tuple(flip_indices), distances
                )
                if candidate is None:
                    continue
                previous = solutions.get(candidate.frame)
                if previous is None or (len(candidate.flip_indices), candidate.cost) < (
                    len(previous.flip_indices), previous.cost
                ):
                    solutions[candidate.frame] = candidate
        if evaluations >= maximum_candidate_evaluations:
            break

    if len(solutions) == 1:
        selected = next(iter(solutions.values()))
        status = "ORIGINAL_CRC_PASS" if not selected.flip_indices else "AUTO_REPAIRED_UNIQUE"
        reason = (
            "The original decisions satisfy HDLC alignment, the hard anchor, and CRC."
            if status == "ORIGINAL_CRC_PASS"
            else "Exactly one bounded soft-bit candidate satisfies HDLC alignment, the hard anchor, and CRC."
        )
        return AutoRepairResult(
            status=status,
            reason=reason,
            original_stuffing_aligned=original_aligned,
            original_crc_pass=original_crc,
            hard_anchor_match=True,
            closing_flag_found=True,
            selected_flag_start=selected.flag_start,
            destuffed_bit_count=EXPECTED_FRAME_BITS,
            excluded_indices=selected.excluded_indices,
            flip_indices=selected.flip_indices,
            source_indices=selected.source_indices,
            candidate_pool=pool,
            candidate_evaluations=evaluations,
            crc_solution_count=1,
            repair_cost=selected.cost,
            frame=selected.frame,
        )

    if len(solutions) > 1:
        status = "AMBIGUOUS_CRC_SOLUTIONS"
        reason = (
            f"{len(solutions)} distinct bounded candidates passed CRC; no flips were auto-applied."
        )
    elif original_aligned:
        status = "STUFFING_ALIGNED_CRC_FAIL"
        reason = (
            "Strict HDLC de-stuffing reaches the closing flags at 2192 bits, but no unique <=2-flip CRC solution exists."
        )
    else:
        status = "UNRECOVERABLE"
        reason = (
            "No bounded <=2-flip candidate satisfied frame length, hard anchor, closing flags, and CRC."
        )

    # Safe fallback: when strict no-flip de-stuffing is structurally aligned,
    # pre-fill only those protocol-mandated exclusions.  Never auto-apply an
    # unverified flip in an ambiguous/no-solution result.
    fallback_excluded = original.excluded_indices if original_aligned and original else tuple()
    fallback_sources = original.source_indices if original_aligned and original else tuple()
    return AutoRepairResult(
        status=status,
        reason=reason,
        original_stuffing_aligned=original_aligned,
        original_crc_pass=original_crc,
        hard_anchor_match=original_anchor,
        closing_flag_found=True,
        selected_flag_start=original_flag,
        destuffed_bit_count=len(original.bits) if original is not None else 0,
        excluded_indices=fallback_excluded,
        flip_indices=tuple(),
        source_indices=fallback_sources,
        candidate_pool=pool,
        candidate_evaluations=evaluations,
        crc_solution_count=len(solutions),
        repair_cost=None,
        frame=original_frame if original_aligned else None,
    )
