"""Hybrid burst detection for SCION-X CS16 complex-IQ recordings.

The front end combines the original adaptive IQ-power detector with a
carrier-offset-tolerant matched correlation against the raw HDLC sequence
``4 x 0x7E + BN0CU -> BN0SCX + Control/PID``.  It only selects candidate
windows; :class:`scionx_burst_scan.StreamingBurstScanner` remains responsible
for demodulation, route arbitration, raw-HDLC output and optional soft repair.

Large CS16 files are read through ``numpy.memmap`` in bounded chunks.  The GNU
Radio stream is therefore not duplicated in memory when ``input_file`` is
available.
"""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass
from pathlib import Path
from typing import Protocol

import numpy as np
from scipy.signal import fftconvolve, find_peaks, firwin

from scionx_burst_scan import (
    FLAG_BITS,
    SCIONX_HEADER_BITS,
    StreamingBurstScanner,
)


class _SampleSource(Protocol):
    sample_count: int

    def read(self, start: int, stop: int) -> np.ndarray:
        """Return normalized complex64 IQ samples in ``[start, stop)``."""


class _Cs16FileSource:
    def __init__(self, filename: str) -> None:
        self.filename = Path(filename)
        self._raw = np.memmap(self.filename, dtype="<i2", mode="r")
        self.sample_count = int(self._raw.size // 2)

    def read(self, start: int, stop: int) -> np.ndarray:
        first = max(0, min(self.sample_count, int(start)))
        last = max(first, min(self.sample_count, int(stop)))
        raw = np.asarray(self._raw[2 * first:2 * last], dtype=np.float32)
        return ((raw[0::2] + 1j * raw[1::2]) / 32767.0).astype(np.complex64)


class _ArraySource:
    def __init__(self, samples: np.ndarray) -> None:
        self._samples = np.asarray(samples, dtype=np.complex64).ravel()
        self.sample_count = int(self._samples.size)

    def read(self, start: int, stop: int) -> np.ndarray:
        first = max(0, min(self.sample_count, int(start)))
        last = max(first, min(self.sample_count, int(stop)))
        return self._samples[first:last].copy()


@dataclass
class Cs16BurstDetection:
    start: int
    stop: int
    origin: int
    method: str
    score: float
    noise_floor_db: float
    power_excess_db: float = float("nan")
    filter_method: str = ""
    polarity: str = ""


def _stuff_payload_bits(bits: np.ndarray) -> np.ndarray:
    """Apply HDLC zero stuffing to frame bits (flags are handled separately)."""
    output: list[int] = []
    ones = 0
    for value in np.asarray(bits, dtype=np.uint8).ravel():
        bit = int(value)
        output.append(bit)
        if bit:
            ones += 1
            if ones == 5:
                output.append(0)
                ones = 0
        else:
            ones = 0
    return np.asarray(output, dtype=np.uint8)


def _header_template(sample_rate: float, symbol_rate: float) -> np.ndarray:
    """Return a Gaussian-smoothed, non-integer-SPS SCION-X header template."""
    bits = np.concatenate((np.tile(FLAG_BITS, 4), _stuff_payload_bits(SCIONX_HEADER_BITS)))
    samples_per_symbol = float(sample_rate) / float(symbol_rate)
    sample_count = max(1, int(np.ceil(bits.size * samples_per_symbol)))
    bit_index = np.floor(
        np.arange(sample_count, dtype=np.float64) / samples_per_symbol
    ).astype(np.int64)
    bit_index = np.minimum(bit_index, bits.size - 1)
    levels = 2.0 * bits[bit_index].astype(np.float64) - 1.0

    # The known modulation uses Gaussian BT=0.5.  A short smoothing kernel is
    # intentionally used: acquisition must remain tolerant of receiver filter
    # differences while still matching the transmitted transition shape.
    bt = 0.5
    sigma = max(0.5, np.sqrt(np.log(2.0)) * samples_per_symbol / (2.0 * np.pi * bt))
    radius = max(2, int(np.ceil(4.0 * sigma)))
    positions = np.arange(-radius, radius + 1, dtype=np.float64)
    taps = np.exp(-0.5 * (positions / sigma) ** 2)
    taps /= np.sum(taps)
    shaped = np.convolve(levels, taps, mode="same")
    shaped -= float(np.mean(shaped))
    return shaped


def _normalised_match(signal: np.ndarray, template: np.ndarray) -> np.ndarray:
    """Zero-mean normalized cross-correlation using FFT convolution."""
    values = np.asarray(signal, dtype=np.float64).ravel()
    pattern = np.asarray(template, dtype=np.float64).ravel()
    if values.size < pattern.size or pattern.size == 0:
        return np.empty(0, dtype=np.float64)
    pattern = pattern - float(np.mean(pattern))
    numerator = fftconvolve(values, pattern[::-1], mode="valid")
    cumulative = np.concatenate(([0.0], np.cumsum(values)))
    cumulative_sq = np.concatenate(([0.0], np.cumsum(values * values)))
    length = pattern.size
    window_sum = cumulative[length:] - cumulative[:-length]
    window_sq = cumulative_sq[length:] - cumulative_sq[:-length]
    variance_energy = np.maximum(window_sq - window_sum * window_sum / length, 0.0)
    denominator = np.sqrt(np.sum(pattern * pattern)) * np.sqrt(variance_energy)
    return numerator / (denominator + 1e-12)


def _robust_z(values: np.ndarray) -> np.ndarray:
    source = np.asarray(values, dtype=np.float64)
    if source.size == 0:
        return source.copy()
    centre = float(np.median(source))
    mad = float(np.median(np.abs(source - centre)))
    standard_deviation = float(np.std(source))
    scale = max(1.4826 * mad, 0.25 * standard_deviation, 1e-9)
    return (source - centre) / scale


def _power_regions(
    source: _SampleSource,
    sample_rate: float,
    *,
    detection_margin_db: float,
    block_size: int,
    pretrigger_ms: float,
    end_hold_ms: float,
    minimum_burst_ms: float,
    maximum_burst_ms: float,
    chunk_samples: int,
) -> tuple[list[Cs16BurstDetection], float]:
    """Reproduce the adaptive block-power detector without buffering the file."""
    history: deque[float] = deque(maxlen=1000)
    pretrigger_blocks = max(1, int(round(pretrigger_ms * 1e-3 * sample_rate / block_size)))
    end_hold_blocks = max(1, int(round(end_hold_ms * 1e-3 * sample_rate / block_size)))
    minimum_samples = int(round(minimum_burst_ms * 1e-3 * sample_rate))
    maximum_samples = int(round(maximum_burst_ms * 1e-3 * sample_rate))
    regions: list[Cs16BurstDetection] = []
    quiet_powers: list[float] = []

    active = False
    below_count = 0
    active_start = 0
    active_noise = float("nan")
    active_max_excess = float("-inf")
    block_index = 0
    pending = np.empty(0, dtype=np.complex64)
    read_chunk = max(block_size, int(chunk_samples // block_size) * block_size)

    def finish_region(stop: int) -> None:
        nonlocal active, below_count, active_max_excess
        stop = min(source.sample_count, max(active_start, int(stop)))
        if stop - active_start >= minimum_samples:
            regions.append(
                Cs16BurstDetection(
                    start=active_start,
                    stop=stop,
                    origin=active_start,
                    method="POWER_RMS",
                    score=float(active_max_excess),
                    noise_floor_db=float(active_noise),
                    power_excess_db=float(active_max_excess),
                )
            )
        active = False
        below_count = 0
        active_max_excess = float("-inf")

    for chunk_start in range(0, source.sample_count, read_chunk):
        values = source.read(chunk_start, min(source.sample_count, chunk_start + read_chunk))
        if pending.size:
            values = np.concatenate((pending, values))
        complete = values.size // block_size
        if complete:
            blocks = values[:complete * block_size].reshape(complete, block_size)
            power_db_values = 10.0 * np.log10(
                np.mean(np.abs(blocks) ** 2, axis=1).astype(np.float64) + 1e-20
            )
            for power_db in power_db_values:
                power_db = float(power_db)
                noise_floor = float(np.median(history)) if history else power_db
                excess = power_db - noise_floor
                is_high = len(history) >= 30 and excess > float(detection_margin_db)
                if not active:
                    if is_high:
                        active = True
                        below_count = 0
                        active_start = max(0, (block_index - pretrigger_blocks) * block_size)
                        active_noise = noise_floor
                        active_max_excess = excess
                    else:
                        history.append(power_db)
                        quiet_powers.append(power_db)
                else:
                    active_max_excess = max(active_max_excess, excess)
                    if power_db <= noise_floor + float(detection_margin_db):
                        below_count += 1
                    else:
                        below_count = 0
                    active_samples = (block_index + 1) * block_size - active_start
                    if below_count >= end_hold_blocks or active_samples >= maximum_samples:
                        finish_region((block_index + 1) * block_size)
                block_index += 1
        pending = values[complete * block_size:].copy()

    if pending.size:
        padded = np.pad(pending, (0, block_size - pending.size))
        power_db = 10.0 * np.log10(float(np.mean(np.abs(padded) ** 2)) + 1e-20)
        noise_floor = float(np.median(history)) if history else power_db
        excess = power_db - noise_floor
        if active:
            active_max_excess = max(active_max_excess, excess)
        block_index += 1
    if active:
        finish_region(source.sample_count)

    global_noise_floor = float(np.median(quiet_powers)) if quiet_powers else float("nan")
    return regions, global_noise_floor


def _header_detections(
    source: _SampleSource,
    sample_rate: float,
    symbol_rate: float,
    *,
    z_threshold: float,
    lowpass_hz: float,
    enable_fir_path: bool,
    pretrigger_ms: float,
    posttrigger_ms: float,
    chunk_seconds: float,
    noise_floor_db: float,
) -> list[Cs16BurstDetection]:
    """Search demodulated IQ for a known raw-HDLC header in bounded chunks."""
    template = _header_template(sample_rate, symbol_rate)
    channel_taps = firwin(81, lowpass_hz, fs=sample_rate, window=("kaiser", 6.76))
    before = int(round(pretrigger_ms * 1e-3 * sample_rate))
    after = int(round(posttrigger_ms * 1e-3 * sample_rate))
    minimum_gap = max(1, int(round(0.35 * sample_rate)))
    core_size = max(int(round(chunk_seconds * sample_rate)), minimum_gap)
    padding = template.size + channel_taps.size + 8
    raw_hits: list[Cs16BurstDetection] = []

    for core_start in range(0, source.sample_count, core_size):
        core_stop = min(source.sample_count, core_start + core_size)
        read_start = max(0, core_start - padding)
        read_stop = min(source.sample_count, core_stop + padding)
        iq = source.read(read_start, read_stop)
        if iq.size < template.size + 2:
            continue
        discriminator_hz = (
            np.angle(iq[1:] * np.conj(iq[:-1])).astype(np.float64)
            * sample_rate
            / (2.0 * np.pi)
        )
        raw_filtered = fftconvolve(discriminator_hz, channel_taps, mode="same")
        routes: list[tuple[str, np.ndarray]] = [("RAW", raw_filtered)]
        if enable_fir_path:
            routes.append((
                "FIR_025_050_025",
                np.convolve(raw_filtered, [0.25, 0.50, 0.25], mode="same"),
            ))

        for filter_method, values in routes:
            correlation = _normalised_match(values, template)
            z = _robust_z(correlation)
            for polarity, polarity_z in (("NORMAL", z), ("INVERTED", -z)):
                peaks, properties = find_peaks(
                    polarity_z,
                    height=float(z_threshold),
                    distance=minimum_gap,
                )
                heights = properties.get("peak_heights", np.empty(0))
                for peak, score in zip(peaks, heights):
                    origin = read_start + int(peak)
                    if origin < core_start or origin >= core_stop:
                        continue
                    raw_hits.append(
                        Cs16BurstDetection(
                            start=max(0, origin - before),
                            stop=min(source.sample_count, origin + after),
                            origin=origin,
                            method=f"HEADER_{filter_method}_{polarity}",
                            score=float(score),
                            noise_floor_db=noise_floor_db,
                            filter_method=filter_method,
                            polarity=polarity,
                        )
                    )

    # The same physical header can be detected by both filter routes and nearby
    # correlation lobes.  Keep the strongest observation inside 350 ms.
    selected: list[Cs16BurstDetection] = []
    for item in sorted(raw_hits, key=lambda candidate: candidate.score, reverse=True):
        if any(abs(item.origin - previous.origin) <= minimum_gap for previous in selected):
            continue
        selected.append(item)
    selected.sort(key=lambda item: item.start)
    return selected


def detect_cs16_bursts(
    source: _SampleSource,
    sample_rate: float,
    symbol_rate: float,
    *,
    detection_margin_db: float = 4.0,
    header_z_threshold: float = 5.0,
    enable_header_detector: bool = True,
    enable_fir_path: bool = True,
    lowpass_hz: float = 6500.0,
    block_size: int = 1000,
    pretrigger_ms: float = 100.0,
    end_hold_ms: float = 350.0,
    minimum_burst_ms: float = 80.0,
    maximum_burst_ms: float = 5000.0,
    detector_posttrigger_ms: float = 500.0,
    header_chunk_seconds: float = 5.0,
) -> list[Cs16BurstDetection]:
    """Return merged power and known-header candidate windows."""
    power, global_noise_floor = _power_regions(
        source,
        sample_rate,
        detection_margin_db=detection_margin_db,
        block_size=block_size,
        pretrigger_ms=pretrigger_ms,
        end_hold_ms=end_hold_ms,
        minimum_burst_ms=minimum_burst_ms,
        maximum_burst_ms=maximum_burst_ms,
        chunk_samples=max(block_size, int(round(header_chunk_seconds * sample_rate))),
    )
    headers = (
        _header_detections(
            source,
            sample_rate,
            symbol_rate,
            z_threshold=header_z_threshold,
            lowpass_hz=lowpass_hz,
            enable_fir_path=enable_fir_path,
            pretrigger_ms=pretrigger_ms,
            posttrigger_ms=detector_posttrigger_ms,
            chunk_seconds=header_chunk_seconds,
            noise_floor_db=global_noise_floor,
        )
        if enable_header_detector
        else []
    )

    detections = list(headers)
    for energy in power:
        overlapping = [
            item for item in detections
            if item.start < energy.stop and energy.start < item.stop
        ]
        if overlapping:
            best = max(overlapping, key=lambda item: item.score)
            best.start = min(best.start, energy.start)
            best.stop = max(best.stop, energy.stop)
            best.method += "+POWER_RMS"
            best.power_excess_db = energy.power_excess_db
            if np.isfinite(energy.noise_floor_db):
                best.noise_floor_db = energy.noise_floor_db
            continue
        detections.append(energy)

    detections.sort(key=lambda item: item.start)
    return detections


class Cs16HybridBurstScanner:
    """Use hybrid CS16 acquisition, then reuse the shared demodulation core."""

    def __init__(
        self,
        *args: object,
        input_file: str = "",
        enable_header_detector: bool = True,
        header_z_threshold: float = 5.0,
        header_chunk_seconds: float = 5.0,
        detector_pretrigger_ms: float = 100.0,
        detector_posttrigger_ms: float = 500.0,
        **kwargs: object,
    ) -> None:
        kwargs["input_mode"] = "IQ"
        kwargs["emit_decoded_candidate_files"] = False
        self.decoder = StreamingBurstScanner(*args, **kwargs)
        self.input_file = str(input_file)
        self.enable_header_detector = bool(enable_header_detector)
        self.header_z_threshold = float(header_z_threshold)
        self.header_chunk_seconds = float(header_chunk_seconds)
        self.detector_pretrigger_ms = float(detector_pretrigger_ms)
        self.detector_posttrigger_ms = float(detector_posttrigger_ms)
        self._file_source_available = Path(self.input_file).is_file()
        self._chunks: list[np.ndarray] = []
        self._finished = False
        self._received_samples = 0

        if self.header_z_threshold <= 0 or self.header_chunk_seconds <= 0:
            raise ValueError("CS16 header Z threshold and chunk duration must be positive")

    def feed(self, samples: np.ndarray) -> None:
        if self._finished:
            raise RuntimeError("Cannot feed CS16 samples after finish()")
        values = np.asarray(samples, dtype=np.complex64).ravel()
        self._received_samples += int(values.size)
        # With a valid source path, finish() reopens the original CS16 through
        # memmap.  Buffering is retained only for unit tests/custom flowgraphs.
        if values.size and not self._file_source_available:
            self._chunks.append(values.copy())

    def finish(self) -> None:
        if self._finished:
            return
        self._finished = True
        if self._file_source_available:
            source: _SampleSource = _Cs16FileSource(self.input_file)
        elif self._chunks:
            source = _ArraySource(np.concatenate(self._chunks))
        else:
            return
        self._chunks.clear()

        detections = detect_cs16_bursts(
            source,
            self.decoder.sample_rate,
            self.decoder.symbol_rate,
            detection_margin_db=self.decoder.detection_margin_db,
            header_z_threshold=self.header_z_threshold,
            enable_header_detector=self.enable_header_detector,
            enable_fir_path=self.decoder.enable_fir_path,
            block_size=self.decoder.block_size,
            pretrigger_ms=self.detector_pretrigger_ms,
            end_hold_ms=(
                self.decoder.end_hold_blocks
                * self.decoder.block_size
                * 1000.0
                / self.decoder.sample_rate
            ),
            minimum_burst_ms=self.decoder.minimum_samples * 1000.0 / self.decoder.sample_rate,
            maximum_burst_ms=self.decoder.maximum_samples * 1000.0 / self.decoder.sample_rate,
            detector_posttrigger_ms=self.detector_posttrigger_ms,
            header_chunk_seconds=self.header_chunk_seconds,
        )
        print(
            "SCION-X CS16 hybrid detector: {} candidate(s) [{}]".format(
                len(detections),
                ", ".join(
                    "{:.3f}s {} score={:.2f}".format(
                        item.start / self.decoder.sample_rate,
                        item.method,
                        item.score,
                    )
                    for item in detections
                ),
            )
        )
        if not detections:
            self.decoder._start_outputs()
            return

        for item in detections:
            segment = source.read(item.start, item.stop)
            header_detected = item.method.startswith("HEADER_")
            self.decoder.process_explicit_burst(
                segment,
                item.start,
                noise_floor_db=item.noise_floor_db,
                force_extended_search=header_detected,
                force_full_search=header_detected,
                detector_method=item.method,
                detector_score=item.score,
            )
        self.decoder.finish()
