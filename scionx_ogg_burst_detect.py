"""SCION-X OGG front-end with matched-header and spectral-energy detection.

The detector is adapted from ``decode_scionx_ogg.m`` and
``run_scionx_destuff_ver2.py``.  It deliberately stops at burst selection:
the selected audio is sent to :class:`StreamingBurstScanner`, whose formal
``GNU_GFSK_Raw_HEX.txt`` output retains HDLC flags, stuffed zeroes and FCS for
``decode_gnuradio_hex.m``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from scionx_burst_scan import StreamingBurstScanner


FLAGS4 = bytes([0x7E] * 4)
HEADER_ADDR = bytes.fromhex("849c6086aa4060849c60a686b0e1")


def _bits_lsb_first(data: bytes) -> np.ndarray:
    values = np.frombuffer(data, dtype=np.uint8)
    return ((values[:, None] >> np.arange(8, dtype=np.uint8)) & 1).reshape(-1)


def _moving_mean(values: np.ndarray, window_size: int) -> np.ndarray:
    """Centred shrinking-window mean, matching MATLAB ``movmean`` edges."""
    source = np.asarray(values, dtype=np.float64).ravel()
    if source.size == 0:
        return source.copy()
    before = window_size // 2
    after = window_size - before - 1
    cumulative = np.concatenate(([0.0], np.cumsum(source)))
    index = np.arange(source.size)
    low = np.maximum(0, index - before)
    high = np.minimum(source.size, index + after + 1)
    return (cumulative[high] - cumulative[low]) / (high - low)


def restore_baseline(
    samples: np.ndarray,
    iterations: int = 7,
    half_window: int = 1000,
) -> np.ndarray:
    """Decision-directed baseline restoration used by the legacy OGG path."""
    values = np.asarray(samples, dtype=np.float64).ravel()
    if values.size == 0:
        return values.copy()
    amplitude = float(np.median(np.abs(values)))
    baseline = np.zeros_like(values)
    for _ in range(max(0, int(iterations))):
        compensated = values - baseline
        decisions = np.sign(compensated)
        decisions[decisions == 0] = 1
        residual = values - amplitude * decisions
        baseline = _moving_mean(residual, 2 * max(1, int(half_window)) + 1)
    return values - baseline


def _header_template(sample_rate: float, symbol_rate: float) -> np.ndarray:
    samples_per_symbol = sample_rate / symbol_rate
    rounded = int(round(samples_per_symbol))
    if rounded < 1 or abs(samples_per_symbol - rounded) > 1e-6:
        raise ValueError(
            "OGG header detector requires an integer samples/symbol ratio; "
            f"received {samples_per_symbol:.9g}"
        )
    bits = np.concatenate((_bits_lsb_first(FLAGS4), _bits_lsb_first(HEADER_ADDR)))
    return np.repeat(2.0 * bits.astype(np.float64) - 1.0, rounded)


def _normalised_match(signal: np.ndarray, template: np.ndarray) -> np.ndarray:
    if signal.size < template.size:
        return np.empty(0, dtype=np.float64)
    values = np.asarray(signal, dtype=np.float64)
    pattern = np.asarray(template, dtype=np.float64)
    numerator = np.correlate(values, pattern, mode="valid")
    cumulative = np.concatenate(([0.0], np.cumsum(values * values)))
    energy = cumulative[pattern.size:] - cumulative[:-pattern.size]
    denominator = np.sqrt(np.sum(pattern * pattern)) * np.sqrt(np.maximum(energy, 0.0))
    return numerator / (denominator + 1e-12)


def _energy_regions(
    samples: np.ndarray,
    sample_rate: float,
    threshold: float,
) -> list[tuple[int, int, float]]:
    """Return guarded 300 Hz..12 kHz short-time spectral-energy regions."""
    window = max(128, int(round(0.020 * sample_rate)))
    hop = max(32, int(round(0.010 * sample_rate)))
    if samples.size < window:
        return []
    positions = np.arange(0, samples.size - window + 1, hop, dtype=np.int64)
    frequency = np.fft.rfftfreq(window, d=1.0 / sample_rate)
    band = (frequency >= 300.0) & (frequency <= min(12000.0, sample_rate * 0.45))
    taper = np.hanning(window)
    offsets = np.arange(window, dtype=np.int64)
    power = np.empty(positions.size, dtype=np.float64)
    for first in range(0, positions.size, 256):
        batch_positions = positions[first:first + 256]
        frames = samples[batch_positions[:, None] + offsets[None, :]].copy()
        frames -= np.mean(frames, axis=1, keepdims=True)
        spectrum = np.fft.rfft(frames * taper[None, :], axis=1)
        power[first:first + batch_positions.size] = np.sum(
            np.abs(spectrum[:, band]) ** 2,
            axis=1,
        )
    median = float(np.median(power))
    mad = float(np.median(np.abs(power - median)))
    scale = max(1.4826 * mad, median * 1e-3, 1e-15)
    robust_z = (power - median) / scale
    active = np.flatnonzero(robust_z >= float(threshold))
    if active.size == 0:
        return []

    maximum_gap = max(1, int(round(0.030 * sample_rate / hop)))
    groups: list[list[int]] = [[int(active[0])]]
    for value in active[1:]:
        index = int(value)
        if index - groups[-1][-1] <= maximum_gap:
            groups[-1].append(index)
        else:
            groups.append([index])

    minimum_duration = int(round(0.035 * sample_rate))
    guard = int(round(0.100 * sample_rate))
    regions: list[tuple[int, int, float]] = []
    for group in groups:
        start = int(positions[group[0]])
        stop = int(positions[group[-1]] + window)
        if stop - start < minimum_duration:
            continue
        score = float(np.max(robust_z[group[0]:group[-1] + 1]))
        regions.append((max(0, start - guard), min(samples.size, stop + guard), score))
    return regions


@dataclass
class BurstDetection:
    start: int
    stop: int
    origin: int
    method: str
    score: float
    energy_z: float = float("nan")


def detect_bursts(
    compensated: np.ndarray,
    sample_rate: float,
    symbol_rate: float,
    *,
    z_threshold: float = 10.5,
    energy_z_threshold: float = 7.0,
    pretrigger_ms: float = 100.0,
    posttrigger_ms: float = 500.0,
) -> list[BurstDetection]:
    """Hybrid known-header correlation with spectral-energy fallback."""
    template = _header_template(sample_rate, symbol_rate)
    correlation = _normalised_match(compensated, template)
    standard_deviation = float(np.std(correlation)) if correlation.size else 0.0
    z = correlation / max(standard_deviation, 1e-12)
    minimum_gap = int(round(0.35 * sample_rate))
    before = int(round(pretrigger_ms * 1e-3 * sample_rate))
    after = int(round(posttrigger_ms * 1e-3 * sample_rate))

    ranked: list[tuple[int, float, str]] = []
    for polarity, values in (("HEADER", z), ("HEADER_INVERTED", -z)):
        indices = np.flatnonzero(values >= float(z_threshold))
        for index in indices:
            ranked.append((int(index), float(values[index]), polarity))
    ranked.sort(key=lambda item: item[1], reverse=True)

    detections: list[BurstDetection] = []
    for origin, score, method in ranked:
        if any(abs(origin - item.origin) <= minimum_gap for item in detections):
            continue
        detections.append(
            BurstDetection(
                start=max(0, origin - before),
                stop=min(compensated.size, origin + after),
                origin=origin,
                method=method,
                score=score,
            )
        )

    regions = _energy_regions(compensated, sample_rate, energy_z_threshold)
    for start, stop, energy_z in regions:
        overlapping = [
            item for item in detections
            if item.start < stop and start < item.stop
        ]
        if overlapping:
            for item in overlapping:
                if not np.isfinite(item.energy_z):
                    item.energy_z = energy_z
            continue
        detections.append(
            BurstDetection(
                start=start,
                stop=stop,
                origin=(start + stop) // 2,
                method="ENERGY_REGION",
                score=energy_z,
                energy_z=energy_z,
            )
        )

    detections.sort(key=lambda item: item.start)
    return detections


class OggHybridBurstScanner:
    """Buffer OGG audio, detect bursts, then reuse the common GFSK scanner."""

    def __init__(
        self,
        *args: object,
        header_z_threshold: float = 10.5,
        energy_z_threshold: float = 7.0,
        baseline_iterations: int = 7,
        baseline_half_window: int = 1000,
        detector_pretrigger_ms: float = 100.0,
        detector_posttrigger_ms: float = 500.0,
        **kwargs: object,
    ) -> None:
        kwargs["input_mode"] = "AUDIO_DISCRIMINATOR"
        kwargs["emit_decoded_candidate_files"] = False
        self.decoder = StreamingBurstScanner(*args, **kwargs)
        self.header_z_threshold = float(header_z_threshold)
        self.energy_z_threshold = float(energy_z_threshold)
        self.baseline_iterations = int(baseline_iterations)
        self.baseline_half_window = int(baseline_half_window)
        self.detector_pretrigger_ms = float(detector_pretrigger_ms)
        self.detector_posttrigger_ms = float(detector_posttrigger_ms)
        self._chunks: list[np.ndarray] = []
        self._finished = False

    def feed(self, samples: np.ndarray) -> None:
        if self._finished:
            raise RuntimeError("Cannot feed OGG samples after finish()")
        values = np.asarray(samples, dtype=np.float32).ravel()
        if values.size:
            self._chunks.append(values.copy())

    def finish(self) -> None:
        if self._finished:
            return
        self._finished = True
        if not self._chunks:
            return
        audio = np.concatenate(self._chunks).astype(np.float64, copy=False)
        self._chunks.clear()
        compensated = restore_baseline(
            audio,
            iterations=self.baseline_iterations,
            half_window=self.baseline_half_window,
        )
        detections = detect_bursts(
            compensated,
            self.decoder.sample_rate,
            self.decoder.symbol_rate,
            z_threshold=self.header_z_threshold,
            energy_z_threshold=self.energy_z_threshold,
            pretrigger_ms=self.detector_pretrigger_ms,
            posttrigger_ms=self.detector_posttrigger_ms,
        )
        print(
            "SCION-X OGG hybrid detector: {} candidate(s) [{}]".format(
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
            # Initialise an empty report so the pipeline can distinguish a
            # completed zero-detection run from a stale previous result.
            self.decoder._start_outputs()
            return

        quiet_mask = np.ones(compensated.size, dtype=bool)
        for item in detections:
            quiet_mask[item.start:item.stop] = False
        quiet = compensated[quiet_mask]
        noise_floor_db = (
            10.0 * np.log10(float(np.mean(quiet * quiet)) + 1e-20)
            if quiet.size
            else float("nan")
        )
        for item in detections:
            segment = compensated[item.start:item.stop].astype(np.complex64)
            self.decoder.process_explicit_burst(
                segment,
                item.start,
                noise_floor_db=noise_floor_db,
                force_extended_search=True,
                detector_method=item.method,
                detector_score=item.score,
                detector_energy_z=item.energy_z,
            )
        # Explicit bursts are finalized one-by-one above.  Calling finish is
        # still required so the shared decoder can run its strongest-N
        # no-preamble fallback searches and finalize the report.
        self.decoder.finish()
