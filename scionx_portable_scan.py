"""File-based CS16/OGG entry point without GNU Radio or conda.

Shares the hybrid detectors, candidate scoring and soft repair with the GRC
front ends. MATLAB remains responsible for formal HDLC deframing and TLM.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


def load_settings(filename: Path, mode: str) -> dict:
    data = json.loads(filename.read_text(encoding="utf-8-sig"))
    if data.get("schema_version") != 1:
        raise ValueError("Unsupported decoder_settings.json schema_version")
    settings = dict(data["common"])
    settings.update(data[mode])
    return settings


def check_environment() -> None:
    import numpy
    import scipy
    import soundfile
    from scionx_burst_scan import _load_advanced_decoder
    from scionx_cs16_burst_detect import Cs16HybridBurstScanner
    from scionx_ogg_burst_detect import OggHybridBurstScanner

    _load_advanced_decoder()
    if "OGG" not in soundfile.available_formats():
        raise RuntimeError("SoundFile/libsndfile does not support OGG")
    if any(name == "gnuradio" or name.startswith("gnuradio.") for name in sys.modules):
        raise RuntimeError("Portable front end unexpectedly imported GNU Radio")
    root = Path(__file__).resolve().parent
    # Instantiate without touching output files, so invalid settings fail early.
    for mode, cls in (("cs16", Cs16HybridBurstScanner), ("ogg", OggHybridBurstScanner)):
        cls(output_txt="unused.txt", report_csv="unused.csv",
            **load_settings(root / "decoder_settings.json", mode))
    print(json.dumps({"python": sys.executable, "version": sys.version.split()[0],
                      "numpy": numpy.__version__, "scipy": scipy.__version__,
                      "soundfile": soundfile.__version__,
                      "libsndfile": soundfile.__libsndfile_version__,
                      "soft_repair": "AVAILABLE", "gnuradio_imported": False}))


def main() -> int:
    root = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--iq-file", type=Path)
    group.add_argument("--ogg-file", type=Path)
    parser.add_argument("--output-txt-path", default="Frames/GNU_GFSK_Raw_HEX.txt")
    parser.add_argument("--report-csv-path", default="GNU_Burst_Scan_Report.csv")
    parser.add_argument("--settings", type=Path, default=root / "decoder_settings.json")
    parser.add_argument("--check-environment", action="store_true")
    args = parser.parse_args()
    if args.check_environment:
        check_environment()
        return 0
    source = args.iq_file or args.ogg_file
    if source is None:
        parser.error("one of --iq-file or --ogg-file is required")
    if not source.is_file():
        raise FileNotFoundError(f"Recording not found: {source}")
    mode = "cs16" if args.iq_file else "ogg"
    settings = load_settings(args.settings, mode)
    for key in ("output_txt", "report_csv", "input_file", "input_mode", "emit_decoded_candidate_files"):
        if key in settings:
            raise ValueError(f"Setting {key} is controlled by the pipeline, not the settings file")
    outputs = dict(output_txt=args.output_txt_path, report_csv=args.report_csv_path)
    print(f"SCION-X portable Python: {mode.upper()} | {source}", flush=True)
    if mode == "cs16":
        from scionx_cs16_burst_detect import Cs16HybridBurstScanner
        if source.stat().st_size == 0 or source.stat().st_size % 4:
            raise ValueError("CS16 must contain nonempty interleaved little-endian int16 I/Q pairs (4 bytes/sample)")
        # The shared CS16 scanner reads the original file via bounded memmap chunks.
        scanner = Cs16HybridBurstScanner(input_file=str(source), **outputs, **settings)
        scanner.finish()
    else:
        import soundfile as sf
        from scionx_ogg_burst_detect import OggHybridBurstScanner
        with sf.SoundFile(str(source), mode="r") as audio:
            if audio.channels != 1:
                raise ValueError(f"OGG must be mono discriminator audio; found {audio.channels} channels")
            expected = int(settings["sample_rate"])
            if audio.samplerate != expected:
                raise ValueError(f"OGG sample rate is {audio.samplerate} Hz, expected {expected} Hz")
            if len(audio) == 0:
                raise ValueError("OGG contains no audio samples")
            scanner = OggHybridBurstScanner(**outputs, **settings)
            for block in audio.blocks(blocksize=65536, dtype="float32", always_2d=False):
                scanner.feed(block)
        scanner.finish()
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"SCION-X portable decoder failed: {type(exc).__name__}: {exc}", file=sys.stderr)
        raise SystemExit(1)
