"""WAV loading helpers. Converts float/IEEE WAV via ffmpeg when the stdlib wave module cannot read them."""

from __future__ import annotations

import subprocess
import wave
from pathlib import Path

import numpy as np

EDGE_TRIM_SEC = 10.0


def read_wav(path: Path) -> tuple[np.ndarray, int]:
    """Load mono float32 audio in range ~[-1, 1], DC-removed."""
    with wave.open(str(path), "rb") as w:
        ch, sw, sr, n = w.getnchannels(), w.getsampwidth(), w.getframerate(), w.getnframes()
        raw = w.readframes(n)
    if sw == 3:
        b = np.frombuffer(raw, dtype=np.uint8).reshape(-1, 3)
        vals = b[:, 0].astype(np.int32) | (b[:, 1].astype(np.int32) << 8) | (b[:, 2].astype(np.int32) << 16)
        vals = vals - (1 << 24) * ((vals & 0x800000) > 0)
        x = vals.astype(np.float32) / (1 << 23)
    elif sw == 2:
        x = np.frombuffer(raw, dtype=np.int16).astype(np.float32) / 32768.0
    elif sw == 4:
        x = np.frombuffer(raw, dtype=np.int32).astype(np.float32) / (1 << 31)
    elif sw == 1:
        x = (np.frombuffer(raw, dtype=np.uint8).astype(np.float32) - 128.0) / 128.0
    else:
        raise RuntimeError(f"Unsupported sample width {sw} for {path}")
    if ch > 1:
        x = x.reshape(-1, ch).mean(axis=1)
    return x - np.mean(x), sr


def ensure_readable_wav(path: Path, tmp_dir: Path) -> Path:
    """Return a PCM WAV path readable by the stdlib wave module."""
    try:
        with wave.open(str(path), "rb") as w:
            w.getsampwidth()
        return path
    except wave.Error:
        pass
    out = tmp_dir / f"{path.stem}_pcm24.wav"
    tmp_dir.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        ["ffmpeg", "-y", "-i", str(path), "-c:a", "pcm_s24le", str(out)],
        check=True,
        capture_output=True,
    )
    return out
