"""
Sibilance-focused metrics — complements broadband fidelity scores.

Broadband IMD can miss HF splatter that sounds like distorted sibilance when one
transfer retains energy above 8 kHz and the other does not.
"""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy import signal

from audio_io import EDGE_TRIM_SEC

SIB_LO_HZ = 4000
SIB_HI_HZ = 8000


def hf_imd_proxy(x: np.ndarray, sr: int, lo: float = 3000, hi: float = 8000) -> dict[str, float]:
    """IMD products landing in the sibilance band from mid + HF partials."""
    win, hop = 16384, 8192
    vals = []
    trim = int(EDGE_TRIM_SEC * sr)
    xc = x[trim:-trim]
    for s in range(0, len(xc) - win, hop):
        seg = xc[s : s + win] * np.hanning(win)
        freqs = np.fft.rfftfreq(win, 1 / sr)
        power = np.abs(np.fft.rfft(seg)) ** 2 + 1e-20
        floor = float(np.median(power[(freqs >= 500) & (freqs <= 8000)]))
        lows, highs = [], []
        for band_lo, band_hi in ((400, 2500), (lo, hi)):
            m = (freqs >= band_lo) & (freqs <= band_hi)
            p, f = power[m], freqs[m]
            for i in range(2, len(p) - 2):
                if p[i] > p[i - 1] and p[i] > p[i + 1] and p[i] > floor * 4:
                    (lows if band_hi <= 2500 else highs).append(float(f[i]))
        lows = sorted(set(lows))[:4]
        highs = sorted(set(highs))[:4]
        imd = []
        for f1 in lows:
            for f2 in highs:
                for ft in (f2 + f1, abs(f2 - f1), f2 + 2 * f1):
                    if lo <= ft <= hi:
                        m = (freqs >= ft - 30) & (freqs <= ft + 30)
                        if np.any(m):
                            imd.append(10 * math.log10(float(np.max(power[m])) / floor))
        if imd:
            vals.append(max(imd))
    if not vals:
        return {}
    a = np.array(vals)
    return {"hf_imd_median": float(np.median(a)), "hf_imd_p95": float(np.percentile(a, 95))}


def sibilance_buzz_db(x: np.ndarray, sr: int) -> dict[str, float]:
    """Non-harmonic vs harmonic energy in 4–8 kHz on sibilant-heavy windows."""
    win, hop = 8192, 4096
    scores = []
    trim = int(EDGE_TRIM_SEC * sr)
    xc = x[trim:-trim]
    for s in range(0, len(xc) - win, hop):
        seg = xc[s : s + win] * np.hanning(win)
        freqs = np.fft.rfftfreq(win, 1 / sr)
        power = np.abs(np.fft.rfft(seg)) ** 2
        sib = (freqs >= SIB_LO_HZ) & (freqs <= SIB_HI_HZ)
        if np.sum(power[sib]) < np.sum(power) * 0.05:
            continue
        f0 = freqs[(freqs >= 300) & (freqs <= 3500)][
            np.argmax(power[(freqs >= 300) & (freqs <= 3500)])
        ]
        if f0 < 200:
            continue
        harm_e = inharm_e = 0.0
        for fk, pk in zip(freqs[sib], power[sib]):
            n = round(fk / f0)
            if n < 1:
                continue
            err = abs(fk - n * f0) / f0
            if err < 0.03:
                harm_e += pk
            elif err > 0.06:
                inharm_e += pk
        if harm_e > 0:
            scores.append(10 * math.log10(inharm_e / (harm_e + 1e-20)))
    if not scores:
        return {}
    a = np.array(scores)
    return {"sib_buzz_median": float(np.median(a)), "sib_buzz_p95": float(np.percentile(a, 95))}


def detect_sibilant_frames(x: np.ndarray, sr: int, win_sec: float = 0.05) -> list[int]:
    """Frame start indices where 4–8 kHz dominates (sibilant events)."""
    win = int(win_sec * sr)
    hop = win // 2
    trim = int(EDGE_TRIM_SEC * sr)
    frames = []
    for s in range(trim, len(x) - win - trim, hop):
        seg = x[s : s + win]
        f = np.fft.rfftfreq(win, 1 / sr)
        p = np.abs(np.fft.rfft(seg * np.hanning(win))) ** 2
        sib = float(np.sum(p[(f >= SIB_LO_HZ) & (f < SIB_HI_HZ)]))
        mid = float(np.sum(p[(f >= 500) & (f < 3000)])) + 1e-20
        if sib / mid > 0.15 and sib > np.sum(p) * 0.12:
            frames.append(s)
    return frames


def sibilant_hf_splatter(x: np.ndarray, sr: int, frames: list[int] | None = None) -> dict[str, float]:
    """
    During sibilant frames: energy in upper bands relative to 4–8 kHz core.
    Higher dB = more ultrahigh splatter on sibilants (often sounds 'distorted').
    """
    win = int(0.05 * sr)
    if frames is None:
        frames = detect_sibilant_frames(x, sr)
    if not frames:
        return {"sibilant_frame_count": 0}
    out = {"sibilant_frame_count": len(frames)}
    for lo, hi, key in [
        (6500, 8000, "sib_6p5_8k_vs_core_db"),
        (8000, 12000, "sib_8_12k_vs_core_db"),
        (12000, 20000, "sib_12_20k_vs_core_db"),
    ]:
        if hi > sr / 2:
            continue
        vals = []
        for s in frames:
            seg = x[s : s + win]
            f = np.fft.rfftfreq(win, 1 / sr)
            p = np.abs(np.fft.rfft(seg * np.hanning(win))) ** 2
            core = float(np.sum(p[(f >= SIB_LO_HZ) & (f < SIB_HI_HZ)])) + 1e-20
            hi_e = float(np.sum(p[(f >= lo) & (f < hi)]))
            vals.append(10 * math.log10(hi_e / core))
        out[f"{key}_median"] = float(np.median(vals))
        out[f"{key}_p95"] = float(np.percentile(vals, 95))
    return out


def analyze_sibilance(x: np.ndarray, sr: int) -> dict[str, float]:
    frames = detect_sibilant_frames(x, sr)
    return {
        **hf_imd_proxy(x, sr),
        **sibilance_buzz_db(x, sr),
        **sibilant_hf_splatter(x, sr, frames),
    }


def plot_sibilance(out_dir: Path, label_a: str, label_b: str, xa: np.ndarray, xb: np.ndarray, sr: int) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    for ax, (label, x), color in zip(axes, [(label_a, xa), (label_b, xb)], ["#1565c0", "#c62828"]):
        n = min(len(x), sr * 120)
        f, pxx = signal.welch(x[:n], sr, nperseg=16384)
        ax.semilogy(f, pxx + 1e-20, color=color, lw=1.2)
        ax.axvspan(SIB_LO_HZ, SIB_HI_HZ, alpha=0.12, color="orange")
        ax.set_xlim(2000, 16000)
        ax.set_xlabel("Hz")
        ax.set_ylabel("PSD")
        ax.set_title(label)
        ax.grid(alpha=0.2)
    fig.suptitle("HF spectrum — raw (no extra LP), first 120 s", fontsize=11)
    fig.tight_layout()
    fig.savefig(out_dir / "hf_psd_raw.png", dpi=160)
    plt.close(fig)

    trim = int(EDGE_TRIM_SEC * sr)
    seg_a, seg_b = xa[trim : trim + sr * 30], xb[trim : trim + sr * 30]
    n_fft, hop = 2048, 512

    def spec(x):
        chunks = []
        for i in range(0, len(x) - n_fft, hop):
            chunks.append(np.abs(np.fft.rfft(x[i : i + n_fft] * np.hanning(n_fft)))[: n_fft // 2])
        return 20 * np.log10(np.array(chunks).T + 1e-8)

    Sa, Sb = spec(seg_a), spec(seg_b)
    m = min(Sa.shape[1], Sb.shape[1])
    fig, axes = plt.subplots(3, 1, figsize=(14, 8), sharex=True)
    extent = [0, m * hop / sr, 0, sr / 2000]
    for ax, data, title in zip(
        axes, [Sa[:, :m], Sb[:, :m], Sa[:, :m] - Sb[:, :m]], [label_a, label_b, f"{label_a} − {label_b}"]
    ):
        im = ax.imshow(data, aspect="auto", origin="lower", extent=extent, cmap="magma", vmin=-70, vmax=-25)
        ax.set_ylabel("kHz")
        ax.set_ylim(2, 12)
        ax.set_title(title, fontsize=9)
        fig.colorbar(im, ax=ax, fraction=0.02)
    axes[-1].set_xlabel("Time (s)")
    fig.suptitle("30 s after trim — HF focus (2–12 kHz)")
    fig.tight_layout()
    fig.savefig(out_dir / "hf_spectrogram_30s.png", dpi=140)
    plt.close(fig)
