"""
Sliding-window fidelity metrics for optical soundtrack transfers.

Focus: cross-modulation, inharmonic distortion, HF harshness, clicks, AM complexity.
Deliberately de-emphasises the 24 Hz frame imprint (handled elsewhere).
"""

from __future__ import annotations

import math

import matplotlib.pyplot as plt
import numpy as np

from audio_io import EDGE_TRIM_SEC

WINDOW_SEC = 2.0
HOP_SEC = 1.0
MOD_EXCLUDE_HZ = (22.0, 26.0)


def _trim_stats(t: np.ndarray, v: np.ndarray) -> dict[str, float]:
    m = (t >= EDGE_TRIM_SEC) & (t <= t.max() - EDGE_TRIM_SEC)
    if not np.any(m):
        m = np.ones_like(t, dtype=bool)
    a = v[np.isfinite(v)]
    b = v[m & np.isfinite(v)]
    if len(b) == 0:
        b = a
    return {
        "median_full": float(np.nanmedian(a)),
        "p95_full": float(np.nanpercentile(a, 95)),
        "max_full": float(np.nanmax(a)),
        "median_trim": float(np.nanmedian(b)),
        "p95_trim": float(np.nanpercentile(b, 95)),
        "max_trim": float(np.nanmax(b)),
    }


def _find_peaks(freqs: np.ndarray, power: np.ndarray, min_hz: float, max_hz: float, top_n: int = 6):
    b = (freqs >= min_hz) & (freqs <= max_hz)
    f, p = freqs[b], power[b]
    if len(p) < 5:
        return []
    floor = float(np.median(p))
    thresh = floor * 10 ** (6 / 10)
    peaks = []
    for i in range(2, len(p) - 2):
        if p[i] >= thresh and p[i] >= p[i - 1] and p[i] >= p[i + 1] and p[i] >= p[i - 2] and p[i] >= p[i + 2]:
            peaks.append((float(f[i]), float(p[i])))
    peaks.sort(key=lambda x: x[1], reverse=True)
    return peaks[:top_n]


def _band_energy_db(freqs: np.ndarray, power: np.ndarray, f0: float, bw: float, ref_floor: float) -> float:
    m = (freqs >= f0 - bw) & (freqs <= f0 + bw)
    if not np.any(m):
        return float("nan")
    return 10 * math.log10(float(np.max(power[m])) / ref_floor)


def imd_crossmod_proxy(freqs: np.ndarray, power: np.ndarray) -> float:
    """Energy at sum/difference products of strong low + high partials (dB above floor)."""
    floor = float(np.median(power[(freqs >= 100) & (freqs <= 12000)])) or 1e-20
    low = _find_peaks(freqs, power, 120, 2500, top_n=5)
    high = _find_peaks(freqs, power, 2500, min(12000, freqs[-1] * 0.45), top_n=4)
    if not low or not high:
        return float("nan")
    imd_vals = []
    for f1, _ in low:
        for f2, _ in high:
            for f_tgt in (f2 + f1, abs(f2 - f1), f2 + 2 * f1, abs(f2 - 2 * f1)):
                if f_tgt < freqs[-1] * 0.98:
                    imd_vals.append(_band_energy_db(freqs, power, f_tgt, bw=max(40, f_tgt * 0.02), ref_floor=floor))
    return float(np.nanmax(imd_vals)) if imd_vals else float("nan")


def inharmonic_residual_db(freqs: np.ndarray, power: np.ndarray) -> float:
    """Mean level of spectral peaks that do not lock to harmonics of the strongest partial."""
    band = (freqs >= 150) & (freqs <= 6000)
    f, p = freqs[band], power[band]
    if len(p) < 20:
        return float("nan")
    floor = float(np.median(p))
    peaks = _find_peaks(f, p, 150, 6000, top_n=12)
    if len(peaks) < 2:
        return float("nan")
    f0 = peaks[0][0]
    if f0 < 80:
        return float("nan")
    inharm = []
    for fk, pk in peaks[1:]:
        n = round(fk / f0)
        if n < 1:
            continue
        if abs(fk - n * f0) / f0 * 100 > 4.0:
            inharm.append(10 * math.log10(pk / floor))
    return float(np.nanmean(inharm)) if inharm else 0.0


def hf_harshness_db(freqs: np.ndarray, power: np.ndarray) -> float:
    """6–14 kHz energy vs 0.3–2.5 kHz. Higher (less negative) = harsher sibilance/HF."""
    def band_pwr(lo, hi):
        m = (freqs >= lo) & (freqs <= hi)
        return float(np.sum(power[m])) if np.any(m) else 1e-20

    return 10 * math.log10(band_pwr(6000, min(14000, freqs[-1] * 0.95)) / band_pwr(300, 2500))


def spectral_flatness(freqs: np.ndarray, power: np.ndarray) -> float:
    """0.4–8 kHz flatness. Higher = noisier / less tonal."""
    m = (freqs >= 400) & (freqs <= 8000)
    p = power[m]
    p = p[p > 0]
    if len(p) < 10:
        return float("nan")
    return float(np.exp(np.mean(np.log(p))) / np.mean(p))


def noise_floor_db(freqs: np.ndarray, power: np.ndarray) -> float:
    """p90/p10 spectral ratio in 0.5–8 kHz. Higher = better SNR proxy."""
    m = (freqs >= 500) & (freqs <= 8000)
    p = power[m]
    if len(p) < 10:
        return float("nan")
    return 10 * math.log10(float(np.percentile(p, 90)) / (float(np.percentile(p, 10)) + 1e-20))


def impulse_index(seg: np.ndarray) -> float:
    """Click-like sample jumps per 10k samples."""
    d = np.abs(np.diff(seg))
    if len(d) == 0:
        return 0.0
    thr = max(float(np.median(d)) * 8.0, float(np.percentile(d, 99.5)))
    return float(np.mean(d > thr)) * 1e4


def modulation_complexity_db(seg: np.ndarray, sr: int) -> float:
    """Broadband AM sideband energy 2–200 Hz, excluding 22–26 Hz frame band."""
    sm = max(16, int(sr * 0.008))
    env = np.convolve(np.abs(seg), np.ones(sm, dtype=np.float32) / sm, mode="same")
    env = env - np.mean(env)
    hann = np.hanning(len(env)).astype(np.float32)
    E = np.abs(np.fft.rfft(env * hann)) ** 2 + 1e-20
    ef = np.fft.rfftfreq(len(env), 1.0 / sr)
    use = (ef >= 2) & (ef <= 200) & ~((ef >= MOD_EXCLUDE_HZ[0]) & (ef <= MOD_EXCLUDE_HZ[1]))
    if not np.any(use):
        return float("nan")
    floor = float(np.median(E[(ef >= 2) & (ef <= 30)]))
    return 10 * math.log10(float(np.sum(E[use])) / (floor * np.sum(use) + 1e-20))


def aliasing_proxy_db(freqs: np.ndarray, power: np.ndarray) -> float:
    """Very-high vs mid-band spectral energy."""
    nyq = freqs[-1]
    top = (freqs >= nyq * 0.85) & (freqs <= nyq * 0.98)
    mid = (freqs >= 1000) & (freqs <= 4000)
    if not np.any(top) or not np.any(mid):
        return float("nan")
    return 10 * math.log10(float(np.sum(power[top])) / (float(np.sum(power[mid])) + 1e-20))


def analyze_window(seg: np.ndarray, sr: int) -> dict[str, float]:
    win = len(seg)
    hann = np.hanning(win).astype(np.float32)
    freqs = np.fft.rfftfreq(win, 1.0 / sr)
    power = (np.abs(np.fft.rfft(seg * hann)) ** 2) + 1e-20
    rms = float(np.sqrt(np.mean(seg**2)) + 1e-20)
    return {
        "imd_db": imd_crossmod_proxy(freqs, power),
        "inharm_db": inharmonic_residual_db(freqs, power),
        "hf_harsh_db": hf_harshness_db(freqs, power),
        "flatness": spectral_flatness(freqs, power),
        "noise_floor_db": noise_floor_db(freqs, power),
        "impulse_idx": impulse_index(seg),
        "crest_factor": float(np.max(np.abs(seg)) / rms),
        "mod_complex_db": modulation_complexity_db(seg, sr),
        "alias_db": aliasing_proxy_db(freqs, power),
    }


def analyze_file(x: np.ndarray, sr: int) -> dict:
    """Slide 2 s windows / 1 s hop; return time series, per-metric stats, composite stress."""
    win = int(WINDOW_SEC * sr)
    hop = int(HOP_SEC * sr)
    if len(x) < win:
        return {}
    keys = list(analyze_window(x[:win], sr).keys())
    series = {k: [] for k in keys}
    times = []
    for s in range(0, len(x) - win + 1, hop):
        m = analyze_window(x[s : s + win], sr)
        times.append((s + win / 2) / sr)
        for k in keys:
            series[k].append(m[k])
    t = np.array(times)
    stats = {k: _trim_stats(t, np.array(series[k], dtype=float)) for k in keys}

    def p95(k):
        return stats[k]["p95_trim"]

    composite = (
        0.30 * (p95("imd_db") if np.isfinite(p95("imd_db")) else 0)
        + 0.20 * (p95("inharm_db") if np.isfinite(p95("inharm_db")) else 0)
        + 0.15 * (p95("hf_harsh_db") if np.isfinite(p95("hf_harsh_db")) else 0)
        + 0.10 * (p95("mod_complex_db") if np.isfinite(p95("mod_complex_db")) else 0)
        + 0.10 * (p95("noise_floor_db") if np.isfinite(p95("noise_floor_db")) else 0)
        + 0.05 * (p95("alias_db") if np.isfinite(p95("alias_db")) else 0)
        + 0.05 * (p95("impulse_idx") if np.isfinite(p95("impulse_idx")) else 0)
        + 0.05 * ((p95("flatness") - 0.1) * 40 if np.isfinite(p95("flatness")) else 0)
    )
    return {
        "duration_sec": float(len(x) / sr),
        "sr": sr,
        "windows": len(t),
        "time_s": t.tolist(),
        "series": series,
        "stats": stats,
        "composite_stress_trim": float(composite),
    }


def summarize_fidelity(result: dict) -> dict[str, float]:
    out = {"composite_stress_trim": result["composite_stress_trim"]}
    for k, st in result["stats"].items():
        out[f"{k}_median_trim"] = st["median_trim"]
        out[f"{k}_p95_trim"] = st["p95_trim"]
    return out


def plot_timelines(stem: str, result: dict, out_path) -> None:
    t = np.array(result["time_s"])
    fig, axes = plt.subplots(4, 1, figsize=(14, 10), sharex=True)
    panels = [
        ("imd_db", "IMD / cross-mod proxy (dB)", "tab:purple"),
        ("inharm_db", "Inharmonic residual (dB)", "tab:orange"),
        ("hf_harsh_db", "HF harshness (dB)", "tab:red"),
        ("mod_complex_db", "AM complexity excl. 24 Hz (dB)", "tab:green"),
    ]
    for ax, (key, title, col) in zip(axes, panels):
        v = np.array(result["series"][key], dtype=float)
        ax.plot(t, v, color=col, lw=0.9)
        ax.axhline(result["stats"][key]["median_trim"], ls="--", color="gray", lw=0.8)
        ax.set_title(title, fontsize=9)
        ax.grid(alpha=0.2)
    axes[-1].set_xlabel("Time (s)")
    fig.suptitle(f"Audio fidelity — {stem}", fontsize=11)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
