#!/usr/bin/env python3
"""
A/B comparison of two digitizations of the same optical soundtrack element.

Applies time alignment, optional 8 kHz low-pass (1960s cinema bandwidth), and RMS
gain matching before computing restoration-oriented metrics.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal

from audio_io import EDGE_TRIM_SEC, ensure_readable_wav, read_wav
from fidelity_metrics import analyze_file, plot_timelines, summarize_fidelity
from sibilance_metrics import analyze_sibilance, plot_sibilance

LP_CUTOFF_HZ = 8000.0
LP_ORDER = 4
ALIGN_DS_HZ = 6000
EDGE_TRIM = EDGE_TRIM_SEC


def lowpass(x: np.ndarray, sr: int, cutoff: float = LP_CUTOFF_HZ, order: int = LP_ORDER) -> np.ndarray:
    wn = min(cutoff / (sr / 2), 0.999)
    sos = signal.butter(order, wn, btype="low", output="sos")
    return signal.sosfiltfilt(sos, x).astype(np.float32)


def align_signals(xa: np.ndarray, xb: np.ndarray, sr: int) -> tuple[np.ndarray, np.ndarray, int]:
    """Cross-correlate band-limited downsampled audio; return overlap and lag in samples (B vs A)."""
    ds = ALIGN_DS_HZ
    da = np.interp(np.linspace(0, len(xa), int(len(xa) * ds / sr)), np.arange(len(xa)), xa)
    db = np.interp(np.linspace(0, len(xb), int(len(xb) * ds / sr)), np.arange(len(xb)), xb)
    sos = signal.butter(2, [100 / (ds / 2), min(2800, ds * 0.45) / (ds / 2)], btype="band", output="sos")
    da = signal.sosfiltfilt(sos, da)
    db = signal.sosfiltfilt(sos, db)
    n = min(len(da), len(db))
    seg = int(min(n, ds * 120))
    start = (n - seg) // 2
    a = (da[start : start + seg] - np.mean(da[start : start + seg])) / (np.std(da[start : start + seg]) + 1e-12)
    b = (db[start : start + seg] - np.mean(db[start : start + seg])) / (np.std(db[start : start + seg]) + 1e-12)
    corr = signal.correlate(a, b, mode="full", method="fft")
    lags = np.arange(-len(a) + 1, len(a))
    lag = int(round(lags[np.argmax(corr)] * sr / ds))
    a0, b0 = (0, lag) if lag >= 0 else (-lag, 0)
    n_overlap = min(len(xa) - a0, len(xb) - b0)
    return xa[a0 : a0 + n_overlap], xb[b0 : b0 + n_overlap], lag


def rms(x: np.ndarray) -> float:
    return float(np.sqrt(np.mean(x**2)) + 1e-20)


def match_rms_sr(xa: np.ndarray, xb: np.ndarray, sr: int, trim_sec: float = EDGE_TRIM) -> tuple[np.ndarray, np.ndarray, float, float]:
    t = int(trim_sec * sr)
    if len(xa) <= 2 * t:
        gb = rms(xa) / rms(xb)
        return xa, xb * gb, 1.0, gb
    mid = slice(t, len(xa) - t)
    gb = rms(xa[mid]) / rms(xb[mid])
    return xa, xb * gb, 1.0, gb


def envelope_stats(x: np.ndarray, sr: int) -> dict[str, float]:
    analytic = signal.hilbert(x)
    env_db = 20 * np.log10(np.abs(analytic) + 1e-12)
    trim = int(EDGE_TRIM * sr)
    e = env_db[trim:-trim] if len(env_db) > 2 * trim else env_db
    return {
        "mod_depth": float(np.std(analytic) / (np.mean(np.abs(analytic)) + 1e-12)),
        "dynamic_range_db": float(np.percentile(e, 99) - np.percentile(e, 10)),
    }


def rumble_db(x: np.ndarray, sr: int) -> float:
    f, pxx = signal.welch(x[: min(len(x), sr * 4)], sr, nperseg=8192)

    def band(lo, hi):
        m = (f >= lo) & (f < hi)
        return float(np.sum(pxx[m])) + 1e-20

    return 10 * math.log10(band(20, 80) / band(300, 2000))


def flutter_proxy_db(x: np.ndarray, sr: int) -> float:
    win = 16384
    if len(x) < win:
        return float("nan")
    freqs = np.fft.rfftfreq(win, 1.0 / sr)
    power = np.abs(np.fft.rfft(x[:win] * np.hanning(win))) ** 2 + 1e-20
    band = (freqs >= 300) & (freqs <= 2500)
    i0 = int(np.argmax(power[band]))
    f0 = float(freqs[band][i0])
    fund = float(power[band][i0])
    side = 0.0
    for d in (2, 3, 4, 6, 8, 12, 24):
        for sign in (-1, 1):
            ft = f0 + sign * d
            m = (freqs >= ft - 1.5) & (freqs <= ft + 1.5)
            if np.any(m):
                side += float(np.max(power[m]))
    return 10 * math.log10(side / fund) if fund > 0 else float("nan")


def snr_modulation_db(x: np.ndarray, sr: int) -> float:
    env = np.abs(signal.hilbert(x))
    sm = max(16, int(sr * 0.02))
    env_s = np.convolve(env, np.ones(sm) / sm, mode="same")
    noise = x - env_s * np.sign(x) * (np.abs(x) / (env + 1e-12))
    trim = int(EDGE_TRIM * sr)
    sig = env_s[trim:-trim] if len(env_s) > 2 * trim else env_s
    noi = noise[trim:-trim] if len(noise) > 2 * trim else noise
    return 20 * math.log10(rms(sig) / (rms(noi) + 1e-20))


def spectral_profile(x: np.ndarray, sr: int) -> dict[str, float]:
    f, pxx = signal.welch(x[: min(len(x), sr * 30)], sr, nperseg=16384)
    total = float(np.sum(pxx)) + 1e-20

    def pct(lo, hi):
        return 100 * float(np.sum(pxx[(f >= lo) & (f < hi)])) / total

    cum = np.cumsum(pxx) / total
    return {
        "pct_above_8k_hz": pct(8000, sr / 2),
        **{f"rolloff_{p}pct_hz": float(f[np.searchsorted(cum, p / 100)]) for p in (85, 95, 99)},
    }


def aggregate_extra_metrics(x: np.ndarray, sr: int) -> dict[str, float]:
    return {
        **envelope_stats(x, sr),
        "rumble_db": rumble_db(x, sr),
        "flutter_proxy_db": flutter_proxy_db(x, sr),
        "snr_modulation_db": snr_modulation_db(x, sr),
        "rms": rms(x),
        "peak": float(np.max(np.abs(x))),
        "crest_factor": float(np.max(np.abs(x)) / rms(x)),
        **spectral_profile(x, sr),
    }


def plot_comparison(out_dir: Path, label_a: str, label_b: str, metrics_a: dict, metrics_b: dict, xa: np.ndarray, xb: np.ndarray, sr: int) -> None:
    panels = [
        ("composite_stress_trim", "Composite stress", True),
        ("imd_db_p95_trim", "IMD / cross-mod", True),
        ("hf_harsh_db_p95_trim", "HF harshness", True),
        ("noise_floor_db_p95_trim", "Spectral SNR", False),
        ("rumble_db", "Rumble", True),
        ("sib_8_12k_vs_core_db_median", "Sib splatter 8–12 kHz", True),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(12, 7))
    for ax, (key, title, lower_better) in zip(axes.flat, panels):
        va, vb = metrics_a.get(key, float("nan")), metrics_b.get(key, float("nan"))
        bars = ax.bar([0, 1], [va, vb], color=["#1565c0", "#c62828"], width=0.6)
        ax.set_xticks([0, 1])
        ax.set_xticklabels([label_a.split()[0], label_b.split()[0]], fontsize=8)
        ax.set_title(title, fontsize=9)
        if np.isfinite(va) and np.isfinite(vb):
            w = 0 if (va < vb) == lower_better else 1
            bars[w].set_edgecolor("#2e7d32")
            bars[w].set_linewidth(2.5)
        ax.grid(axis="y", alpha=0.2)
    fig.suptitle("Restoration metrics (8 kHz LP, gain-matched)", fontsize=11)
    fig.tight_layout()
    fig.savefig(out_dir / "metrics_comparison.png", dpi=160)
    plt.close(fig)

    n = min(len(xa), len(xb), sr * 60)
    f, pa = signal.welch(xa[:n], sr, nperseg=16384)
    _, pb = signal.welch(xb[:n], sr, nperseg=16384)
    fig, ax = plt.subplots(figsize=(12, 4))
    ax.semilogy(f, pa + 1e-20, label=label_a, lw=1.2)
    ax.semilogy(f, pb + 1e-20, label=label_b, lw=1.2, alpha=0.85)
    ax.axvline(LP_CUTOFF_HZ, color="gray", ls="--", lw=0.8)
    ax.set_xlim(20, 12000)
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("PSD")
    ax.set_title("Welch PSD after 8 kHz low-pass")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(out_dir / "spectrum_overlay.png", dpi=160)
    plt.close(fig)


def score_restoration(metrics: dict) -> float:
    weights = {
        "composite_stress_trim": 0.20,
        "imd_db_p95_trim": 0.12,
        "hf_harsh_db_p95_trim": 0.12,
        "noise_floor_db_p95_trim": -0.10,
        "snr_modulation_db": -0.08,
        "rumble_db": 0.10,
        "flutter_proxy_db": 0.08,
        "sib_8_12k_vs_core_db_median": 0.10,
        "impulse_idx_p95_trim": 0.05,
    }
    score = 0.0
    for k, w in weights.items():
        v = metrics.get(k)
        if v is not None and np.isfinite(v):
            score += w * v
    return score


def run_compare(
    path_a: Path,
    label_a: str,
    path_b: Path,
    label_b: str,
    out_dir: Path,
    apply_lp: bool = True,
    match_gain: bool = True,
) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    tmp_dir = out_dir / "_tmp_pcm"

    xa, sr_a = read_wav(ensure_readable_wav(path_a, tmp_dir))
    xb, sr_b = read_wav(ensure_readable_wav(path_b, tmp_dir))
    if sr_a != sr_b:
        raise ValueError(f"Sample rate mismatch: {sr_a} vs {sr_b}")
    sr = sr_a

    xa_al, xb_al, lag = align_signals(xa, xb, sr)
    raw_a, raw_b = aggregate_extra_metrics(xa_al, sr), aggregate_extra_metrics(xb_al, sr)

    xa_lp = lowpass(xa_al, sr) if apply_lp else xa_al.copy()
    xb_lp = lowpass(xb_al, sr) if apply_lp else xb_al.copy()
    ga, gb = 1.0, 1.0
    if match_gain:
        xa_lp, xb_lp, ga, gb = match_rms_sr(xa_lp, xb_lp, sr)

    fid_a = summarize_fidelity(analyze_file(xa_lp, sr))
    fid_b = summarize_fidelity(analyze_file(xb_lp, sr))
    sib_a = analyze_sibilance(xa_al, sr)  # raw — sibilance splatter is pre-LP phenomenon
    sib_b = analyze_sibilance(xb_al, sr)
    extra_a = aggregate_extra_metrics(xa_lp, sr)
    extra_b = aggregate_extra_metrics(xb_lp, sr)
    metrics_a = {**fid_a, **extra_a, **sib_a}
    metrics_b = {**fid_b, **extra_b, **sib_b}

    score_a, score_b = score_restoration(metrics_a), score_restoration(metrics_b)
    winner = label_a if score_a < score_b else label_b

    plot_comparison(out_dir, label_a, label_b, metrics_a, metrics_b, xa_lp, xb_lp, sr)
    plot_sibilance(out_dir / "sibilance", label_a, label_b, xa_al, xb_al, sr)

    ra, rb = analyze_file(xa_lp, sr), analyze_file(xb_lp, sr)
    plot_timelines(label_a, ra, out_dir / f"{path_a.stem}_fidelity_timeline.png")
    plot_timelines(label_b, rb, out_dir / f"{path_b.stem}_fidelity_timeline.png")

    fig, axes = plt.subplots(2, 2, figsize=(11, 7))
    for ax, (key, title) in zip(
        axes.flat,
        [("imd_db", "IMD"), ("inharm_db", "Inharmonic"), ("hf_harsh_db", "HF harshness"), ("noise_floor_db", "SNR proxy")],
    ):
        ta, tb = np.array(ra["time_s"]), np.array(rb["time_s"])
        ax.plot(ta / ta.max(), ra["series"][key], label=label_a, lw=0.9)
        ax.plot(tb / tb.max(), rb["series"][key], label=label_b, lw=0.9, alpha=0.85)
        ax.set_title(title)
        ax.legend(fontsize=7)
        ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(out_dir / "fidelity_timeline_AB.png", dpi=160)
    plt.close(fig)

    rows = [
        ("Composite stress (↓)", "composite_stress_trim", True),
        ("IMD p95 (↓)", "imd_db_p95_trim", True),
        ("HF harshness p95 (↓)", "hf_harsh_db_p95_trim", True),
        ("Spectral SNR p95 (↑)", "noise_floor_db_p95_trim", False),
        ("Rumble 20–80 Hz (↓)", "rumble_db", True),
        ("Flutter proxy (↓)", "flutter_proxy_db", True),
        ("Sib splatter 8–12 kHz (↓)", "sib_8_12k_vs_core_db_median", True),
        ("Sib splatter 12–20 kHz (↓)", "sib_12_20k_vs_core_db_median", True),
        ("% energy >8 kHz (↓)", "pct_above_8k_hz", True),
    ]
    md = [
        f"# Restoration digitization A/B: {label_a} vs {label_b}",
        "",
        "## Preprocessing",
        "",
        f"| Step | Value |",
        f"|---|---|",
        f"| Alignment lag | {lag} samples ({lag/sr:.3f} s), B relative to A |",
        f"| Overlap | {len(xa_al)/sr:.1f} s |",
        f"| 8 kHz low-pass | {'Butterworth order ' + str(LP_ORDER) if apply_lp else 'off'} |",
        f"| Gain scale B | ×{gb:.4f} (A reference) |",
        "",
        "## Raw levels (aligned, before gain match)",
        "",
        f"| | RMS | Peak | Crest | >8 kHz % |",
        f"|---|---:|---:|---:|---:|",
        f"| {label_a} | {raw_a['rms']:.5f} | {raw_a['peak']:.4f} | {raw_a['crest_factor']:.1f} | {raw_a['pct_above_8k_hz']:.3f} |",
        f"| {label_b} | {raw_b['rms']:.5f} | {raw_b['peak']:.4f} | {raw_b['crest_factor']:.1f} | {raw_b['pct_above_8k_hz']:.3f} |",
        "",
        "## Metrics (after LP + gain match, except sibilance splatter = raw)",
        "",
        "| Metric | " + label_a + " | " + label_b + " | Better |",
        "|---|---:|---:|---|",
    ]
    for title, key, lower in rows:
        va, vb = metrics_a.get(key, float("nan")), metrics_b.get(key, float("nan"))
        better = label_a if (va < vb) == lower else label_b
        md.append(f"| {title} | {va:.3f} | {vb:.3f} | **{better}** |")

    md.extend([
        "",
        "## Ranking",
        "",
        f"- {label_a}: **{score_a:.2f}**",
        f"- {label_b}: **{score_b:.2f}**",
        f"- Suggested starting point: **{winner}** (lower score = cleaner)",
        "",
        "See `README.md` in the `analysis/` folder for metric definitions.",
        "",
        "## Outputs",
        "",
        "- `metrics_comparison.png`, `spectrum_overlay.png`",
        "- `fidelity_timeline_AB.png`, per-file timelines",
        "- `sibilance/hf_psd_raw.png`, `sibilance/hf_spectrogram_30s.png`",
        "- `comparison_summary.json`",
    ])

    summary = {
        "file_a": str(path_a),
        "file_b": str(path_b),
        "preprocessing": {"lag_samples": lag, "gain_scale_b": gb, "lp_hz": LP_CUTOFF_HZ if apply_lp else None},
        "raw_levels": {"a": raw_a, "b": raw_b},
        "metrics": {"a": metrics_a, "b": metrics_b},
        "scores": {"a": score_a, "b": score_b, "suggested": winner},
    }
    (out_dir / "comparison_summary.json").write_text(json.dumps(summary, indent=2, default=float))
    (out_dir / "restoration_AB_report.md").write_text("\n".join(md) + "\n")
    print(f"Report → {out_dir / 'restoration_AB_report.md'}")
    print(f"Suggested: {winner} ({score_a:.2f} vs {score_b:.2f})")
    return summary


def main():
    p = argparse.ArgumentParser(description="Compare two optical soundtrack digitizations")
    p.add_argument("file_a", type=Path)
    p.add_argument("file_b", type=Path)
    p.add_argument("--label-a", default=None)
    p.add_argument("--label-b", default=None)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--no-lp", action="store_true")
    p.add_argument("--no-gain-match", action="store_true")
    args = p.parse_args()
    run_compare(
        args.file_a,
        args.label_a or args.file_a.stem,
        args.file_b,
        args.label_b or args.file_b.stem,
        args.out,
        apply_lp=not args.no_lp,
        match_gain=not args.no_gain_match,
    )


if __name__ == "__main__":
    main()
