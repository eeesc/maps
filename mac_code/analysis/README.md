# Optical soundtrack digitization comparison

Objective A/B analysis for two transfers of the **same** optical sound element (e.g. duplicate positive, variable area). Designed to answer: *which digitization is a better starting point for digital restoration?*

Metrics cover cross-modulation, SNR, contrast/modulation, rumble, flutter, clicks, and **sibilance-specific HF splatter** — the last of which broadband scores can miss.

---

## Quick start

```bash
cd /path/to/digitizations_tests

python3 analysis/compare.py \
  PEKARUV_CISAR_DNZ6394_LR01.wav \
  DUPNEGSPINNER_PC_01_2026-06-29_13-45-49_resampled.wav \
  --label-a "Sondor PEKARUV LR01" \
  --label-b "ITST DupNeg Spinner" \
  --out comparison_LR01
```

**Requirements:** Python 3.10+, `numpy`, `scipy`, `matplotlib`. `ffmpeg` on PATH for IEEE-float WAV files.

```bash
pip install numpy scipy matplotlib
```

---

## Folder layout

```
analysis/
├── README.md              ← this file
├── compare.py             ← main entry point — run this
├── audio_io.py            ← WAV loading
├── fidelity_metrics.py    ← sliding-window distortion / SNR metrics
├── sibilance_metrics.py   ← HF sibilance-specific analysis
└── requirements.txt
```

---

## What `compare.py` does (pipeline)

1. **Load** both WAV files (`audio_io.py`). Float WAV → temporary PCM via ffmpeg if needed.
2. **Align** in time via cross-correlation on band-limited downsampled audio (handles different start trims / leader).
3. **8 kHz low-pass** (4th-order Butterworth) on both — models 1960s cinema playback bandwidth. Use `--no-lp` to skip.
4. **Gain match** — scales file B to match file A RMS on the overlap, trimming ±10 s from each end. Use `--no-gain-match` to skip.
5. **Compute metrics** — see tables below.
6. **Write report** — `restoration_AB_report.md`, plots, `comparison_summary.json`.

---

## Module reference

### `audio_io.py`

| Function | Purpose |
|---|---|
| `read_wav(path)` | Load mono float32, DC removed. Supports 8/16/24/32-bit PCM. |
| `ensure_readable_wav(path, tmp_dir)` | If stdlib `wave` cannot open the file (e.g. IEEE float), transcode to PCM via ffmpeg. |
| `EDGE_TRIM_SEC` | Seconds trimmed from start/end of timelines (default 10). |

### `fidelity_metrics.py`

Sliding-window analysis: **2 s windows, 1 s hop**. Reports median and p95 over the trimmed timeline.

| Metric | Key | Direction | What it measures |
|---|---|---|---|
| **Composite stress** | `composite_stress_trim` | ↓ better | Weighted blend of distortion proxies |
| **IMD / cross-mod** | `imd_db` | ↓ | Intermodulation at f₂±f₁, f₂±2f₁ from strong low + high partials. Classic variable-area nonlinearity. |
| **Inharmonic residual** | `inharm_db` | ↓ | Spectral peaks not locked to harmonics — scrape, heterodynes, dirt |
| **HF harshness** | `hf_harsh_db` | ↓ | 6–14 kHz vs 0.3–2.5 kHz energy. **Listen correlate: sibilance fizz** |
| **Spectral flatness** | `flatness` | ↓ | 0.4–8 kHz noise-like vs tonal character |
| **Noise floor / SNR** | `noise_floor_db` | ↑ | p90/p10 spectral ratio in 0.5–8 kHz |
| **Click index** | `impulse_idx` | ↓ | Large sample-to-sample jumps per 10k samples |
| **AM complexity** | `mod_complex_db` | ↓ | Broadband amplitude modulation 2–200 Hz (excludes 24 Hz frame band) |
| **Aliasing proxy** | `alias_db` | ↓ | Very-high vs mid-band energy |
| **Crest factor** | `crest_factor` | ↓ | Peak / RMS per window |

`analyze_file(x, sr)` returns time series + stats. `plot_timelines()` writes a 4-panel PNG.

### `sibilance_metrics.py`

**Why separate?** Broadband IMD can rate two transfers similarly while one retains ultrahigh splatter on sibilants that sounds distorted. Analysed on **raw aligned audio** (before 8 kHz LP) because the phenomenon is energy above 8 kHz during sibilant events.

| Metric | Key | Direction | What it measures |
|---|---|---|---|
| **HF IMD** | `hf_imd_p95` | ↓ | Cross-mod products landing in 4–8 kHz |
| **Sibilance buzz** | `sib_buzz_p95` | ↓ | Non-harmonic vs harmonic energy in sibilant windows |
| **Sib splatter 6.5–8 kHz** | `sib_6p5_8k_vs_core_db_median` | ↓ | Upper sibilance vs 4–8 kHz core, on sibilant frames |
| **Sib splatter 8–12 kHz** | `sib_8_12k_vs_core_db_median` | ↓ | **Key metric for "distorted sibilance"** — ultrahigh fizz |
| **Sib splatter 12–20 kHz** | `sib_12_20k_vs_core_db_median` | ↓ | Digitization noise / unfiltered HF on sibilants |

`detect_sibilant_frames()` finds 50 ms windows where 4–8 kHz dominates over mids.

### `compare.py` (extra metrics inline)

| Metric | Key | Direction | What it measures |
|---|---|---|---|
| **Rumble** | `rumble_db` | ↓ | 20–80 Hz vs 300–2 kHz — tracking / LF contamination |
| **Flutter / wow** | `flutter_proxy_db` | ↓ | FM sidebands ±2–24 Hz around dominant partial |
| **Envelope SNR** | `snr_modulation_db` | ↑ | Hilbert envelope vs residual |
| **Modulation depth** | `mod_depth` | context | Envelope std / mean — variable-area contrast |
| **Dynamic range** | `dynamic_range_db` | ↑ | p99 − p10 of envelope in dB |
| **% energy >8 kHz** | `pct_above_8k_hz` | ↓ | How much HF survives (workflow difference) |

---

## Outputs (in `--out` directory)

| File | Content |
|---|---|
| `restoration_AB_report.md` | Human-readable A/B table + ranking |
| `comparison_summary.json` | All numbers for scripting |
| `metrics_comparison.png` | Bar chart of key metrics |
| `spectrum_overlay.png` | PSD after 8 kHz LP |
| `fidelity_timeline_AB.png` | IMD, inharmonic, HF harshness, SNR over time |
| `*_fidelity_timeline.png` | Per-file 4-panel timelines |
| `sibilance/hf_psd_raw.png` | Raw HF spectrum 2–16 kHz |
| `sibilance/hf_spectrogram_30s.png` | HF spectrogram comparison |

---

## Interpreting results

### When transfers differ in gain and bandwidth

Typical situation (as with PEKARUV vs DupNeg Spinner):

- One transfer is **louder** and may already have **8 kHz low-pass** applied.
- The other retains trace energy **above 8 kHz**.

The pipeline compensates by LP + gain matching before most metrics. Sibilance splatter is measured on **raw** audio because that's where the audible difference lives.

### What to trust

| Question | Best metrics | Always verify by ear |
|---|---|---|
| Cross-modulation on loud passages | IMD p95, composite stress | ✓ |
| Hiss / noise floor | noise_floor_db, spectral flatness | ✓ quiet scenes |
| Sibilance quality | **hf_harsh_db**, **sib_8_12k_vs_core_db** | ✓ dialogue |
| Rumble / tracking | rumble_db | ✓ |
| Transport stability | flutter_proxy_db | ✓ steady tones |

### Caveats

- Metrics are **proxies**, not perceptual scores. Low waveform correlation between workflows is normal (different contrast curves).
- **THD proxy** on heavily filtered material can be meaningless (negative dB).
- **Suggested starting point** is a weighted heuristic — override with listening.
- Does **not** measure 24 Hz frame imprint (see `itst_tests/wav_artifact_analysis.py` for that).

---

## Example: LR01 findings (Jun 2026)

| Aspect | Sondor PEKARUV | ITST DupNeg |
|---|---|---|
| Rumble | Better (−20 dB) | Worse |
| Broadband IMD | Slightly better | — |
| Spectral SNR | — | Better |
| **Sibilance / HF** | More 8–20 kHz splatter | Cleaner (pre-LP'd) |
| **Listening** | Sibilants can sound grittier | Sibilants cleaner |

For restoration targeting 8 kHz cinema playback, **ITST may be preferable if sibilance matters**; **Sondor if LF cleanliness matters more**.

---

## CLI options

```
python3 analysis/compare.py FILE_A FILE_B --out OUTPUT_DIR
  --label-a NAME    Display name for file A (default: filename stem)
  --label-b NAME    Display name for file B
  --no-lp           Skip 8 kHz low-pass
  --no-gain-match   Skip RMS gain matching
```

---

## Related code elsewhere in the repo

| Path | Purpose |
|---|---|
| `itst_tests/wav_fidelity_analysis.py` | Original batch fidelity tool (multiple source folders) |
| `itst_tests/wav_artifact_analysis.py` | 24 Hz frame-imprint analysis |
| `itst_tests/restoration_digitization_compare.py` | Earlier single-file prototype of this comparison |

The `analysis/` folder is the **self-contained, documented** version intended for reuse.
