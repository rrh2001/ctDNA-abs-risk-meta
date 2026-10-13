# Reproducibility Runbook

This document describes how to reproduce every analysis table and figure in
**"Time-resolved absolute recurrence risk and longitudinal ctDNA trajectories
across the colorectal cancer course: a reconstructed-IPD meta-analysis with an
open risk calculator."**

All paths are relative to the repository root. Open `ctdna-abs-risk-meta.Rproj`
in RStudio (or otherwise set the working directory to the repository root)
before running anything.

## 1. Software and dependencies

* **R 4.6.1** (R 4.3+ should work), with packages:
  `survival`, `dplyr`, `tidyr`, `purrr`, `readr`, `ggplot2`, `metafor`,
  `flexsurv`, `future`, `furrr`, `patchwork`, `shiny`, `IPDfromKM`.
* **Python 3.13** (3.10+ should work), with:
  `pandas`, `numpy`, `scipy`, `matplotlib`, `pillow`.
* The figures are rendered with the base/Arial or DejaVu Sans fonts; the
  Python figure scripts auto-detect Windows/macOS/Linux font locations and fall
  back to the default font.

Install R packages if needed:

```r
install.packages(c("survival","dplyr","tidyr","purrr","readr","ggplot2",
                   "metafor","flexsurv","future","furrr","patchwork","shiny"))
install.packages("IPDfromKM")   # optional; used only as an independent
                                 # cross-check of the in-house Guyot algorithm
```

## 2. One-command run

```r
source("run_all.R")
```

This runs the R modules in order (M1 -> M2 -> M3 -> M4). Then run the Python
steps listed in Section 4. The Shiny app under `app/` is not rebuilt by
`run_all.R`; run `python python/setup_app.py` to refresh its data bundle.

## 3. R pipeline (order, inputs, outputs, runtime)

| Step | Script | Reads | Writes (`outputs/`) | Approx. time |
|---|---|---|---|---|
| M1a | `R/M1a_reconstruct_IPD.R` | `data/meta/{curves,nar,km_meta}.csv`, `data/digitized/*.csv` | `M1/M1_ipd_reconstructed.csv` (11,481 rows), `M1_event_counts_reconstructed.csv`, `M1_reconstruction_diagnostics.csv` | 2–4 min |
| M1b | `R/M1b_two_stage_absrisk.R` | M1a outputs, `data/meta/km_meta.csv` | `M1/M1_absrisk_perstudy.csv`, `M1_absrisk_pooled.csv`, `M1_contrast_pooled.csv`, `M1_core_results_DFS.csv` | 1–2 min |
| M1c | `R/M1c_onestage_RP.R` | M1a/M1b outputs | `M1/M1c_continuous_curves.csv`, `M1c_onestage_vs_twostage.csv` | 2–4 min |
| M1d | `R/M1d_sensitivity_external.R` | M1b outputs, `data/meta/km_meta.csv` | `M1/S3_reconstructed_vs_published_HR.csv`, `S4_leave_one_study_out.csv`, `S4_egger_small_study.csv`, `S5_external_calibration.csv`* | 1–2 min |
| M1e | `R/M1e_meta_regression.R` | M1b outputs, `data/meta/study_moderators.csv` | `M1/S4_meta_regression.csv`, `S4_subgroup_RD.csv` | <1 min |
| M1f | `R/M1f_figures.R` | M1/M2 outputs | `figures/Figure_2_absolute_risk_curves`, `figures/Figure_3_RD_forest_plot` (600 dpi PNG + LZW TIFF) | <1 min |
| M2a | `R/M2a_reconstruct_traj.R` | `data/meta/*` | `M2/M2_ipd_trajectory_reconstructed.csv`, `M2_trajectory_absolute_risk.csv`, `M2_event_counts.csv`, `M2_traj_HR_narrative_EN.csv` | 2–3 min |
| M2b | `R/M2b_figures.R` | M2 outputs | `figures/Figure_4_trajectories` (600 dpi PNG + LZW TIFF; fixed jitter seed) | <1 min |
| M3 | `R/M3_timing.R`, `R/M3_figures.R` | `M1/M1c_continuous_curves.csv` | `M3/M3_risk_timing.csv`; `figures/Figure_5_surveillance_framework` (600 dpi PNG + TIFF). The lead-time table is written by the Python step below | <1 min |
| M4a | `R/M4_absrisk_reconstruction_thresholds.R` | self-contained simulation | `M4/M4_threshold_summary.csv`, `M4_threshold_by_events.csv`, `M4_threshold_gate_dual_evidence.csv` | **20–40 min full**; ~1 min with `TEST_MODE <- TRUE` |
| M4b | `R/M4b_empirical_IPD_validation.R` | third-party appendix IPD* | `M4/M4b_empirical_threshold.csv`* | ~1 min smoke |

\* Steps M1d-P4 and M4b consume cleaned, patient-level appendix IPD from Li
(2022) and Ogaard (2022), which are **not redistributed** (see README, Data
availability). The scripts detect their presence; if the files are absent they
print a skip message and **retain the released tables** (`S5_external_calibration.csv`,
`M4b_empirical_threshold.csv`) instead of failing.

### M4a test mode

For a quick check open `R/M4_absrisk_reconstruction_thresholds.R`, set
`TEST_MODE <- TRUE` (line ~25), and source. The released M4 tables are from the
full 144-grid run with `TEST_MODE <- FALSE` and `B = 1000`; the master seed is
`MASTER_SEED <- 20261013`.

## 4. Python steps

Run after the R pipeline (M6 uses the M1 outputs and the released S5 table).

```bash
python python/M3_leadtime_lit.py               # M3 published lead-time table (self-contained)
python python/M6_calibration_competing.py      # competing-risk + internal/external calibration
python python/make_figure1_prisma.py           # Figure 1 PRISMA flow (600 dpi)
python python/make_figure6_calibration.py      # Figure 6 calibration (600 dpi)
python python/setup_app.py                     # refresh app/data from outputs/M1
```

M6 writes `outputs/M6/`: `M6_competing_risk_sensitivity.csv`,
`M6_internal_calibration_pairs.csv`, `M6_internal_calibration_metrics.csv`,
`M6_calibration_deciles.csv`, `S7_external_calibration_summary.csv`,
`M6_summary.txt`. The numpy seed is `20261014`.

## 5. Upstream ETL (one-time, not needed for reproduction)

`python/upstream_etl/` holds four one-time preparation scripts that built
the English analysis bundles in `data/meta/` from the maintainers' original
digitization workbooks and text extractions. They are archived for provenance
and require the original Chinese-language source workbooks (see
`python/upstream_etl/README.md`). They are **not** part of the reproducible
chain and should not be run against the released repository. The self-contained
M3 lead-time compilation (`M3_leadtime_lit.py`) lives in `python/` and is part
of the reproducible chain.

## 6. Verification performed on the released repository

A clean rerun on the released repository (R 4.6.1, the listed packages, fixed
seeds) was compared cell-by-cell against the released `outputs/`:

* M1 (IPD reconstruction, event counts, two-stage pooled absolute risk,
  contrasts, S3 reconstructed-vs-published HR, S4 leave-one-study-out/Egger/
  meta-regression/subgroups), M2 (trajectory IPD and absolute risk) and M3
  (lead-time tables) reproduced to a maximum numeric difference of **0.000000**.
* M6 (competing-risk sensitivity, internal LOO calibration, external summary)
  reproduced to **0.000000**.
* The one-stage Royston-Parmar spline tables (M1c) showed a maximum point-level
  difference of ~0.015 attributable to platform-dependent MLE/BLAS numerical
  optimization; no reported pooled estimate, interval, or conclusion changes.
* The released outputs contain exactly 35 files; the rerun reproduced the same
  file set. The Shiny data bundle (`app/data/continuous_curves.csv`, 1,424 rows)
  matches `M1c_continuous_curves.csv`.
* Figures 1, 2, 3, 5 and 6 reproduce pixel-for-pixel from their scripts. Figure
  4B uses a fixed jitter seed (`set.seed(20261013)` in `M2b_figures.R`), so two
  consecutive runs are also pixel-identical; the seed only fixes the
  overlap-avoiding horizontal placement and changes no data. All six figures are
  written to `figures/` as 600 dpi PNG and LZW TIFF (12 files).

## 7. Troubleshooting

* **A step cannot find `data/meta/...` or writes to a stray `M1_out/`,
  `M2_out/`, `M3_out/` or `M4_out/` folder:** the working directory is not the
  repository root. Open the `.Rproj` or `setwd()` to the root.
* **"P4/M4b skipped" messages:** expected; the third-party appendix IPD are not
  shipped. The released S5/M4b tables are kept.
* **M4a is slow:** run with `TEST_MODE <- TRUE` first; the full grid is
  intentionally expensive and parallelizes with `future::multisession`.
* **Font warnings in figures:** cosmetic; the scripts fall back to the default
  font and the layout is unchanged.
