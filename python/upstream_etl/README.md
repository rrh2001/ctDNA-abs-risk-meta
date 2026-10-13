# Upstream ETL (one-time, archived for provenance)

The four scripts in this folder are **one-time data-preparation scripts** that
built the English analysis bundles shipped under `data/meta/`. They are
preserved for provenance and audit but are **not part of the sealed
reproduction chain**, which starts from `data/` and the `R/` modules. You do
**not** need to run them to reproduce the paper; running them requires the
maintainers' original, non-redistributed source workbooks.

| Script | Purpose | Inputs (not shipped) |
|---|---|---|
| `M1_prep.py` | Builds the M1 curve/narrative/KM-meta bundles (`curves.csv`, `nar.csv`, `km_meta.csv`) from the master extraction workbook and digitized-curve folder. | original Chinese-language master workbook, digitized curve folder |
| `M2prep.py` | Builds the M2 longitudinal-trajectory extraction bundles. | original Chinese-language master workbook |
| `M2_hr_narrative.py` | Extracts the trajectory hazard-ratio narrative table (`M2_traj_HR_narrative_EN.csv`). | original Chinese-language master workbook |
| `clean_real_ipd.py` | Cleans the published Li 2022 / Ogaard 2022 appendix patient-level data into `ipd_clean_*.csv`. | third-party appendix `.tsv` files (set `REAL_IPD_DIR`) |

## Why some scripts contain Chinese string literals

`M1_prep.py`, `M2prep.py` and `M2_hr_narrative.py` match against the original
worksheets, column headers and group labels of the maintainers' master
workbooks, which were written in Chinese. The Chinese string literals in these
scripts are therefore **deliberate matchers against those source files**, not
user-facing output. Every file committed to the repository (data, outputs,
figures, app, user-facing documentation) is in English.

## Input path placeholders

These scripts read local paths from environment variables (`CTDNA_WORK_DIR`,
`CTDNA_MASTER_XLSX`, `CTDNA_DIGITIZED`, `CTDNA_MATRIX`, `CTDNA_REAL_IPD`) that
default to a `<LOCAL_PATH>/...` placeholder. To re-run them, point those
variables at your own copy of the source files. The original Chinese-language
master workbook and the digitized curves are available from the corresponding
author on reasonable request, subject to the publisher's terms; the third-party
appendix IPD must be obtained directly from the original publications and must
not be redistributed.

## Self-contained scripts are elsewhere

The self-contained, reproducible Python scripts (which need no external
workbooks) live in the parent `python/` directory, not here:
`M3_leadtime_lit.py`, `M6_calibration_competing.py`,
`make_figure1_prisma.py`, `make_figure6_calibration.py`, `setup_app.py`.
