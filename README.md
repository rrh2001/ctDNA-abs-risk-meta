# Time-resolved absolute recurrence risk of ctDNA across the colorectal cancer course

**A reconstructed individual-participant-data meta-analysis with an open, web-based risk calculator**

Ren Ruhua, Gao Ruike, Li Jie*

Department of Oncology, Guang'anmen Hospital, China Academy of Chinese Medical Sciences, Beijing, China

\*Corresponding author: Li Jie — qfm2020jieli@yeah.net

Protocol: PROSPERO **CRD420261516792** (registration and amendment)
**Online calculator:** <https://renruhua.shinyapps.io/ctDNA_risk_calculator/>

> ### Before publication / peer review — placeholders to replace
>
> 1. **Zenodo archive DOI.** Create a GitHub release, enable the repository in
>    Zenodo, and paste the DOI Zenodo mints into this README and `CITATION.cff`
>    (search for `ZENODO_DOI`).
> 2. **Confirm the public GitHub URL** below matches the repository you create:
>    `https://github.com/renruhua/ctdna-abs-risk-meta`.
> 3. The manuscript's Data Availability statement must use the same URLs/DOI.

---

## 1. What this is

Published studies report the prognostic performance of circulating tumour DNA
(ctDNA) almost exclusively as hazard ratios or two-group significance tests.
Clinicians, however, need **absolute, time-resolved recurrence risks** ("what is
the probability of recurrence at 12 / 24 / 36 months for a ctDNA-positive versus
ctDNA-negative patient *at this point* on the treatment pathway?") and these are
not available as a pooled, course-wide summary.

We obtain patient-level time-to-event data from published Kaplan–Meier curves
using the validated Guyot digitisation-reconstruction algorithm, pool
curative-intent settings with a pre-specified reliability gate (Module M4), and
produce:

- **absolute (cumulative incidence) risk** for ctDNA-positive and ctDNA-negative
  patients at every landmark (post-surgery, post-(neo)adjuvant, post-liver
  resection, during surveillance) and the between-group **risk difference**;
- **longitudinal trajectory risks** (sustained-negative, clearance, newly
  emergent, persistent);
- an **open Shiny calculator** that returns absolute risk curves from the
  one-stage Royston–Parmar spline models;
- a **dual-evidence reliability gate** (parametric simulation **and**
  resampling of real published appendix IPD) defining the minimum digitised
  points, number-at-risk table density and event count required for valid
  reconstruction;
- quantitative **competing-risk** sensitivity, internal leave-one-study-out
  calibration and external calibration against two published appendix IPD sets.

This repository contains everything needed to reproduce the analysis from the
digitised curves and the released data bundles.

## 2. Repository layout

```
ctdna-abs-risk-meta/
├── data-raw/
│   └── study_sources.csv          # 24 included studies: citations, DOI/PII, scope
├── data/
│   ├── digitized/                 # 164 digitised KM curves (time, S(t)) per KMID
│   └── meta/                      # English, sealed analysis bundles (M1/M2 inputs)
├── run_all.R / *.Rproj            # One-command R pipeline (open .Rproj, source run_all.R)
├── R/                             # Main, fully reproducible pipeline (12 scripts)
├── python/                        # M3 lead-time, M6 calibration, figures, app setup
│   └── upstream_etl/              # ONE-TIME ETL reading the original Chinese master
│                                  # workbook / third-party IPD (not re-run; see §5)
├── app/                           # Shiny calculator (app.R, tests, bundled data)
├── outputs/                       # All released result tables (M1–M6, S1–S7)
│   ├── M1  M2  M3  M4  M6
│   ├── S1_QUIPS  (risk-of-bias)
│   └── S2_CINeMA (confidence rating)
├── figures/                       # Figures 1–6, 600 dpi PNG + LZW TIFF
├── supplementary/                 # Released supplementary workbook (S6/S7)
└── docs/                          # Reproduction runbook and GitHub/Zenodo upload guide
```

The complete column-level dictionary is in **[CODEBOOK.md](CODEBOOK.md)**.

## 3. Reproducing the analysis

### 3.1 Software

- **R ≥ 4.3** with packages:
  `survival`, `dplyr`, `tidyr`, `purrr`, `readr`, `ggplot2`, `metafor`,
  `flexsurv`, `future`, `furrr`, `patchwork`, `shiny` (for the app), and
  optionally `IPDfromKM` (used only for an additional cross-check; the validated
  Guyot main algorithm runs without it).
- **Python ≥ 3.10** with `pandas`, `numpy`, `scipy`, `matplotlib`, `pillow`
  (for Module M6 and Figures 1 and 6).

### 3.2 Run order

Main pipeline (set the working directory to the repository root). The easiest
way is to open `ctdna-abs-risk-meta.Rproj` and run `source("run_all.R")`, which
executes the R scripts in order, then run the Python steps in rows 12–14. All
scripts read the released English bundles in `data/` and write to `outputs/`:

| Order | Script | Produces |
|---|---|---|
| 1 | `R/M1a_reconstruct_IPD.R` | Guyot reconstructed IPD (`M1_ipd_reconstructed.csv`) |
| 2 | `R/M1b_two_stage_absrisk.R` | two-stage absolute-risk pooling, RD/RR |
| 3 | `R/M1c_onestage_RP.R` | one-stage Royston–Parmar continuous curves |
| 4 | `R/M1d_sensitivity_external.R` | S3 HR validation, LOO, Egger, external calibration |
| 5 | `R/M1e_meta_regression.R` | meta-regression and subgroup analyses |
| 6 | `R/M1f_figures.R` | Figures 2 and 3 |
| 7 | `R/M2a_reconstruct_traj.R` | longitudinal trajectory IPD and absolute risks |
| 8 | `R/M2b_figures.R` | Figure 4 |
| 9 | `R/M3_timing.R`, `R/M3_figures.R` | risk timing table, Figure 5 |
| 10 | `R/M4_absrisk_reconstruction_thresholds.R` | parametric reliability gate (M4a) |
| 11 | `R/M4b_empirical_IPD_validation.R` | real-IPD gate (M4b; needs appendix IPD, see §5) |
| 12 | `python/M3_leadtime_lit.py`, `python/M6_calibration_competing.py` | M3 lead-time table; competing-risk, internal/external calibration |
| 13 | `python/make_figure1_prisma.py`, `python/make_figure6_calibration.py` | Figures 1 and 6 |
| 14 | `python/setup_app.py` | stage `app/data/continuous_curves.csv` |

The detailed runbook (seeds, run times, expected outputs) is in
**[docs/runbook.md](docs/runbook.md)**. A clean rerun of the released repository
reproduced the M1–M3 and M6 tables to a maximum numeric difference of
**0.000000** (fixed seeds); only the one-stage spline curves can move by ≤0.015
across BLAS/platforms, with no change to any reported estimate or conclusion
(runbook §6).

Step-by-step instructions for publishing this snapshot on GitHub and minting
the Zenodo DOI are in **[docs/UPLOAD_to_GitHub_Zenodo.md](docs/UPLOAD_to_GitHub_Zenodo.md)**.

### 3.3 Running the calculator locally

```r
shiny::runApp("app")
```

or simply open the hosted version at
<https://renruhua.shinyapps.io/ctDNA_risk_calculator/>.

## 4. Data availability

**Provided here.** The 164 digitised curves (`data/digitized/`) are coordinates
extracted from the Kaplan–Meier figures of the 24 included public studies; the
sealed analysis bundles (`data/meta/`) and every released result table
(`outputs/`) are included. `data-raw/study_sources.csv` lists every included
study with its DOI / PII / trial registration so the source article can be
located and the extraction audited.

**Not redistributed (third-party patient-level data).** Two studies (Li 2022;
Ogaard 2022) published patient-level appendix data that we use **only** for
validation (M4b empirical reconstruction and the external-calibration checks).
These are third-party data used for academic validation; they are not our data
and are deliberately excluded. They can be reconstructed from the published
appendices or obtained from the corresponding author of the original reports,
subject to the publishers' terms. The cleaning routine
(`python/upstream_etl/clean_real_ipd.py`) is released for transparency but its
input files are not committed.

**Large simulation outputs.** The raw per-iteration M4 grids
(`M4_all_iterations.csv`, `M4b_empirical_all.csv`, ~150 MB total) are excluded
from Git. They are reproducible from the fixed master seed
(`MASTER_SEED = 20261013`) and are available via the Zenodo archive (see
`ZENODO_DOI`) or on request. The aggregated, reported tables are in
`outputs/M4/`.

## 5. Upstream ETL scripts (transparency, not part of the sealed chain)

`python/upstream_etl/` contains one-time extraction scripts that read the
authors' original master extraction workbook (Chinese worksheet/column names at
the time of extraction) and the third-party appendix IPD. They are provided for
provenance; **the sealed, fully reproducible analysis starts from `data/` and
`R/`** and does not require them. Chinese string literals in those scripts are
deliberate matchers against the original source files, not user-facing text.

## 6. Methods standards and reporting

- PRISMA 2020 flow diagram (Figure 1) and checklist (Supplementary Information).
- Risk of bias: QUIPS (S1). Confidence in absolute-risk estimates: CINeMA (S2).
- Reconstructed-vs-published hazard-ratio validation: S3 (52 matched; 50 within
  the published 95% CI).
- Heterogeneity is reported with I², τ² and prediction intervals; influence
  diagnostics, Egger tests, leave-one-study-out and meta-regression are in S4.
- Competing-risk and calibration analyses (S6/S7) follow the manuscript Methods.

## 7. Citation

See [CITATION.cff](CITATION.cff). Please cite the manuscript once assigned and
the Zenodo archive for the code/data snapshot. DOI: `ZENODO_DOI` (to be added
after the first release).

## 8. License

Code is released under the **MIT License** ([LICENSE](LICENSE)). The released
tables and figures are made available for academic use under the same terms; the
underlying studies retain their original copyright. Third-party patient-level
data are not included.

## 9. Disclaimer

The calculator and all outputs are for **research and educational use only**,
do not constitute individual medical advice, and must not replace clinical
judgement. The surveillance timing material in the manuscript is presented as a
**hypothesis-generating framework**, not a clinical guideline.
