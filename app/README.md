# ctDNA Time-Resolved Absolute Recurrence Risk Calculator

A Shiny web app that turns the meta-analysis' continuous absolute-risk curves (one-stage
Royston-Parmar fits on **reconstructed patient-level data** from published Kaplan–Meier curves,
24 studies) into an interactive, time-resolved recurrence-risk tool for colorectal cancer.

## Contents
```
ctDNA_risk_calculator/
├── app.R                      # self-contained Shiny UI + server
├── data/
│   └── continuous_curves.csv  # engine = M1c continuous CIR curves (7 DFS cells)
└── README.md
```

## Run locally (RStudio)
1. Install R (>= 4.2) and RStudio.
2. Install packages once (in the R console):
```r
install.packages(c("shiny", "readr", "dplyr", "ggplot2"))
```
3. Open `app.R` in RStudio and click **Run App** (top-right of the editor), or in the console:
```r
setwd("the/folder/that/contains/app.R")
shiny::runApp()
```

## Deploy to shinyapps.io (free, public URL)
```r
install.packages("rsconnect")
library(rsconnect)
setAccountInfo(name = "<your-account>", token = "<token>", secret = "<secret>")  # from your shinyapps.io dashboard
setwd("the/folder/that/contains/app.R")
deployApp(appName = "ctdna-crc-risk", appTitle = "ctDNA CRC Recurrence Risk Calculator")
```
After the first deploy, subsequent updates are simply `deployApp()`. The app is dependency-light
(no external API, no database), so the free tier is sufficient.

## What it shows
- **Inputs**: disease setting (non-metastatic CRC / locally-advanced rectal cancer / resectable
  colorectal liver metastases), assessment window (pre-operative / post-op MRD 4–12 wk / after
  therapy), ctDNA status, and a month slider.
- **Outputs**: selected-horizon cumulative incidence of recurrence (CIR) with 95% CI,
  positive-vs-negative contrast (risk difference and risk ratio), a continuous CIR plot with
  95% bands, and 12/24/36-month landmark table.

## Scope, reliability and data statement
- Only cells with **≥3 independent studies, grade-A curves, complete risk tables** and passing the
  M4 dual-evidence reliability gate (MAE ≤ 0.05, |relative bias| ≤ 0.10, coverage between 0.90 and 0.98) are exposed.
  Advanced-disease and anti-EGFR rechallenge settings (<3 studies) are deliberately excluded.
- Confidence bands are individual-level (narrower than random-effects CIs); the manuscript's
  two-stage random-effects CIs are the primary between-study uncertainty.
- Reconstruction used published Kaplan–Meier curves and risk tables. The published appendix
  individual-level datasets of **Li 2022 and Øgaard 2022 were used for validation only and are not
  redistributed**; obtain them from the original publications.
- **For research and educational use only; not a medical device and not for individual clinical
  decisions.** Protocol: PROSPERO CRD420261516792.
