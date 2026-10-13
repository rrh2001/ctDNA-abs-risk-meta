# =============================================================================
# Master run script
#
# Open ctdna-abs-risk-meta.Rproj in RStudio (this sets the working directory to
# the repository root), then run:  source("run_all.R")
#
# Reproducibility
#   * Deterministic steps reproduce the released tables in outputs/ exactly
#     (fixed seeds; the released M1/M2/M3/S3/S4 tables were verified to match a
#     clean rerun to the 6th decimal).
#   * The one-stage Royston-Parmar spline (M1c) is a numerical MLE fit; across
#     different BLAS/operating systems a continuous-curve point can move by at
#     most ~0.015, with no change to any reported estimate or conclusion.
#
# Long / optional steps
#   * R/M4_absrisk_reconstruction_thresholds.R runs the full 144-grid simulation
#     with TEST_MODE <- FALSE (20-40 min). Set TEST_MODE <- TRUE for a ~1 minute
#     smoke run; the released outputs/M4 tables are from the full run.
#   * M1d P4 and M4b use third-party published appendix IPD that are NOT
#     redistributed in this repository. When those files are absent the steps
#     print a skip message and retain the released S5 / M4b threshold tables.
#
# Figures: R/M1f_figures.R, R/M2b_figures.R and R/M3_figures.R write the
# submission-ready figures directly to figures/ (Figure_2-5; 600 dpi PNG and
# LZW TIFF). Figure_1 and Figure_6 are produced by the Python steps below.
#
# Python steps (run after the R pipeline; Python 3 with pandas/numpy/scipy/
# matplotlib/pillow):
#   python python/M3_leadtime_lit.py
#   python python/M6_calibration_competing.py
#   python python/make_figure1_prisma.py
#   python python/make_figure6_calibration.py
#   python python/setup_app.py
# See docs/runbook.md for the full order and expected outputs.
# =============================================================================

scripts <- c(
  "R/M1a_reconstruct_IPD.R",
  "R/M1b_two_stage_absrisk.R",
  "R/M1c_onestage_RP.R",
  "R/M1d_sensitivity_external.R",
  "R/M1e_meta_regression.R",
  "R/M1f_figures.R",
  "R/M2a_reconstruct_traj.R",
  "R/M2b_figures.R",
  "R/M3_timing.R",
  "R/M3_figures.R",
  "R/M4_absrisk_reconstruction_thresholds.R",
  "R/M4b_empirical_IPD_validation.R"
)

for (s in scripts) {
  message("\n==================  source ", s, "  ==================")
  tryCatch(source(s), error = function(e) message("!! ERROR in ", s, ": ", conditionMessage(e)))
}
message("\nR pipeline finished. Next run the Python steps listed in docs/runbook.md.")
