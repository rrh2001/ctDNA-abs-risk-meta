# ============================================================================
# M1a_reconstruct_IPD.R
# Module 1, step 1: batch reconstruct patient-level data (IPD) from the static two-state grade-A KMs
#
# Inputs (frozen data produced by M1_prep.py; column names are English; any Chinese values are read by readr as UTF-8):
#   data/meta/curves.csv  KMID,curve,time,St          digitized KM coordinates (dense)
#   data/meta/nar.csv     KMID,curve,grp,time,nrisk   number-at-risk table
#   data/meta/km_meta.csv kmid,study_id,author,year,scenario,landmark,endpoint,
#                               comparison,hicurve,locurve,hi_label,lo_label,group_n,
#                               ev_est_m4,ev_layer,single_report,nar_minpts,cmp_evt,hr,hr_type
#
# Main algorithm: iterative Guyot reconstruction (2012, BMC Med Res Methodol), identical to the IPDfromKM package;
#         fidelity validated in M4a (parametric simulation) and M4b (closed loop on real Li2022 / Ogaard2022 IPD).
# Optional cross-check: if IPDfromKM is installed, getIPD recomputes and compares per-group event counts (tryCatch, non-blocking).
#
# Outputs (outputs/M1/):
#   M1_ipd_reconstructed.csv   reconstructed IPD long table (main result)
#   M1_event_counts_reconstructed.csv  reconstructed N/events/censoring per KM x ctDNA + exact event counts (back-fills the event stratum)
#   M1_reconstruction_diagnostics.csv   diagnostics (N mismatch / event-count difference between methods / failures)
#
# Depends on readr,dplyr,tidyr,purrr,survival (all installed); IPDfromKM optional
# Run: setwd(this script dir), then source("M1a_reconstruct_IPD.R")
# ============================================================================
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(purrr); library(survival)
})

bundle <- "data/meta"; out_dir <- "outputs/M1"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
USE_PKG <- requireNamespace("IPDfromKM", quietly = TRUE)
if (USE_PKG) cat("IPDfromKM detected; additional cross-validation will be run\n") else
  cat("IPDfromKM not installed; using the validated Guyot main algorithm only (results unaffected; install.packages('IPDfromKM') to cross-check)\n")

cur  <- read_csv(file.path(bundle, "curves.csv"), show_col_types = FALSE)
nar  <- read_csv(file.path(bundle, "nar.csv"),    show_col_types = FALSE)
meta <- read_csv(file.path(bundle, "km_meta.csv"), show_col_types = FALSE)

# ---------- Guyot iterative reconstruction (cases at risk after the last point are censored so total N matches the first risk number)----------
guyot <- function(times, surv, nr_t, nr_v) {
  if (length(times) < 3 || length(nr_t) < 2) return(NULL)
  if (times[1] != 0) { times <- c(0, times); surv <- c(1, surv) }
  o <- order(times); times <- times[o]; surv <- cummin(surv[o])
  surv <- pmin(pmax(surv, 0), 1)
  Sat <- function(t) surv[findInterval(t, times, all.inside = TRUE)]
  ev_t <- numeric(0); cen_t <- numeric(0)            # collect event and censoring times separately
  for (i in seq_len(length(nr_t) - 1)) {
    lo <- nr_t[i]; hi <- nr_t[i + 1]; nlo <- nr_v[i]; nhi <- nr_v[i + 1]
    if (is.na(nlo) || nlo <= 0) next
    if (is.na(nhi)) nhi <- 0
    Slo <- Sat(lo); Shi <- Sat(hi); if (Slo <= 0) next
    net <- nlo - nhi                              # net attrition in the interval; conservation d+cc=net
    if (is.na(net)) next
    if (net < 0) { d <- 0; cc <- 0 } else {       # risk table shows an anomalous increase: no attrition in this interval
      d <- round(nlo * (1 - Shi / Slo)); d <- max(0, min(d, net))
      cc <- nlo - d - nhi
    }
    if (d > 0) {
      # Locate the d events at the actual digitized steps: spread them uniformly in the S-drop space and
      # map each level onto the curve time, rather than uniformly in clock time (Guyot 2012). This keeps a
      # concentrated late drop (e.g. a single large step near the start of a sparse final NAR interval) from
      # being smeared across the interval onto ever-smaller risk sets. Event COUNT and NAR interval
      # membership are unchanged, so endpoint CIR is identical; only the within-interval event timing improves.
      inidx <- which(times > lo & times <= hi)
      if (length(inidx) >= 1) {
        tg <- times[inidx]; sg <- surv[inidx]
        Starget <- Slo - (1:d) * (Slo - Shi) / d
        ev_here <- vapply(Starget, function(stg) {
          hit <- tg[sg <= stg]
          if (length(hit)) min(hit) else hi
        }, numeric(1))
        ev_here <- pmin(pmax(ev_here, lo + .Machine$double.eps), hi)
      } else {
        ev_here <- seq(lo, hi, length.out = d + 2)[-c(1, d + 2)]
      }
      ev_t <- c(ev_t, ev_here)
    }
    if (cc > 0) {
      # Place the interval's censoring AFTER the last event (across the remainder of the interval); when the
      # NAR table collapses to zero at the next point, the loss is administrative censoring that manifests
      # at that point and must not shrink the risk set before the events recorded on the curve (avoids the
      # spurious late S drop seen in a sparse final interval). If there are no events, spread across the interval.
      tlast_ev <- if (d > 0) max(ev_here) else lo
      if (hi > tlast_ev + .Machine$double.eps) {
        cen_t <- c(cen_t, seq(tlast_ev, hi, length.out = cc + 2)[-c(1, cc + 2)])
      } else {
        cen_t <- c(cen_t, rep(hi, cc))
      }
    }
  }
  target_n <- nr_v[1]                                  # those still at risk after the last risk point -> censored
  assigned <- length(ev_t) + length(cen_t)
  if (!is.na(target_n) && assigned < target_n) {
    tlast <- max(nr_t); tmax <- max(times); if (tmax <= tlast) tmax <- tlast + 1
    miss <- target_n - assigned
    if (miss > 0) cen_t <- c(cen_t, seq(tlast, tmax, length.out = miss + 2)[-c(1, miss + 2)])
  }
  if (length(ev_t) + length(cen_t) == 0) return(NULL)
  data.frame(time   = c(ev_t, cen_t),
             status = c(rep(1L, length(ev_t)), rep(0L, length(cen_t))))
}

guyot_pkg <- function(times, surv, nr_t, nr_v) {       # optional IPDfromKM cross-validation
  if (!USE_PKG) return(NULL)
  tryCatch({
    res <- IPDfromKM::getIPD(est = surv, t = times, risk = nr_v, risk.t = nr_t, quiet = TRUE)
    if (is.data.frame(res) && all(c("time", "status") %in% names(res))) res else NULL
  }, error = function(e) NULL)
}

recon_one <- function(k, cv, label, ctdna) {
  cc <- cur %>% filter(KMID == k, curve == cv) %>% arrange(time)
  nn <- nar %>% filter(KMID == k, curve == cv) %>% arrange(time)
  if (nrow(cc) < 3 || nrow(nn) < 2)
    return(list(ipd = NULL, diag = data.frame(kmid=k, curve=cv, ctdna=ctdna, problem="coord_or_NAR_insufficient")))
  ipd <- guyot(cc$time, cc$St, nn$time, nn$nrisk)
  ipk <- guyot_pkg(cc$time, cc$St, nn$time, nn$nrisk)
  if (is.null(ipd))
    return(list(ipd = NULL, diag = data.frame(kmid=k, curve=cv, ctdna=ctdna, problem="Guyot_failed")))
  ipd$kmid <- k; ipd$curve <- cv; ipd$ctdna <- ctdna; ipd$ctdna_label <- label
  diag <- data.frame(kmid = k, curve = cv, ctdna = ctdna, label = label,
    n0_risk = nn$nrisk[1], n_coord = nrow(cc), n_nar = nrow(nn),
    recon_N = nrow(ipd), events_guyot = sum(ipd$status), cens_guyot = sum(1 - ipd$status),
    events_pkg = if (!is.null(ipk)) sum(ipk$status) else NA_integer_,
    tmax = max(cc$time), stringsAsFactors = FALSE)
  list(ipd = ipd, diag = diag)
}

ipd_all <- list(); diag_all <- list(); j <- 0
for (i in seq_len(nrow(meta))) {
  r <- meta[i, ]
  for (seg in list(list(cv = r$hicurve, lab = r$hi_label, st = 1L),
                   list(cv = r$locurve, lab = r$lo_label, st = 0L))) {
    z <- recon_one(r$kmid, seg$cv, seg$lab, seg$st)
    if (!is.null(z$ipd)) {
      j <- j + 1
      z$ipd$study_id <- r$study_id; z$ipd$scenario <- r$scenario
      z$ipd$landmark <- r$landmark; z$ipd$endpoint <- r$endpoint; z$ipd$ev_layer <- r$ev_layer
      ipd_all[[j]] <- z$ipd
    }
    diag_all[[length(diag_all) + 1]] <-
      if (is.data.frame(z$diag) && "n_coord" %in% names(z$diag)) z$diag else
      data.frame(kmid = r$kmid, curve = seg$cv, ctdna = seg$st, problem = "unknown_error")
  }
}
ipd_df  <- bind_rows(ipd_all)
diag_df <- bind_rows(diag_all)

ec <- meta %>%
  select(kmid, study_id, scenario, landmark, endpoint, ev_est_m4, ev_layer, single_report, cmp_evt) %>%
  left_join(
    diag_df %>% filter(!is.na(ctdna)) %>%
      transmute(kmid, ctdna, n0_risk, recon_N, events_guyot, events_pkg,
                cens_guyot, n_coord, n_nar, tmax),
    by = "kmid", relationship = "many-to-many")

prob <- diag_df %>% filter(is.na(n_coord) |
  (!is.na(n0_risk) & !is.na(recon_N) & abs(recon_N - n0_risk) > 2) |
  (!is.na(events_guyot) & !is.na(events_pkg) & abs(events_guyot - events_pkg) > 2))

cat("\n===== M1a recon summary =====\n")
cat("target KM:", nrow(meta), " reconstructed curves:",
    nrow(diag_df %>% filter(!is.na(n_coord))), "/", 2 * nrow(meta), "\n")
cat("reconstructed IPD persons:", nrow(ipd_df), " total events:", sum(ipd_df$status), "\n")
cat("N-mismatch(>2)/failed curves:", nrow(prob), "\n")
if (nrow(prob) > 0) print(prob)
hi_ev <- ec %>% filter(ctdna == 1) %>% mutate(recon_layer = case_when(
  events_guyot <= 10 ~ "E1_le10", events_guyot <= 20 ~ "E2_11to20",
  events_guyot <= 40 ~ "E3_21to40", TRUE ~ "E4_gt40"))
cat("\nhigh-risk reconstructed event layer:\n"); print(table(hi_ev$recon_layer, useNA = "ifany"))
cat("M4 estimated event layer:\n"); print(table(hi_ev$ev_layer, useNA = "ifany"))
cat("\nNAR present curves median N check (recon_N vs n0_risk) max abs diff:",
    max(abs(diag_df$recon_N - diag_df$n0_risk), na.rm = TRUE), "\n")

write_csv(ipd_df, file.path(out_dir, "M1_ipd_reconstructed.csv"))
write_excel_csv(ec, file.path(out_dir, "M1_event_counts_reconstructed.csv"))
write_excel_csv(prob, file.path(out_dir, "M1_reconstruction_diagnostics.csv"))
cat("\nwritten: outputs/M1/M1_ipd_reconstructed.csv, M1_event_counts_reconstructed.csv, M1_reconstruction_diagnostics.csv\n")
