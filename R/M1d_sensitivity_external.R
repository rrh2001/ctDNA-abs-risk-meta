# ============================================================================
# M1d_sensitivity_external.R
# Module 1 robustness and external validity:
#   P1 leave-one-study-out (stability of pooled RD)
#   P2 publication bias (Egger test, cells with k>=5)
#   P3 reconstructed HR vs published values (log-log approximation, agreement)
#   P4 external calibration on Li2022 real IPD (postACT scenario); Ogaard kept as an advanced-disease external anchor
# Inputs: outputs/M1/M1_absrisk_perstudy.csv, data/meta/km_meta.csv,
#       working dir ipd_clean_Li2022_postchemo.csv / ipd_clean_Ogaard2022_L1_postOP.csv / _L3_postACT.csv
# Outputs: outputs/M1/S4_leave_one_study_out.csv, S4_egger_small_study.csv,
#       S3_reconstructed_vs_published_HR.csv, S5_external_calibration.csv
# Depends on survival,metafor,dplyr,tidyr,readr
# ============================================================================
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(survival); library(metafor)
})
options(stringsAsFactors = FALSE)
OUT <- "outputs/M1"; EVAL <- c(12, 24, 36)

per <- read_csv(file.path(OUT, "M1_absrisk_perstudy.csv"), show_col_types = FALSE)
# pair per kmid
w <- per %>% select(kmid, study_id, scen_code, win, endp_code, tp, ctdna, CIR, seCIR) %>%
  pivot_wider(names_from = ctdna, values_from = c(CIR, seCIR), names_sep = "") %>%
  filter(!is.na(CIR1), !is.na(CIR0)) %>%
  mutate(RD = CIR1 - CIR0, vRD = seCIR1^2 + seCIR0^2,
         # logHR on log-log scale (reconstructed)
         lHR = log(-log(pmax(1 - CIR1, 1e-4))) - log(-log(pmax(1 - CIR0, 1e-4))),
         vHR = (seCIR1 / (-(1 - CIR1) * log(pmax(1 - CIR1, 1e-4))))^2 +
               (seCIR0 / (-(1 - CIR0) * log(pmax(1 - CIR0, 1e-4))))^2)

# ---------------- P1 leave-one-study-out ----------------
loo <- list(); a <- 0
for (gidx in (w %>% distinct(scen_code, win, endp_code, tp) %>% row.names())) {
  key <- w %>% distinct(scen_code, win, endp_code, tp)
  for (r in seq_len(nrow(key))) {
    k <- key[r,]
    dd <- w %>% filter(scen_code == k$scen_code, win == k$win, endp_code == k$endp_code, tp == k$tp) %>%
      filter(is.finite(vRD), vRD > 0)
    if (n_distinct(dd$study_id) < 4) next
    full <- rma(yi = dd$RD, vi = dd$vRD, method = "REML")
    for (s in unique(dd$study_id)) {
      d2 <- dd %>% filter(study_id != s) %>% filter(is.finite(vRD), vRD > 0)
      if (n_distinct(d2$study_id) < 3) next
      f2 <- rma(yi = d2$RD, vi = d2$vRD, method = "REML")
      a <- a + 1
      loo[[a]] <- data.frame(scen_code = k$scen_code, win = k$win, endp_code = k$endp_code,
        tp = k$tp, dropped = s, k_remaining = n_distinct(d2$study_id),
        RD_full = as.numeric(coef(full)), RD_loo = as.numeric(coef(f2)),
        delta = as.numeric(coef(f2) - coef(full)))
    }
  }
}
loo <- bind_rows(loo); write_csv(loo, file.path(OUT, "S4_leave_one_study_out.csv"))
cat("P1 leave-one-out rows:", nrow(loo),
    " max|delta RD|:", if (nrow(loo)) round(max(abs(loo$delta)), 3), "\n")

# ---------------- P2 Egger publication bias ----------------
eg <- list(); e <- 0
for (r in seq_len(nrow(w %>% distinct(scen_code, win, endp_code, tp)))) {
  k <- (w %>% distinct(scen_code, win, endp_code, tp))[r,]
  dd <- w %>% filter(scen_code == k$scen_code, win == k$win, endp_code == k$endp_code, tp == k$tp) %>%
    filter(is.finite(vRD), vRD > 0)
  if (n_distinct(dd$study_id) < 5) next
  fit <- rma(yi = dd$RD, vi = dd$vRD, method = "REML")
  rt <- regtest(fit)
  e <- e + 1
  eg[[e]] <- data.frame(scen_code = k$scen_code, win = k$win, endp_code = k$endp_code,
    tp = k$tp, k = n_distinct(dd$study_id),
    Egger_z = as.numeric(rt$zval), Egger_p = as.numeric(rt$pval))
}
eg <- bind_rows(eg); write_csv(eg, file.path(OUT, "S4_egger_small_study.csv"))
cat("P2 Egger cells:", nrow(eg), "\n"); print(as.data.frame(eg), row.names = FALSE)

# ---------------- P3 reconstructed HR vs published values ----------------
km <- read_csv("data/meta/km_meta.csv", show_col_types = FALSE) %>%
  select(kmid, author, year, hr, hr_type)
hrec <- w %>% mutate(hr_reconstructed = exp(lHR), hr_lcl = exp(lHR - 1.96 * sqrt(pmax(vHR, 0))),
                     hr_ucl = exp(lHR + 1.96 * sqrt(pmax(vHR, 0)))) %>%
  group_by(kmid) %>% arrange(desc(tp)) %>% slice(1) %>% ungroup() %>%
  left_join(km, by = "kmid") %>%
  mutate(hr_published = suppressWarnings(as.numeric(hr)),
         within_CI = ifelse(!is.na(hr_published),
                            hr_published >= hr_lcl & hr_published <= hr_ucl, NA))
write_csv(hrec, file.path(OUT, "S3_reconstructed_vs_published_HR.csv"))
hc <- hrec %>% filter(is.finite(hr_published), hr_published > 0, is.finite(lHR))
cat("P3 HR pairs with published value:", nrow(hc),
    " within reconstructed CI:", sum(hc$within_CI, na.rm = TRUE),
    " log-HR cor:", round(cor(log(hc$hr_published), hc$lHR), 3), "\n")

# ---------------- P4 external calibration on real IPD ----------------
true_cir <- function(file, tcol, evalt) {
  d <- read_csv(file, show_col_types = FALSE)
  names(d)[names(d) == tcol] <- "t"; names(d)[names(d) == "ctDNA"] <- "ctdna"
  d <- d %>% filter(!is.na(t), !is.na(ctdna))
  out <- list()
  for (g in c(0, 1)) {
    dg <- d %>% filter(ctdna == g)
    sf <- survfit(Surv(t, event) ~ 1, data = dg)
    ss <- summary(sf, times = evalt, extend = TRUE)
    for (j in seq_along(evalt))
      if (evalt[j] <= max(dg$t))
        out[[length(out) + 1]] <- data.frame(ctdna = g, tp = evalt[j],
                                             observed = 1 - ss$surv[j])
  }
  bind_rows(out)
}
# Third-party published appendix IPD are NOT redistributed; P4 runs only when
# the cleaned files are present locally (see README, Data availability). The
# released S5 table is otherwise retained unchanged.
ipd_files <- c("ipd_clean_Li2022_postchemo.csv",
               "ipd_clean_Ogaard2022_L1_postOP.csv",
               "ipd_clean_Ogaard2022_L3_postACT.csv")
if (all(file.exists(ipd_files))) {
  li <- true_cir(ipd_files[1], "time_month", c(12, 24)) %>%
    mutate(dataset = "Li2022_postchemo_stageIII", mapped = "MRD postACT")
  og1 <- true_cir(ipd_files[2], "fu_from_lm", c(12, 24, 36)) %>%
    mutate(dataset = "Ogaard2022_L1_postOP_mCRC", mapped = "advanced_anchor_excluded_from_pool")
  og3 <- true_cir(ipd_files[3], "fu_from_lm", c(12, 24)) %>%
    mutate(dataset = "Ogaard2022_L3_postACT_mCRC", mapped = "advanced_anchor_excluded_from_pool")
  pool <- read_csv(file.path(OUT, "M1_absrisk_pooled.csv"), show_col_types = FALSE)
  pred <- pool %>% filter(endp_code == "DFS", scen_code == "MRD", win == "postACT") %>%
    transmute(mapped = "MRD postACT", ctdna, tp, predicted = pool, pred_lo = lo, pred_hi = hi)
  ext <- bind_rows(li, og1, og3) %>%
    left_join(pred, by = c("mapped", "ctdna", "tp")) %>%
    mutate(abs_diff = abs(observed - predicted))
  write_csv(ext, file.path(OUT, "S5_external_calibration.csv"))
  cat("\nP4 external calibration:\n"); print(as.data.frame(ext), row.names = FALSE)
} else {
  cat("\nP4 skipped: third-party published appendix IPD not present in the working\n",
      "directory; the released outputs/M1/S5_external_calibration.csv is retained\n",
      "(see README, Data availability).\n")
}
cat("\nM1d done.\n")
