# ============================================================================
# M1b_two_stage_absrisk.R
# Module 1, step 2: two-stage random-effects pooling -> time-resolved absolute recurrence-risk map
#
# Approach:
#   Stage 1: for each reconstructed KM (kmid) x ctDNA state, fit a KM on the reconstructed IPD and take CIR=1-S at 12/24/36 months with Greenwood SE;
#   Stage 2: random-effects (REML) pooling on the logit(CIR) scale, by scenario (scen_code) x window (win) x endpoint (endp_code)
#          x ctDNA x timepoint; no point estimate when independent studies (study_id) in a cell are <3 (hard rule).
#   Contrast: within-study paired ctDNA-positive vs negative: absolute risk difference RD (primary) and relative risk RR.
#   Window postopMRD = L1 (post-op <=1 month) + L2 (post-op 1-3 months), combined into the post-op 4-12 week MRD window;
#         if multiple KMs share a window/study, keep the main KM with the largest sample to avoid double counting.
#
# Inputs: outputs/M1/M1_ipd_reconstructed.csv, data/meta/km_meta.csv
# Outputs (outputs/M1/):
#   M1_absrisk_perstudy.csv   CIR/SE/n/events/tmax per KM x state x timepoint
#   M1_absrisk_pooled.csv     pooled absolute risk (positive/negative) with I2/tau2/prediction interval/insufficiency flag
#   M1_contrast_pooled.csv    pooled RD and RR (positive vs negative)
#   M1_quick_absrisk.png      quick faceted plot (publication figure refined later)
# Depends on survival,metafor,dplyr,tidyr,readr,ggplot2
# ============================================================================
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(survival); library(metafor); library(ggplot2)
})
options(stringsAsFactors = FALSE)
out_dir <- "outputs/M1"; EVAL <- c(12, 24, 36); MIN_K <- 3

ipd  <- read_csv(file.path(out_dir, "M1_ipd_reconstructed.csv"), show_col_types = FALSE)
meta <- read_csv("data/meta/km_meta.csv", show_col_types = FALSE) %>%
  select(kmid, study_id, scen_code, win, endp_code)
ipd <- ipd %>% left_join(meta, by = c("kmid", "study_id"))
clip01 <- function(p) pmin(pmax(p, 1e-4), 1 - 1e-4)

# ---------- Stage 1: CIR per KM x state ----------
keys <- ipd %>% distinct(kmid, study_id, scen_code, win, endp_code, ctdna)
ps <- list(); jj <- 0
for (i in seq_len(nrow(keys))) {
  k <- keys[i, ]
  dd <- ipd %>% filter(kmid == k$kmid, ctdna == k$ctdna)
  if (nrow(dd) < 5) next
  sf <- survfit(Surv(time, status) ~ 1, data = dd)
  tmax <- max(dd$time)
  ss <- summary(sf, times = EVAL, extend = TRUE)
  for (j in seq_along(EVAL)) {
    if (EVAL[j] > tmax) next                       # insufficient follow-up, do not extrapolate
    jj <- jj + 1
    ps[[jj]] <- data.frame(kmid = k$kmid, study_id = k$study_id, scen_code = k$scen_code,
      win = k$win, endp_code = k$endp_code, ctdna = k$ctdna, tp = EVAL[j],
      n = nrow(dd), events = sum(dd$status), tmax = tmax,
      CIR = 1 - ss$surv[j], seCIR = ss$std.err[j])
  }
}
per <- bind_rows(ps)
# multiple KMs per window/study: keep the KM with largest n for each study x window x endpoint x state x timepoint, avoiding double counting
per <- per %>% group_by(study_id, scen_code, win, endp_code, ctdna, tp) %>%
  arrange(desc(n), .by_group = TRUE) %>% slice(1) %>% ungroup()
write_csv(per, file.path(out_dir, "M1_absrisk_perstudy.csv"))

# ---------- Stage 2: logit random-effects pooling of absolute risk ----------
pool_one <- function(d) {
  k <- n_distinct(d$study_id)
  n_tot <- sum(d$n); ev_tot <- sum(d$events)
  if (k < MIN_K)
    return(data.frame(k = k, n_tot = n_tot, events_tot = ev_tot,
                      pool = NA, lo = NA, hi = NA, pi_lo = NA, pi_hi = NA,
                      i2 = NA, tau2 = NA, insufficient = TRUE))
  p <- clip01(d$CIR); se <- pmax(d$seCIR, 1e-4)
  yi <- log(p / (1 - p)); vi <- (se / (p * (1 - p)))^2
  fit <- tryCatch(rma(yi = yi, vi = vi, method = "REML"),
                  error = function(e) tryCatch(rma(yi = yi, vi = vi, method = "DL"), error = function(e2) NULL))
  if (is.null(fit))
    return(data.frame(k = k, n_tot = n_tot, events_tot = ev_tot,
                      pool = NA, lo = NA, hi = NA, pi_lo = NA, pi_hi = NA,
                      i2 = NA, tau2 = NA, insufficient = TRUE))
  b <- as.numeric(fit$b); vb <- as.numeric(fit$vb)
  z <- qnorm(.975); pi <- if (k >= 4) qt(.975, k - 2) * sqrt(as.numeric(fit$tau2) + vb) else NA
  data.frame(k = k, n_tot = n_tot, events_tot = ev_tot,
             pool = plogis(b), lo = plogis(b - z * sqrt(vb)), hi = plogis(b + z * sqrt(vb)),
             pi_lo = if (!is.na(pi)) plogis(b - pi) else NA,
             pi_hi = if (!is.na(pi)) plogis(b + pi) else NA,
             i2 = as.numeric(fit$I2), tau2 = as.numeric(fit$tau2), insufficient = FALSE)
}
pooled <- per %>% group_by(scen_code, win, endp_code, ctdna, tp) %>%
  group_modify(~ pool_one(.x)) %>% ungroup()
write_csv(pooled, file.path(out_dir, "M1_absrisk_pooled.csv"))

# ---------- positive vs negative paired contrast: RD (primary) and RR ----------
w <- per %>% select(kmid, study_id, scen_code, win, endp_code, tp, ctdna, CIR, seCIR, n) %>%
  pivot_wider(names_from = ctdna, values_from = c(CIR, seCIR, n), names_sep = "")
con <- list(); cc <- 0
for (i in seq_len(nrow(w))) {
  r <- w[i,]
  if (is.na(r$CIR1) || is.na(r$CIR0) || is.na(r$n1) || is.na(r$n0)) next
  s1 <- max(r$seCIR1, 1e-4); s0 <- max(r$seCIR0, 1e-4)
  p1c <- (r$CIR1 * r$n1 + 0.5) / (r$n1 + 1)          # Gart 0.5 continuity correction to prevent unstable RR in sparse groups
  p0c <- (r$CIR0 * r$n0 + 0.5) / (r$n0 + 1)
  cc <- cc + 1
  con[[cc]] <- data.frame(study_id = r$study_id, scen_code = r$scen_code, win = r$win,
    endp_code = r$endp_code, tp = r$tp,
    RD = r$CIR1 - r$CIR0, vRD = s1^2 + s0^2,
    lRR = log(p1c / p0c), vRR = (s1 / p1c)^2 + (s0 / p0c)^2)
}
con <- bind_rows(con)
con_one <- function(d) {
  k <- n_distinct(d$study_id)
  if (k < MIN_K) return(data.frame(k = k, RD = NA, RDlo = NA, RDhi = NA, RR = NA, RRlo = NA, RRhi = NA, insufficient = TRUE))
  fRD <- tryCatch(rma(yi = d$RD, vi = d$vRD, method = "REML"), error = function(e) NULL)
  fRR <- tryCatch(rma(yi = d$lRR, vi = d$vRR, method = "REML"), error = function(e) NULL)
  z <- qnorm(.975)
  data.frame(k = k,
    RD = if (!is.null(fRD)) as.numeric(fRD$b) else NA,
    RDlo = if (!is.null(fRD)) as.numeric(fRD$b) - z * sqrt(as.numeric(fRD$vb)) else NA,
    RDhi = if (!is.null(fRD)) as.numeric(fRD$b) + z * sqrt(as.numeric(fRD$vb)) else NA,
    RR = if (!is.null(fRR)) exp(as.numeric(fRR$b)) else NA,
    RRlo = if (!is.null(fRR)) exp(as.numeric(fRR$b) - z * sqrt(as.numeric(fRR$vb))) else NA,
    RRhi = if (!is.null(fRR)) exp(as.numeric(fRR$b) + z * sqrt(as.numeric(fRR$vb))) else NA,
    insufficient = FALSE)
}
contrast <- con %>% group_by(scen_code, win, endp_code, tp) %>%
  group_modify(~ con_one(.x)) %>% ungroup()
write_csv(contrast, file.path(out_dir, "M1_contrast_pooled.csv"))

# ---------- main table (DFS recurrence/disease-free endpoints; three curative-intent scenarios x main windows x timepoints)----------
main <- pooled %>% filter(endp_code == "DFS", scen_code %in% c("LARC", "MRD", "CRLM"),
                          win %in% c("preop", "postopMRD", "postACT", "postNAC")) %>%
  arrange(scen_code, win, tp, ctdna)
cat("\n===== M1b pooled absolute risk (DFS), pos=1 / neg=0 =====\n")
print(as.data.frame(main %>% mutate(across(c(pool, lo, hi, pi_lo, pi_hi), ~ round(.x, 3)))), row.names = FALSE)
cat("\n===== M1b ctDNA pos-vs-neg contrast (DFS) =====\n")
cm <- contrast %>% filter(endp_code == "DFS", scen_code %in% c("LARC", "MRD", "CRLM"),
                          win %in% c("preop", "postopMRD", "postACT", "postNAC")) %>% arrange(scen_code, win, tp)
print(as.data.frame(cm %>% mutate(across(c(RD, RDlo, RDhi, RR, RRlo, RRhi), ~ round(.x, 3)))), row.names = FALSE)

# quick faceted plot (publication figure refined later; white background)
plt <- main %>% mutate(ctDNA_state = ifelse(ctdna == 1, "ctDNA+", "ctDNA-"),
                       grp = paste(scen_code, win, sep = " / "))
p <- ggplot(plt, aes(factor(tp), pool * 100, fill = ctDNA_state)) +
  geom_col(position = position_dodge(.8), width = .7) +
  geom_errorbar(aes(ymin = lo * 100, ymax = hi * 100), position = position_dodge(.8), width = .2) +
  facet_wrap(~ grp, scales = "free_y") +
  scale_fill_manual(values = c("ctDNA+" = "#C0392B", "ctDNA-" = "#2C3E50")) +
  labs(x = "Time (months)", y = "Cumulative incidence of recurrence/event, CIR (%)", fill = NULL,
       title = "Pooled absolute recurrence risk by ctDNA status (DFS, two-stage random-effects on reconstructed IPD)") +
  theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = .5, size = 11))
ggsave(file.path(out_dir, "M1_quick_absrisk.png"), p, width = 11, height = 8, dpi = 300)
cat("\nwritten: M1_absrisk_perstudy.csv / M1_absrisk_pooled.csv / M1_contrast_pooled.csv / M1_quick_absrisk.png\n")
