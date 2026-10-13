# ============================================================================
# M1c_onestage_RP.R
# Module 1, step 3: one-stage Royston-Parmar flexible-spline pooling -> continuous-time absolute-risk curves
#               (Shiny calculator engine), cross-validated against M1b two-stage discrete points
#
# Method: within each scenario x window x endpoint x ctDNA cell, fit flexsurvspline on the reconstructed IPD
#       (Royston-Parmar on the log cumulative hazard scale; internal knots k adapt to event count, with fallback to lower degrees);
#       use cluster=study_id for cluster-robust SE; predict continuous (0..tmax) CIR=1-S with 95% CI.
#       no output when independent studies <3 per cell (hard rule, consistent with M1b).
#       Note: the one-stage model gives continuous smooth curves (for the calculator); M1b two-stage pooling is the primary result,
#       and the two are cross-validated at 12/24/36 months as mutual sensitivity analyses.
#
# Inputs: outputs/M1/M1_ipd_reconstructed.csv, data/meta/km_meta.csv,
#       outputs/M1/M1_absrisk_pooled.csv (for cross-validation)
# Outputs (outputs/M1/):
#   M1c_continuous_curves.csv   continuous absolute-risk curves (calculator engine)
#   M1c_onestage_vs_twostage.csv one-stage vs two-stage comparison at 12/24/36 months
#   M1c_continuous_curves.png   continuous-curve plot (white background)
# Depends on survival,flexsurv,dplyr,tidyr,readr,ggplot2
# ============================================================================
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(survival)
  library(flexsurv); library(ggplot2)
})
options(stringsAsFactors = FALSE)
out_dir <- "outputs/M1"; MIN_K <- 3; EVAL <- c(12, 24, 36)
SCEN_MAIN <- c("LARC", "MRD", "CRLM"); WIN_MAIN <- c("preop", "postopMRD", "postNAC", "postACT")

ipd  <- read_csv(file.path(out_dir, "M1_ipd_reconstructed.csv"), show_col_types = FALSE)
meta <- read_csv("data/meta/km_meta.csv", show_col_types = FALSE) %>%
  select(kmid, study_id, scen_code, win, endp_code)
ipd <- ipd %>% left_join(meta, by = c("kmid", "study_id")) %>%
  filter(endp_code == "DFS", scen_code %in% SCEN_MAIN, win %in% WIN_MAIN)

fit_curve <- function(dd) {
  kstud <- n_distinct(dd$study_id); nev <- sum(dd$status)
  if (kstud < MIN_K || nev < 8) return(NULL)
  kc <- if (nev >= 80) 2 else 1
  # continuous curves are the calculator engine: individual-level flexible splines; random-effects pooling and uncertainty use M1b two-stage as the primary result
  fit <- tryCatch(flexsurvspline(Surv(time, status) ~ 1, data = dd, scale = "hazard", k = kc),
                  error = function(e) tryCatch(flexsurvspline(Surv(time, status) ~ 1,
                    data = dd, scale = "hazard", k = 1), error = function(e2) NULL))
  if (is.null(fit)) return(NULL)
  tmax <- max(dd$time)
  grid <- seq(0.5, tmax, by = 0.5)
  ss <- summary(fit, type = "survival", t = grid, tidy = TRUE)  # columns: time, est, lcl, ucl
  data.frame(time_month = ss$time, CIR = 1 - ss$est,
             CIRlo = 1 - ss$ucl, CIRhi = 1 - ss$lcl,
             n_persons = nrow(dd), events = nev, k_studies = kstud,
             spline_df = kc, tmax = tmax)
}

cells <- ipd %>% distinct(scen_code, win, endp_code, ctdna)
cc <- list(); ii <- 0
for (i in seq_len(nrow(cells))) {
  z <- cells[i, ]
  dd <- ipd %>% filter(scen_code == z$scen_code, win == z$win,
                       endp_code == z$endp_code, ctdna == z$ctdna)
  f <- fit_curve(dd)
  if (is.null(f)) next
  ii <- ii + 1
  cc[[ii]] <- data.frame(scen_code = z$scen_code, win = z$win, endp_code = z$endp_code,
                         ctdna = z$ctdna, f, stringsAsFactors = FALSE)
}
curv <- bind_rows(cc)
write_csv(curv, file.path(out_dir, "M1c_continuous_curves.csv"))

# ---------- one-stage vs two-stage cross-validation at 12/24/36 months ----------
ts <- read_csv(file.path(out_dir, "M1_absrisk_pooled.csv"), show_col_types = FALSE)
one <- curv %>% mutate(tp = round(time_month)) %>% filter(tp %in% EVAL) %>%
  group_by(scen_code, win, endp_code, ctdna, tp) %>%
  arrange(abs(time_month - tp), .by_group = TRUE) %>% slice(1) %>% ungroup() %>%
  transmute(scen_code, win, endp_code, ctdna, tp,
            onestage = CIR, os_lo = CIRlo, os_hi = CIRhi, k_studies)
cmp <- ts %>% filter(endp_code == "DFS") %>%
  select(scen_code, win, ctdna, tp, twostage = pool, ts_lo = lo, ts_hi = hi,
         k_two = k, insufficient) %>%
  left_join(one, by = c("scen_code", "win", "ctdna", "tp")) %>%
  mutate(abs_diff = abs(onestage - twostage))
write_csv(cmp, file.path(out_dir, "M1c_onestage_vs_twostage.csv"))
cat("\n===== M1c one-stage vs two-stage (DFS, matched points) =====\n")
print(as.data.frame(cmp %>% filter(!insufficient) %>%
  mutate(across(c(twostage, onestage, abs_diff, os_lo, os_hi), ~ round(.x, 3)))), row.names = FALSE)
v <- cmp %>% filter(!insufficient) %>% pull(abs_diff)
cat("\nmatched points:", sum(!is.na(v)),
    "  mean|diff|:", round(mean(v, na.rm = TRUE), 3),
    "  max|diff|:", round(max(v, na.rm = TRUE), 3), "\n")

# ---------- continuous-curve plot (white background; positive red, negative blue)----------
lab <- c("preop" = "Pre-op baseline", "postopMRD" = "Post-op MRD window (4-12 wk)",
         "postNAC" = "Post-neoadjuvant", "postACT" = "Post adjuvant/systemic")
pltd <- curv %>% mutate(ctDNA_state = ifelse(ctdna == 1, "ctDNA+", "ctDNA-"),
                        panel = paste(scen_code, lab[win], sep = " / "))
p <- ggplot(pltd, aes(time_month, CIR * 100, color = ctDNA_state, fill = ctDNA_state)) +
  geom_ribbon(aes(ymin = CIRlo * 100, ymax = CIRhi * 100), alpha = .12, linetype = 0) +
  geom_line(linewidth = .9) +
  facet_wrap(~ panel, scales = "free_y") +
  scale_color_manual(values = c("ctDNA+" = "#C0392B", "ctDNA-" = "#2C3E50")) +
  scale_fill_manual(values = c("ctDNA+" = "#C0392B", "ctDNA-" = "#2C3E50")) +
  labs(x = "Months since surgery / enrolment", y = "Cumulative incidence of recurrence/event, CIR (%)", color = NULL, fill = NULL,
       title = "Pooled continuous absolute recurrence risk by ctDNA status (one-stage Royston-Parmar, reconstructed IPD)") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"), plot.title = element_text(hjust = .5, size = 11),
        legend.position = "bottom")
ggsave(file.path(out_dir, "M1c_continuous_curves.png"), p, width = 11, height = 8, dpi = 300)
cat("written: M1c_continuous_curves.csv/.png, M1c_onestage_vs_twostage.csv\n")
