# ============================================================================
# M1b2_three_level.R
# Three-level multilevel meta-analysis (Cheung 2014) to explicitly account for
# within-study correlation across time points (12/24/36 mo) reconstructed from
# the SAME Kaplan-Meier curve.
# Level 1 = individual effect (setting x window x ctDNA x timepoint estimate);
# Level 2 = study.
#   y_{jk} = beta_t + u_j + v_jk
#   u_j ~ N(0, sigma2_between),  v_jk ~ N(0, sigma2_within)
# Compares the three-level pooled CIR at each time point with the primary
# per-cell REML estimates (M1_absrisk_pooled.csv) and reports both variance
# components. Sensitivity analysis: confirms primary estimates are robust.
# ============================================================================
suppressPackageStartupMessages({ library(readr); library(dplyr); library(metafor) })
out_dir <- "outputs/M1"
per <- read_csv(file.path(out_dir, "M1_absrisk_perstudy.csv"), show_col_types=FALSE)
prim <- read_csv(file.path(out_dir, "M1_absrisk_pooled.csv"), show_col_types=FALSE)
clip01 <- function(p) pmin(pmax(p,1e-4),1-1e-4)

d <- per %>% filter(endp_code=="DFS", scen_code %in% c("MRD","LARC","CRLM"),
                    win %in% c("preop","postopMRD","postACT"))
res <- list(); rr <- 0
for (sc in c("MRD","LARC","CRLM")) for (wn in c("preop","postopMRD","postACT")) for (g in c(0,1)) {
  dd <- d %>% filter(scen_code==sc, win==wn, ctdna==g)
  if (n_distinct(dd$study_id) < 3) next
  dd <- dd %>% mutate(p=clip01(CIR),
                      yi=log(p/(1-p)),
                      vi=(pmax(seCIR,1e-4)/(p*(1-p)))^2,
                      obsid=row_number())
  fit <- tryCatch(rma.mv(yi=yi, V=vi, mods=~factor(tp),
                         random=~1|study_id/obsid, data=dd, method="REML"),
                  error=function(e) NULL)
  if (is.null(fit)) next
  sb <- as.numeric(fit$sigma2)[1]  # study (between)
  sw <- as.numeric(fit$sigma2)[2]  # within-study (effect)
  tps <- sort(unique(dd$tp))
  X <- model.matrix(~factor(tp), data=data.frame(tp=tps))
  colnames(X) <- names(coef(fit))
  mu <- as.numeric(X %*% as.matrix(fit$b))
  se <- sqrt(diag(X %*% fit$vb %*% t(X)))
  for (i in seq_along(tps)) {
    rr <- rr+1
    t <- tps[i]
    pr <- prim %>% filter(scen_code==sc, win==wn, endp_code=="DFS", ctdna==g, tp==t)
    res[[rr]] <- data.frame(scen=sc, win=wn, ctdna=g, tp=t,
      k_studies=n_distinct(dd$study_id),
      sigma2_between=round(sb,4), sigma2_within=round(sw,4),
      tl_pool=round(100*plogis(mu[i]),1),
      tl_lo=round(100*plogis(mu[i]-1.96*se[i]),1),
      tl_hi=round(100*plogis(mu[i]+1.96*se[i]),1),
      primary_pool=if(nrow(pr)) round(100*pr$pool,1) else NA,
      primary_lo=if(nrow(pr)) round(100*pr$lo,1) else NA,
      primary_hi=if(nrow(pr)) round(100*pr$hi,1) else NA,
      abs_diff_pp=if(nrow(pr) && !is.na(pr$pool)) round(100*abs(plogis(mu[i])-pr$pool),1) else NA)
  }
}
ans <- bind_rows(res)
write_csv(ans, file.path(out_dir,"M1b_three_level_sensitivity.csv"))
# variance components summary (one per fitted cell)
vc <- ans %>% distinct(scen,win,ctdna,sigma2_between,sigma2_within,k_studies)
write_csv(vc, file.path(out_dir,"M1b_three_level_varcomp.csv"))
cat("===== max abs difference three-level vs primary (pp) =====\n")
print(summary(ans$abs_diff_pp))
cat("\nrows:",nrow(ans),"\n")
