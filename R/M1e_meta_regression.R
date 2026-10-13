# ============================================================================
# M1e_meta_regression.R
# Module 1 sources of heterogeneity: study-level meta-regression / subgroups (post-op MRD window DFS, 12-month primary, most studies)
#   Outcome: per-study positive vs negative RD (absolute risk difference)
#   Moderators: informed (ctDNA result disclosed), metastatic (metastatic vs non-metastatic), multicenter, prospective,
#             followup_m (median follow-up, months, centred)
#   per-variable univariable + joint (methodological-quality) models; report beta, p, tau2, R2; subgroup pooled RD
# Inputs: outputs/M1/M1_absrisk_perstudy.csv, data/meta/study_moderators.csv
# Outputs: outputs/M1/S4_meta_regression.csv, S4_subgroup_RD.csv
# Depends on metafor,dplyr,tidyr,readr
# ============================================================================
suppressPackageStartupMessages({library(readr);library(dplyr);library(tidyr);library(metafor)})
options(stringsAsFactors=FALSE)
OUT <- "outputs/M1"

per <- read_csv(file.path(OUT,"M1_absrisk_perstudy.csv"), show_col_types=FALSE)
mod <- read_csv("data/meta/study_moderators.csv", show_col_types=FALSE)

# study-level paired RD; keep the main KM (largest n) for multiple KMs in the same study/cell, avoiding double counting
w <- per %>% select(kmid,study_id,scen_code,win,endp_code,tp,ctdna,CIR,seCIR,n) %>%
  pivot_wider(names_from=ctdna, values_from=c(CIR,seCIR,n), names_sep="") %>%
  filter(!is.na(CIR1),!is.na(CIR0)) %>%
  mutate(RD=CIR1-CIR0, vRD=seCIR1^2+seCIR0^2,
         nmax=pmax(n1,n0,na.rm=TRUE)) %>%
  filter(is.finite(vRD),vRD>0) %>%
  group_by(scen_code,win,endp_code,tp,study_id) %>%
  arrange(desc(nmax),.by_group=TRUE) %>% slice(1) %>% ungroup() %>%
  left_join(mod, by="study_id")

ana <- w %>% filter(endp_code=="DFS", win=="postopMRD", tp==12)
cat("postopMRD 12m independent studies:", n_distinct(ana$study_id), "\n")
print(table(ana$scen_code)); print(table(ana$informed, useNA="ifany"))

run_models <- function(d, label){
  base <- rma(yi=d$RD, vi=d$vRD, method="REML")
  tau0 <- as.numeric(base$tau2)
  res <- list()
  covs <- c("informed","metastatic","multicenter","prospective","followup_c")
  d$followup_c <- scale(d$followup_m)[,1]
  for(cv in covs){
    dd <- d %>% filter(!is.na(.data[[cv]]))
    if(n_distinct(dd[[cv]])<2 || nrow(dd)<5) next
    f <- rma(yi=dd$RD, vi=dd$vRD, mods=as.formula(paste("~",cv)), data=dd, method="REML")
    b <- coef(f)[2]; ci <- f$ci.lb[2]; cu <- f$ci.ub[2]; p <- f$pval[2]
    R2 <- 100*max(0,(tau0-as.numeric(f$tau2))/max(tau0,1e-6))
    res[[length(res)+1]] <- data.frame(window=label, moderator=cv, k=nrow(dd),
      beta=as.numeric(b), ci_lb=as.numeric(ci), ci_ub=as.numeric(cu), p=as.numeric(p),
      tau2=as.numeric(f$tau2), R2_percent=R2)
  }
  # joint methodological-quality model (not collinear with scenario)
  dd <- d %>% filter(!is.na(informed),!is.na(multicenter),!is.na(prospective))
  if(nrow(dd)>=6){
    fj <- rma(yi=dd$RD, vi=dd$vRD, mods=~informed+multicenter+prospective, data=dd, method="REML")
    for(j in 2:length(coef(fj)))
      res[[length(res)+1]] <- data.frame(window=label,
        moderator=paste0("joint:",names(coef(fj))[j]), k=nrow(dd),
        beta=as.numeric(coef(fj)[j]), ci_lb=as.numeric(fj$ci.lb[j]),
        ci_ub=as.numeric(fj$ci.ub[j]), p=as.numeric(fj$pval[j]),
        tau2=as.numeric(fj$tau2),
        R2_percent=100*max(0,(tau0-as.numeric(fj$tau2))/max(tau0,1e-6)))
  }
  bind_rows(res)
}
`%||%` <- function(a,b) if(is.null(a)) b else a

r1 <- run_models(ana, "postopMRD_12m")
a2 <- w %>% filter(endp_code=="DFS", win=="postopMRD", tp==24)
r2 <- if(n_distinct(a2$study_id)>=6) run_models(a2,"postopMRD_24m") else data.frame()
mr <- bind_rows(r1,r2)
write_csv(mr, file.path(OUT,"S4_meta_regression.csv"))
cat("\n===== meta regression (RD, pp) beta*100 =====\n")
print(as.data.frame(mr %>% mutate(across(c(beta,ci_lb,ci_ub),~round(.x*100,1)))), row.names=FALSE)

# subgroup pooled RD (by assay / metastatic status)
sub <- list(); s<-0
for(cv in c("informed","metastatic")){
  for(lvl in na.omit(unique(ana[[cv]]))){
    dd <- ana %>% filter(.data[[cv]]==lvl)
    if(n_distinct(dd$study_id)<3) next
    f <- rma(yi=dd$RD, vi=dd$vRD, method="REML")
    s<-s+1
    sub[[s]] <- data.frame(moderator=cv, level=lvl, k=n_distinct(dd$study_id),
      RD=as.numeric(coef(f)), lo=f$ci.lb, hi=f$ci.ub, I2=f$I2)
  }
}
sub <- bind_rows(sub); write_csv(sub, file.path(OUT,"S4_subgroup_RD.csv"))
cat("\n===== subgroup pooled RD (postopMRD 12m) =====\n")
print(as.data.frame(sub %>% mutate(across(c(RD,lo,hi),~round(.x*100,1)))), row.names=FALSE)
cat("\nM1e done.\n")
