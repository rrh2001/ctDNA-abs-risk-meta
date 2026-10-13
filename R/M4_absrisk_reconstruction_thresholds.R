# ============================================================================
# M4_absrisk_reconstruction_thresholds.R
# Module 4 (methodological safeguard): threshold simulation for reconstructing ABSOLUTE RISK from digitized KM curves
# Goal: use simulation to answer - "to reconstruct 12/24/36-month absolute recurrence risk, how many 
#       digitized points per curve and how dense a risk table are needed for reliable reconstruction?" Back-fill the A/B matrix thresholds accordingly.
#
# Why the old 36_bias_simulation_grid_v3.R cannot be reused:
#   old grid: fixed 31 equally spaced dense points + full n.risk, evaluating only regression coefficients beta_tt/beta_ctdna;
#           n=500-5000, HR<=4.5 -- this answers"whether HR(t) has systematic bias under ideal digitization"。
#   this script: (1) treats digitized points per curve k, risk-table density and digitization noise as grid variables;
#           (2) evaluates bias/MAE/RMSE of absolute risk CIR=1-S(t) at 12/24/36 months,
#             95% CI coverage/calibration; (3) extends n and ctDNA prevalence down to realistic small samples (positive groups often have only
#             10-60 patients with single-digit events); (4) covers HR 4/7/12 (realistic ctDNA range).
# Follows ADEMP; parallelized with future/furrr; run the TEST_MODE smoke test first, then the full run.
# Depends on survival, dplyr, tidyr, purrr, readr, ggplot2, future, furrr
# ============================================================================
suppressPackageStartupMessages({
  library(survival); library(dplyr); library(tidyr); library(purrr)
  library(readr); library(ggplot2); library(future); library(furrr)
})

out_dir <- "outputs/M4"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Test switch: TRUE=quick smoke (B=15, 2 cells, ~1 min); FALSE=full (B=500, 320 main cells, several hours)
TEST_MODE   <- FALSE
MASTER_SEED <- 20261012
EVAL_TIMES  <- c(12, 24, 36)          # timepoints (months) for absolute-risk evaluation
FOLLOW_UP   <- 60                     # maximum follow-up (months)
RELBIAS_TOL <- 0.10                   # |relative bias| <=10% considered acceptable
COV_LO <- 0.90; COV_HI <- 0.98        # acceptable range for 95% CI coverage

if (TEST_MODE) {
  B_CORE <- 15
  core_grid <- expand.grid(n=200, prev=0.20, k_pts=8, risk="medium_6m",
                           noise=0.02, beta0=log(7), beta_tt=0, cens=0.4,
                           stringsAsFactors=FALSE)
  confirm_grid <- core_grid[0,]
  cat("\n[TEST] B=",B_CORE,", 1 main-grid smoke cell\n\n")
} else {
  B_CORE <- 500
  # ---- main grid: under a typical effect (HR=7, no time decay, 40% censoring), vary sample x prevalence x points x risk table x noise ----
  core_grid <- expand.grid(
    n      = c(100, 200, 500, 1000),          # total sample (covers small positive groups)
    prev   = c(0.15, 0.30),                    # ctDNA prevalence
    k_pts  = c(5, 8, 12, 20, 40),              # digitized points per curve
    risk   = c("full", "medium_6m", "sparse_12m", "none"),  # risk-table density
    noise  = c(0, 0.02),                       # digitization readout noise (SD of St)
    beta0  = log(7), beta_tt = 0, cens = 0.40,
    stringsAsFactors = FALSE)
  # ---- confirmation grid: at the provisional threshold (k=12, medium_6m), verify robustness across effect size/decay/censoring ----
  confirm_grid <- expand.grid(
    n=200, prev=0.15, k_pts=12, risk="medium_6m", noise=0.02,
    beta0  = c(log(4), log(7), log(12)),
    beta_tt= c(0, -0.20),
    cens   = c(0.30, 0.50),
    stringsAsFactors = FALSE)
  cat("\n[FULL] B=",B_CORE,", main grid ",nrow(core_grid)," cells + confirmation ",nrow(confirm_grid)," cells\n\n")
}

# ----------------------------------------------------------------------------
# 1. true IPD data-generating mechanism (Weibull, time-dependent HR in the ctDNA group)
#    S0(t)=exp[-(t/scale)^shape]; cumulative hazard in the positive group is multiplied by exp(beta0 + beta_tt*log(t+1))
# ----------------------------------------------------------------------------
sim_true_ipd <- function(n, prev, beta0, beta_tt, cens, seed) {
  set.seed(seed)
  x  <- rbinom(n, 1, prev)
  U  <- runif(n)
  shape <- 1.2; scale <- 40
  te <- numeric(n)
  for (i in seq_len(n)) {
    Ht <- function(t) (t/scale)^shape * exp(beta0*x[i] + beta_tt*x[i]*log(t+1))
    if (Ht(FOLLOW_UP) <= -log(U[i])) te[i] <- FOLLOW_UP
    else te[i] <- uniroot(function(t) Ht(t)+log(U[i]), lower=0.001,
                          upper=FOLLOW_UP, tol=1e-6)$root
  }
  lam_c <- -log(1-cens)/FOLLOW_UP
  tc <- pmin(rexp(n, lam_c), FOLLOW_UP)
  time <- pmin(te, tc); status <- as.integer(te <= tc)
  data.frame(time=time, status=status, x=x)
}

# ----------------------------------------------------------------------------
# 2. Guyot IPD reconstruction (as in the prior script; events uniform within interval, others uniformly censored)
#    times/survivals = digitized points; nrisk_times/nrisk_vals = risk table
# ----------------------------------------------------------------------------
guyot <- function(times, surv, nr_t, nr_v) {
  if (length(times) < 3) return(NULL)
  if (times[1] != 0) { times <- c(0,times); surv <- c(1,surv) }
  o <- order(times); times <- times[o]; surv <- cummin(surv[o])
  surv <- pmin(pmax(surv,0),1)
  Sat <- function(t) surv[findInterval(t, times, all.inside=TRUE)]
  et <- numeric(0); es <- integer(0)
  for (i in seq_len(length(nr_t)-1)) {
    lo<-nr_t[i]; hi<-nr_t[i+1]; nlo<-nr_v[i]; nhi<-nr_v[i+1]
    if (is.na(nlo)||nlo<=0) next
    if (is.na(nhi)) nhi <- 0
    Slo<-Sat(lo); Shi<-Sat(hi); if (Slo<=0) next
    d <- round(nlo*(1-Shi/Slo)); d<-max(0,min(d,nlo)); c<-max(0,nlo-d-nhi)
    if (d>0){z<-seq(lo,hi,length.out=d+2)[-c(1,d+2)];et<-c(et,z);es<-c(es,rep(1L,d))}
    if (c>0){z<-seq(lo,hi,length.out=c+2)[-c(1,c+2)];et<-c(et,z);es<-c(es,rep(0L,c))}
  }
  if (length(et)==0) return(NULL)
  data.frame(time=et,status=es)
}

# ----------------------------------------------------------------------------
# 3. from true IPD -> digitize by (points k, risk-table scheme, noise) -> reconstruct IPD
#    risk: full=numbers at every digitized point; medium_6m=every 6 months; sparse_12m=every 12 months;
#          none=only initial N known (end at-risk estimated from the curve, representing"the no-risk-table"most permissive case)
# ----------------------------------------------------------------------------
digitize_and_reconstruct <- function(d, grp_x, k_pts, risk, noise) {
  dg <- d[d$x==grp_x, ]; N0 <- nrow(dg)
  sf <- survfit(Surv(time,status)~1, data=dg)
  tmax <- max(d$time)
  # digitization times: 0 and k equally spaced points within follow-up
  dt <- c(0, seq(0, tmax, length.out=k_pts+1)[-1]); dt <- unique(round(dt,3))
  ss <- summary(sf, times=dt, extend=TRUE)
  St <- pmin(pmax(ss$surv + rnorm(length(ss$surv),0,noise),0),1)  # digitization readout noise
  St <- cummin(St)
  # risk table
  if (risk=="full") { rt<-dt; rv<-summary(sf,times=dt,extend=TRUE)$n.risk }
  else {
    step <- if(risk=="medium_6m") 6 else if(risk=="sparse_12m") 12 else NA
    if (is.na(step)) { rt<-c(0,tmax); rv<-c(N0, max(0,round(N0*St[length(St)]))) }  # none
    else { rt<-seq(0, floor(tmax/step)*step, by=step); if(!tmax%in%rt) rt<-c(rt,tmax)
           rv<-summary(sf,times=rt,extend=TRUE)$n.risk }
  }
  guyot(dt[-1], St[-1], rt, rv)   # drop t=0 (the function adds 1 internally)
}

# ----------------------------------------------------------------------------
# 4. single iteration: return true vs reconstructed absolute risk and CI coverage for +/- groups at 12/24/36 months
# ----------------------------------------------------------------------------
one_iter <- function(n,prev,k_pts,risk,noise,beta0,beta_tt,cens,seed) {
  d <- tryCatch(sim_true_ipd(n,prev,beta0,beta_tt,cens,seed),error=function(e)NULL)
  if (is.null(d)) return(NULL)
  out <- list(); ii<-0
  for (gx in c(1,0)) {
    sf.t <- survfit(Surv(time,status)~1, data=d[d$x==gx,])
    st.t  <- summary(sf.t, times=EVAL_TIMES, extend=TRUE)
    rec <- tryCatch(digitize_and_reconstruct(d,gx,k_pts,risk,noise),error=function(e)NULL)
    for (j in seq_along(EVAL_TIMES)) {
      ii<-ii+1
      true_cir <- 1-st.t$surv[j]
      if (is.null(rec)) { out[[ii]]<-data.frame(group=gx,time=EVAL_TIMES[j],true_cir=true_cir,
                 rec_cir=NA,lcl=NA,ucl=NA,events_true=sum(d$x==gx & d$status==1),
                 n_true=sum(d$x==gx),converged=FALSE); next }
      sf.r <- survfit(Surv(time,status)~1,data=rec)
      sr <- summary(sf.r, times=EVAL_TIMES[j], extend=TRUE)
      # 95% CI on the log-log scale (survfit standard output)
      out[[ii]]<-data.frame(group=gx,time=EVAL_TIMES[j],true_cir=true_cir,
        rec_cir=1-sr$surv, lcl=1-sr$upper, ucl=1-sr$lower,
        events_true=sum(d$x==gx & d$status==1), n_true=sum(d$x==gx), converged=TRUE)
    }
  }
  bind_rows(out)
}

# ----------------------------------------------------------------------------
# 5. run one grid (parallel)
# ----------------------------------------------------------------------------
run_grid <- function(gd, B, tag) {
  cells <- vector("list", nrow(gd))
  plan(multisession, workers=max(1,parallel::detectCores()-1))
  for (g in seq_len(nrow(gd))) {
    r<-gd[g,]
    res<-future_map_dfr(seq_len(B), function(b){
      seed<-MASTER_SEED + ifelse(tag=="core",g*10000, 900000+g*10000) + b
      z<-tryCatch(one_iter(r$n,r$prev,r$k_pts,r$risk,r$noise,r$beta0,r$beta_tt,r$cens,seed),
                  error=function(e)NULL)
      if(is.null(z))return(NULL)
      z$iter<-b; cbind(z,r,row.names=NULL)
    }, .options=furrr_options(seed=TRUE, packages=c("survival","dplyr")))
    cells[[g]]<-res
    cat(sprintf("[%s] %d/%d done (n=%d prev=%.2f k=%d risk=%s noise=%.2f HR=%.1f)\n",
        tag,g,nrow(gd),r$n,r$prev,r$k_pts,r$risk,r$noise,exp(r$beta0)))
    write_csv(res, file.path(out_dir,sprintf("%s_chunk_%03d.csv",tag,g)))  # checkpoint resume
  }
  plan(sequential); bind_rows(cells)
}

cat("==== Main grid ====\n"); core_all<-run_grid(core_grid,B_CORE,"core")
cat("==== Confirmation grid ====\n"); conf_all<-run_grid(confirm_grid,ifelse(TEST_MODE,B_CORE,500),"confirm")
all<-bind_rows(core_all,conf_all)
write_csv(all,file.path(out_dir,"M4_all_iterations.csv"))

# ----------------------------------------------------------------------------
# 6. performance summary: per (cell x group x time) absolute-risk bias/MAE/RMSE/coverage
# ----------------------------------------------------------------------------
perf <- all %>% filter(converged) %>% mutate(
  err=rec_cir-true_cir, abserr=abs(rec_cir-true_cir),
  relerr=ifelse(true_cir>0.05,(rec_cir-true_cir)/true_cir,NA_real_),
  cover=(true_cir>=lcl & true_cir<=ucl)) %>%
  group_by(across(c(n,prev,k_pts,risk,noise,beta0,beta_tt,cens,group,time))) %>%
  summarise(n_valid=n(), events=mean(events_true), n_grp=mean(n_true),
    true_cir=mean(true_cir), rec_cir=mean(rec_cir),
    bias=mean(err,na.rm=TRUE), relbias=mean(relerr,na.rm=TRUE),
    mae=mean(abserr,na.rm=TRUE), rmse=sqrt(mean(err^2,na.rm=TRUE)),
    coverage=mean(cover,na.rm=TRUE),
    pass=(abs(mean(relerr,na.rm=TRUE))<=RELBIAS_TOL &
          mean(cover,na.rm=TRUE)>=COV_LO & mean(cover,na.rm=TRUE)<=COV_HI),
    .groups="drop")
write_csv(perf,file.path(out_dir,"M4_threshold_performance.csv"))

# threshold summary: pass rate by (points x risk table x noise) across n/prevalence/group/time (focus on the positive group, most clinically relevant)
thr <- perf %>% filter(group==1) %>%
  group_by(k_pts,risk,noise) %>%
  summarise(n_cells=n(), n_cells_pass=sum(pass,na.rm=TRUE),
            pass_rate=mean(pass,na.rm=TRUE), median_MAE=median(mae),
            median_coverage=median(coverage), .groups="drop") %>%
  arrange(noise,risk,k_pts)
write_csv(thr,file.path(out_dir,"M4_threshold_summary.csv"))

# event-stratified thresholds (clinical decisions use"positive-group event count"rather than total n)
thr_e <- perf %>% filter(group==1) %>%
  mutate(ev_stratum=cut(events,breaks=c(-1,5,10,20,40,Inf),
                     labels=c("≤5","6-10","11-20","21-40",">40"))) %>%
  group_by(ev_stratum,k_pts,risk) %>%
  summarise(pass_rate=mean(pass,na.rm=TRUE),median_MAE=median(mae),
            median_relbias=median(relbias,na.rm=TRUE),.groups="drop")
write_csv(thr_e,file.path(out_dir,"M4_threshold_by_events.csv"))

# ----------------------------------------------------------------------------
# 7. figure: points x risk-table pass-rate heatmap (by noise; positive group)
# ----------------------------------------------------------------------------
if (!TEST_MODE) {
  p<-ggplot(perf%>%filter(group==1),aes(factor(k_pts),factor(risk,levels=c("full","medium_6m","sparse_12m","none")),fill=mae))+
    geom_tile(color="white")+geom_text(aes(label=sprintf("%.2f",mae)),size=3)+
    facet_grid(noise~cut(events,breaks=c(-1,10,20,Inf),labels=c("events<=10","11-20",">20")))+
    scale_fill_gradient2(low="#2C7BB6",mid="#FFFFBF",high="#D7191C",midpoint=0.10,limits=c(0,0.25),name="MAE\n(absolute risk)")+
    labs(x="digitized points per curve",y="risk-table density",
         title="Absolute-risk reconstruction error (ctDNA+ group)")+
    theme_minimal(base_size=11)
  ggsave(file.path(out_dir,"M4_heatmap_MAE.png"),p,width=10,height=6,dpi=300)

  pc<-ggplot(perf%>%filter(group==1),aes(factor(k_pts),factor(risk,levels=c("full","medium_6m","sparse_12m","none")),fill=coverage))+
    geom_tile(color="white")+geom_text(aes(label=sprintf("%.2f",coverage)),size=3)+
    facet_grid(noise~cut(events,breaks=c(-1,10,20,Inf),labels=c("events<=10","11-20",">20")))+
    scale_fill_gradient2(low="#D7191C",mid="#FFFFBF",high="#2C7BB6",midpoint=0.95,limits=c(0.75,1),name="95%CI\ncoverage")+
    labs(x="digitized points per curve",y="risk-table density",title="95%CI coverage of reconstructed absolute risk (ctDNA+)")+
    theme_minimal(base_size=11)
  ggsave(file.path(out_dir,"M4_heatmap_coverage.png"),pc,width=10,height=6,dpi=300)
}

# ----------------------------------------------------------------------------
# 8. automatic threshold decision (print + write suggested A/B criteria)
# ----------------------------------------------------------------------------
cat("\n================ M4 thresholds (positive group, across timepoints)================\n")
print(as.data.frame(thr))
# find: among cells with events>=10, the most permissive (fewest points / sparsest risk table) combination with pass rate>=90% -> minimum grade-A configuration
cand <- perf %>% filter(group==1, events>=10, noise==0.02) %>%
  group_by(k_pts,risk) %>% summarise(pass_rate=mean(pass,na.rm=TRUE),.groups="drop") %>%
  filter(pass_rate>=0.90) %>% arrange(k_pts, desc(risk))
rec_A <- if(nrow(cand)>0) cand[1,] else NULL
sink(file.path(out_dir,"M4_recommended_thresholds.txt"))
cat("M4 reconstructed absolute-risk thresholds (confirm after the full simulation)\n\n")
cat("Grade A (pool absolute risk directly) suggested minimum (events>=10, pass rate>=90% with 0.02 digitization noise):\n")
if(!is.null(rec_A)) print(as.data.frame(rec_A)) else cat("No combination passes in the full result; relax the event count or increase points/risk-table requirements; see M4_threshold_by_events.csv\n")
cat("\nGrade B (sensitivity/narrative only): between A and C; Grade C (HR only, no absolute risk): no digitized curve.\n")
cat("\nThe confirmation grid (HR 4/7/12 x decay x censoring) tests threshold generalizability; see M4_threshold_performance.csv.\n")
sink()
cat("\nDone. Results in ",normalizePath(out_dir),"\n")
