# ============================================================================
# M3_timing.R
# Use the M1c continuous absolute-risk curves to quantify "how many post-operative months concentrate recurrence" (dense-surveillance window definition),
# combined with published lead time (median ~7-10 months) to support the surveillance-interval framework.
# Output outputs/M3/M3_risk_timing.csv
# ============================================================================
suppressPackageStartupMessages({library(readr);library(dplyr)})
dir.create("outputs/M3",showWarnings=FALSE,recursive=TRUE)
cc <- read_csv("outputs/M1/M1c_continuous_curves.csv",show_col_types=FALSE)

cells <- cc %>% distinct(scen_code,win,endp_code,ctdna)
res <- list(); ii<-0
for (i in seq_len(nrow(cells))) {
  r<-cells[i,]; g<-cc%>%filter(scen_code==r$scen_code,win==r$win,endp_code==r$endp_code,ctdna==r$ctdna)
  tmax<-max(g$time_month)
  at<-function(tp){ if(tp>tmax) return(NA_real_); z<-g[which.min(abs(g$time_month-tp)),]; z$CIR }
  c6<-at(6); c12<-at(12); c18<-at(18); c24<-at(24); c36<-at(36)
  cref<-if(!is.na(c36)) c36 else c24
  ii<-ii+1
  res[[ii]]<-data.frame(scen_code=r$scen_code,win=r$win,endp_code=r$endp_code,ctdna=r$ctdna,
    CIR6=round(100*c6,1),CIR12=round(100*c12,1),CIR18=round(100*c18,1),
    CIR24=round(100*c24,1),CIR36=round(100*c36,1),
    pct_events_by12=round(100*c12/cref,0),pct_events_by24=round(100*c24/cref,0),
    tmax=round(tmax,1))
}
res<-bind_rows(res)
write_csv(res,"outputs/M3/M3_risk_timing.csv")
cat("done\n")
