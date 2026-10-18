# Leave-one-out sensitivity: drop Lygre (study_id 50) and re-pool cells it enters
suppressPackageStartupMessages({library(readr);library(dplyr);library(metafor)})
setwd(r'{D:\研究生学习\博士-李杰老师课题组\论文\1、ctDNA预测结直肠癌复发效应的时间衰减轨迹\ctDNA-abs-risk-meta}')
per <- read_csv('outputs/M1/M1_absrisk_perstudy.csv', show_col_types=FALSE)
clip01 <- function(p) pmin(pmax(p,1e-4),1-1e-4)
pool <- function(d){
  if(n_distinct(d$study_id)<3) return(c(pool=NA,lo=NA,hi=NA,k=n_distinct(d$study_id)))
  p<-clip01(d$CIR); se<-pmax(d$seCIR,1e-4); yi<-log(p/(1-p)); vi<-(se/(p*(1-p)))^2
  f<-tryCatch(rma(yi=yi,vi=vi,method='REML'),error=function(e)NULL)
  if(is.null(f)) return(c(pool=NA,lo=NA,hi=NA,k=n_distinct(d$study_id)))
  z<-qnorm(.975); c(pool=plogis(f$b),lo=plogis(f$b-z*sqrt(f$vb)),hi=plogis(f$b+z*sqrt(f$vb)),k=n_distinct(d$study_id))
}
cells <- per %>% filter(study_id==50) %>% distinct(scen_code,win,endp_code,ctdna,tp)
rows<-list();ii<-0
for(r in seq_len(nrow(cells))){
  cc<-cells[r,]
  d <- per %>% inner_join(cc,by=c('scen_code','win','endp_code','ctdna','tp'))
  full<-pool(d); drop<-pool(d %>% filter(study_id!=50))
  ii<-ii+1
  rows[[ii]]<-data.frame(cc, k_full=full['k'], pool_full=full['pool'], pool_lo=full['lo'], pool_hi=full['hi'],
                         k_drop=drop['k'], pool_drop=drop['pool'], drop_lo=drop['lo'], drop_hi=drop['hi'],
                         delta_pp=100*(full['pool']-drop['pool']))
}
res<-bind_rows(rows)
write_csv(res,'outputs/M1/M1b_leaveoneout_lygre.csv')
print(as.data.frame(res %>% mutate(across(where(is.numeric),~round(.x,3)))))
