# ============================================================================
# M2b_figures.R  longitudinal-trajectory figures (publication-grade, white background, English)
#   Fig 4A GALAXY (Kotani) non-metastatic DFS, four-trajectory absolute recurrence risk (the only setting with sufficient events)
#   Fig 4B cross-study 12-month trajectory CIR ranking (point shape = single-trajectory eligibility; hollow = events <=10, direction only)
# ============================================================================
suppressPackageStartupMessages({library(readr);library(dplyr);library(ggplot2);library(survival)})
set.seed(20261013)  # fix geom_jitter placement so Figure 4B is exactly reproducible
ipd <- read_csv("outputs/M2/M2_ipd_trajectory_reconstructed.csv",show_col_types=FALSE)
ar  <- read_csv("outputs/M2/M2_trajectory_absolute_risk.csv",show_col_types=FALSE)
out <- "outputs/M2"
save_pub <- function(plot, name, w, h) {
  dir.create("figures", showWarnings=FALSE, recursive=TRUE)
  ggsave(paste0("figures/", name, ".png"), plot, width=w, height=h, dpi=600, bg="white")
  ggsave(paste0("figures/", name, ".tif"), plot, width=w, height=h, dpi=600, bg="white",
         device="tiff", compression="lzw")
}
tcol <- c("sust_neg"="#2C3E50","clearance"="#5D9BD3","emergent"="#E67E22","persistent"="#C0392B","ever_pos"="#95A5A6")
tlab <- c(sust_neg="Sustained negative",clearance="Clearance (+ to -)",
          emergent="Emergent (+ after -)",persistent="Persistent positive",ever_pos="Ever-positive")
pub <- theme_bw(base_size=14)+theme(panel.background=element_rect(fill="white"),
  plot.background=element_rect(fill="white"),panel.grid.minor=element_blank(),
  legend.position="bottom",strip.text=element_text(size=12,face="bold"),
  axis.title=element_text(size=13),text=element_text(family="sans"))

# ---- Fig 4A: GALAXY four trajectories ----
g <- ipd %>% filter(kmid=="3-KM4") %>%
  mutate(trajectory=factor(trajectory,levels=names(tcol)))
sf <- survfit(Surv(time,status)~trajectory,data=g)
ss <- summary(sf)
step <- data.frame(time=ss$time,surv=ss$surv,lo=ss$lower,hi=ss$upper,strata=ss$strata)
step$trajectory <- sub("trajectory=","",step$strata)
step <- step %>% group_by(trajectory) %>%
  group_modify(~data.frame(time=c(0,.x$time),surv=c(1,.x$surv),lo=c(1,.x$lo),hi=c(1,.x$hi))) %>%
  mutate(CIR=1-surv,CIRlo=1-hi,CIRhi=1-lo,trajectory=factor(trajectory,levels=names(tcol)))
pa <- ggplot(step,aes(x=time,y=CIR*100,color=trajectory,fill=trajectory))+
  geom_ribbon(aes(ymin=CIRlo*100,ymax=CIRhi*100),alpha=0.12,linewidth=0,show.legend=FALSE)+
  geom_step(linewidth=1.1)+
  scale_color_manual(values=tcol,labels=tlab,name=NULL)+
  scale_fill_manual(values=tcol,labels=tlab,name=NULL)+
  scale_x_continuous(limits=c(0,26),expand=c(0,0))+
  scale_y_continuous(limits=c(0,100),expand=c(0,0))+
  labs(x="Months since landmark",y="Cumulative incidence of recurrence, CIR (%)",
       title="A  Non-metastatic CRC | serial ctDNA trajectory (GALAXY)")+
  pub+guides(color=guide_legend(nrow=2))

# ---- Fig 4B: cross-study 12-month trajectory ranking ----
b <- ar %>% filter(time==12, trajectory!="ever_pos") %>%
  mutate(trajectory=factor(trajectory,levels=c("sust_neg","clearance","emergent","persistent")),
         qual=ifelse(events>=11,"n>=11 events","n<=10 (direction only)"),
         author=factor(author,levels=rev(unique(author))))
pb <- ggplot(b,aes(x=trajectory,y=CIR*100,color=scen_code,shape=qual))+
  geom_hline(yintercept=0,color="grey80")+
  geom_jitter(width=0.12,height=0,size=3.2,alpha=0.9)+
  scale_x_discrete(labels=c("Sustained\nnegative","Clearance\n(+ to -)","Emergent\n(+ after -)","Persistent\npositive"))+
  scale_color_manual(values=c(MRD="#2C3E50",LARC="#27AE60",CRLM="#C0392B"),name="Disease setting")+
  scale_shape_manual(values=c("n>=11 events"=16,"n<=10 (direction only)"=1),name="Precision")+
  labs(x="Serial ctDNA trajectory",y="12-month cumulative incidence of recurrence, CIR (%)",
       title="B  Across studies | 12-month CIR by trajectory")+
  pub+theme(legend.position="bottom",legend.box="vertical")+
  guides(color=guide_legend(nrow=1,order=1),shape=guide_legend(nrow=1,order=2))

if (requireNamespace("patchwork", quietly=TRUE)) {
  fig <- patchwork::wrap_plots(pa,pb,nrow=1,widths=c(1,1.1))
  save_pub(fig, "Figure_4_trajectories", 15, 6.8)
} else {
  save_pub(pa, "Figure_4A_GALAXY_trajectory", 7, 6)
  save_pub(pb, "Figure_4B_across_studies", 8, 6)
}
cat("done\n")
