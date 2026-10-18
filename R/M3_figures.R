# ============================================================================
# M3_figures.R  surveillance window/interval decision figures (publication-grade, white background, English)
#   Fig 5A non-metastatic post-op MRD +/- absolute recurrence-risk curves with the first-24-month concentration window (96% of recurrences)
#   Fig 5B stratified surveillance calendar (positive: every 3 months for 24 months; negative: every 6 months) with an illustrative 8.7-month lead time
# ============================================================================
suppressPackageStartupMessages({library(readr);library(dplyr);library(ggplot2);library(patchwork)})
save_pub <- function(plot, name, w, h) {
  dir.create("figures", showWarnings=FALSE, recursive=TRUE)
  ggsave(paste0("figures/", name, ".png"), plot, width=w, height=h, dpi=600, bg="white")
  ggsave(paste0("figures/", name, ".tif"), plot, width=w, height=h, dpi=600, bg="white",
         device="tiff", compression="lzw")
}
cc <- read_csv("outputs/M1/M1c_continuous_curves.csv",show_col_types=FALSE)
pub <- theme_bw(base_size=14)+theme(panel.background=element_rect(fill="white"),
  plot.background=element_rect(fill="white"),panel.grid.minor=element_blank(),
  legend.position="bottom",text=element_text(family="sans"),
  axis.title=element_text(size=13))

# ---- Panel A ----
g <- cc %>% filter(scen_code=="MRD",win=="postopMRD",endp_code=="DFS") %>%
  mutate(grp=ifelse(ctdna==1,"ctDNA-positive","ctDNA-negative"))
pa <- ggplot(g,aes(x=time_month,y=CIR*100,color=grp,fill=grp))+
  annotate("rect",xmin=0,xmax=24,ymin=-Inf,ymax=Inf,fill="grey90",alpha=0.5)+
  geom_ribbon(aes(ymin=CIRlo*100,ymax=CIRhi*100),alpha=0.12,linewidth=0,show.legend=FALSE)+
  geom_line(linewidth=1.2)+
  annotate("text",x=12,y=86,label="~96% of recurrences accrue\nin first 24 months (positive)",
           size=3.6,hjust=0.5,color="#C0392B")+
  annotate("text",x=32,y=33,label="negative:\nsteady, low",size=3.4,hjust=0.5,color="#2C3E50")+
  scale_color_manual(values=c("ctDNA-positive"="#C0392B","ctDNA-negative"="#2C3E50"),name=NULL)+
  scale_fill_manual(values=c("ctDNA-positive"="#C0392B","ctDNA-negative"="#2C3E50"),name=NULL)+
  scale_x_continuous(limits=c(0,60),expand=c(0,1))+
  scale_y_continuous(limits=c(0,100),expand=c(0,0))+
  labs(x="Months since post-op MRD landmark",y="Cumulative incidence of recurrence, CIR (%)",
       title="A  Recurrence is front-loaded in the first 24 months")+pub

# ---- Panel B: surveillance calendar ----
pos_t <- sort(unique(c(seq(0,24,3),seq(30,36,6))))
neg_t <- seq(0,36,6)
cal <- bind_rows(
  data.frame(t=pos_t,y=1,grp="ctDNA-positive"),
  data.frame(t=neg_t,y=0,grp="ctDNA-negative"))
pb <- ggplot(cal,aes(x=t,y=y))+
  annotate("rect",xmin=0,xmax=24,ymin=0.55,ymax=1.45,fill="#F4D6D2",alpha=0.35)+
  geom_hline(yintercept=c(0,1),color="grey70")+
  geom_hline(yintercept=0.5,color="grey90",linetype="dashed")+
  geom_point(aes(color=grp),size=3.4)+
  scale_color_manual(values=c("ctDNA-positive"="#C0392B","ctDNA-negative"="#2C3E50"),name=NULL)+
  annotate("text",x=0.2,y=1.18,hjust=0,size=3.5,
           label="MRD-positive: ctDNA + imaging every 3 mo (0-24), then every 6 mo")+
  annotate("text",x=0.2,y=-0.15,hjust=0,size=3.5,
           label="MRD-negative: surveillance every 6 mo")+
  # lead time arrow on positive arm
  annotate("segment",x=12,xend=20.7,y=1.32,yend=1.32,arrow=arrow(length=unit(0.12,"inches"),ends="both"),
           color="#E67E22",linewidth=0.9)+
  annotate("text",x=16.35,y=1.40,label="lead time ~8.7 mo (ctDNA precedes imaging)",
           size=3.2,color="#E67E22",hjust=0.5)+
  annotate("point",x=12,y=1.32,shape=17,color="#C0392B",size=2.6)+
  annotate("text",x=11.2,y=1.24,label="ctDNA+",size=2.8,hjust=1,color="#C0392B")+
  annotate("point",x=20.7,y=1.32,shape=15,color="black",size=2.4)+
  annotate("text",x=21.5,y=1.24,label="imaging",size=2.8,hjust=0,color="black")+
  annotate("text",x=12,y=0.64,label="intensive window",size=3.2,color="#C0392B",alpha=0.8)+
  annotate("text",x=18,y=-0.33,hjust=0.5,size=3.2,fontface="italic",
           label="MRD test at 4-12 weeks post-op sets the pathway")+
  scale_x_continuous(limits=c(0,36),breaks=seq(0,36,3))+
  scale_y_continuous(limits=c(-0.42,1.55),breaks=NULL)+
  labs(x="Months after surgery",y=NULL,title="B  Risk-stratified surveillance schedule")+
  pub+theme(panel.grid.major.y=element_blank(),axis.ticks.y=element_blank())

fig <- wrap_plots(pa,pb,nrow=2,heights=c(1.15,1))
save_pub(fig, "Figure_5_surveillance_framework", 11, 11)
cat("done\n")
