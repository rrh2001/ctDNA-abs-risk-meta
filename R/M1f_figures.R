# ============================================================================
# M1f_figures.R  publication figures (English, white background, minimal text; high resolution, submission-ready)
#   Fig 2 continuous absolute-risk curves (key windows of the three curative scenarios, ctDNA+/-)
#   Fig 3 forest plot of absolute risk difference (RD) (DFS, k>=3)
# ============================================================================
suppressPackageStartupMessages({library(readr);library(dplyr);library(ggplot2)})
OUT <- "outputs/M1"
# Save a submission-ready figure to figures/ as 600 dpi PNG and LZW TIFF.
save_pub <- function(plot, name, w, h) {
  dir.create("figures", showWarnings=FALSE, recursive=TRUE)
  ggsave(paste0("figures/", name, ".png"), plot, width=w, height=h, dpi=600, bg="white")
  ggsave(paste0("figures/", name, ".tif"), plot, width=w, height=h, dpi=600, bg="white",
         device="tiff", compression="lzw")
}
RED <- "#C0392B"; NAVY <- "#2C3E50"
pub <- theme_bw(base_size=13) +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_line(color="grey90",linewidth=.4),
        panel.background=element_rect(fill="white"),
        plot.background=element_rect(fill="white"),
        panel.border=element_rect(color="black",linewidth=.6),
        strip.background=element_rect(fill="white",color="black"),
        strip.text=element_text(face="bold",hjust=0,size=10.5),
        legend.position="bottom", legend.key=element_rect(fill="white"),
        axis.title=element_text(size=12), plot.title=element_blank())

# ---------------- Fig 2 continuous curves ----------------
curv <- read_csv(file.path(OUT,"M1c_continuous_curves.csv"), show_col_types=FALSE)
pan <- data.frame(
  scen_code=c("MRD","MRD","LARC","CRLM","CRLM","CRLM"),
  win=c("postopMRD","postACT","postopMRD","postopMRD","postACT","preop"),
  panel=c("A  Non-metastatic CRC | post-op MRD",
          "B  Non-metastatic CRC | post-ACT",
          "C  LARC | post-op MRD",
          "D  Resectable CRLM | post-op MRD",
          "E  Resectable CRLM | post-ACT",
          "F  Resectable CRLM | pre-op"),
  ord=1:6)
d <- curv %>% inner_join(pan, by=c("scen_code","win")) %>%
  mutate(Status=ifelse(ctdna==1,"ctDNA-positive","ctDNA-negative"),
         panel=factor(panel, levels=pan$panel[order(pan$ord)]))
fig2 <- ggplot(d, aes(time_month, CIR*100, color=Status, fill=Status)) +
  geom_ribbon(aes(ymin=CIRlo*100, ymax=CIRhi*100), alpha=.13, linetype=0) +
  geom_line(linewidth=1.0) +
  facet_wrap(~panel, nrow=2, scales="free_y") +
  scale_color_manual(values=c("ctDNA-positive"=RED,"ctDNA-negative"=NAVY)) +
  scale_fill_manual(values=c("ctDNA-positive"=RED,"ctDNA-negative"=NAVY)) +
  labs(x="Months since landmark / surgery", y="Cumulative incidence of recurrence, CIR (%)",
       color=NULL, fill=NULL) +
  pub + guides(color=guide_legend(nrow=1), fill=guide_legend(nrow=1))
save_pub(fig2, "Figure_2_absolute_risk_curves", 10.5, 6.8)

# ---------------- Fig 3 RD forest plot ----------------
con <- read_csv(file.path(OUT,"M1_contrast_pooled.csv"), show_col_types=FALSE)
poo <- read_csv(file.path(OUT,"M1_absrisk_pooled.csv"), show_col_types=FALSE)
kpos <- poo %>% filter(ctdna==1, endp_code=="DFS") %>%
  transmute(scen_code,win,endp_code,tp,i2)
f <- con %>% filter(endp_code=="DFS", !insufficient) %>%
  left_join(kpos, by=c("scen_code","win","endp_code","tp")) %>%
  mutate(scen_lab=recode(scen_code,LARC="LARC",MRD="Non-metastatic CRC",CRLM="Resectable CRLM"),
         win_lab=recode(win,preop="pre-op",postopMRD="post-op MRD",postNAC="post-NAC",postACT="post-ACT"),
         label=paste0(scen_lab," | ",win_lab," | ",tp," mo"),
         lab=paste0("k=",k,", I\u00b2=",round(i2,0),"%"),
         rd=RD*100, lo=RDlo*100, hi=RDhi*100)
ord <- f %>% arrange(scen_code,win,tp) %>% pull(label) %>% unique()
f <- f %>% mutate(label=factor(label, levels=rev(ord)))
fig3 <- ggplot(f, aes(rd,label)) +
  geom_vline(xintercept=0, color="grey40", linetype="dashed", linewidth=.5) +
  geom_errorbarh(aes(xmin=lo,xmax=hi), height=.18, color=NAVY, linewidth=.7) +
  geom_point(size=3.2, color=RED) +
  geom_text(aes(label=lab, x=108), hjust=1, size=3.1, color="grey15") +
  scale_x_continuous(limits=c(-15,115), breaks=seq(-10,80,10)) +
  labs(x="Absolute risk difference (ctDNA-positive minus negative), percentage points",
       y=NULL, title=NULL) +
  pub + theme(panel.grid.major.y=element_blank(), axis.text.y=element_text(size=10))
save_pub(fig3, "Figure_3_RD_forest_plot", 10.5, 6.2)
cat("figures written:", length(unique(d$panel)), "Fig2 panels;", nrow(f), "Fig3 rows\n")
