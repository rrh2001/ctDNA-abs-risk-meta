# -*- coding: utf-8 -*-
"""Fig6 hi-res: internal LOO calibration, double-column 7.2in @600dpi.
Same drawing as make_figure6_calibration.py; only figsize enlarged to a
double-column width so tight-bbox output is >=4320 px wide at 600 dpi."""
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import os

ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
M6=os.path.join(ROOT,"outputs","M6")
FIG=os.path.join(ROOT,"figures"); os.makedirs(FIG,exist_ok=True)
pr=pd.read_csv(os.path.join(M6,"M6_internal_calibration_pairs.csv"),encoding="utf-8-sig")
ex=pd.read_csv(os.path.join(M6,"S7_external_calibration_summary.csv"),encoding="utf-8-sig")
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10,"axes.edgecolor":"black",
                     "axes.linewidth":1.0,"figure.dpi":120})

fig,ax=plt.subplots(figsize=(7.9,7.5))
neg=pr[pr.ctdna==0]; pos=pr[pr.ctdna==1]
ax.scatter(neg.pred_loo,neg.obs,s=np.sqrt(neg.n)*7,facecolor="#3b6fb6",edgecolor="white",
           linewidth=0.5,alpha=0.65,zorder=3,label="ctDNA-negative (LOO)")
ax.scatter(pos.pred_loo,pos.obs,s=np.sqrt(pos.n)*7,facecolor="#d1495b",edgecolor="white",
           linewidth=0.5,alpha=0.7,zorder=3,label="ctDNA-positive (LOO)")

def trend(d,color):
    d=d.sort_values("pred_loo")
    x=d.pred_loo.values; y=d.obs.values; w=d.n.values.astype(float)
    grid=np.linspace(x.min(),x.max(),60); out=[]
    bw=max(0.16,(x.max()-x.min())*0.22)
    for g in grid:
        k=np.exp(-0.5*((x-g)/bw)**2)
        out.append(np.sum(k*w*y)/np.sum(k*w))
    ax.plot(grid,out,color=color,lw=2.6,zorder=4)
trend(neg,"#1f4e8c"); trend(pos,"#9e2232")

li=ex[ex.dataset.str.contains("Li")]
for _,r in li.iterrows():
    c="#3b6fb6" if r.ctdna==0 else "#d1495b"
    ax.scatter(r.predicted,r.observed,marker="*",s=520,facecolor=c,edgecolor="black",
               linewidth=1.3,zorder=6)
ax.scatter([],[],marker="*",s=300,facecolor="white",edgecolor="black",label="Li2022 external")

ax.plot([0,1],[0,1],color="black",lw=1.3,ls="--",zorder=2)
ax.set_xlim(0,1); ax.set_ylim(0,1)
ax.set_xlabel("Predicted 12/24/36-month cumulative incidence of recurrence\n(leave-one-study-out pooled)")
ax.set_ylabel("Observed reconstructed cumulative incidence")
ax.set_xticks(np.arange(0,1.01,0.2)); ax.set_yticks(np.arange(0,1.01,0.2))
ax.set_aspect("equal")
ax.text(0.03,0.965,"calibration slope 0.91 (95% CI 0.79\u20131.11)\nintercept \u22120.13 (95% CI \u22120.57 to 0.28)",
        transform=ax.transAxes,va="top",ha="left",fontsize=10,
        bbox=dict(boxstyle="round,pad=0.35",fc="white",ec="0.7",lw=0.7))
ax.legend(loc="lower right",frameon=True,framealpha=0.95,edgecolor="0.7",fontsize=9.5,
          handles=[
            Line2D([],[],marker="o",color="w",markerfacecolor="#3b6fb6",markersize=9,label="ctDNA-negative (LOO)"),
            Line2D([],[],marker="o",color="w",markerfacecolor="#d1495b",markersize=9,label="ctDNA-positive (LOO)"),
            Line2D([],[],marker="*",color="w",markerfacecolor="white",markeredgecolor="black",markersize=15,label="Li2022 external")])
plt.tight_layout()
for ext in ["png","tif"]:
    kw = {} if ext == "png" else {"pil_kwargs": {"compression": "tiff_lzw"}}
    plt.savefig(os.path.join(FIG,f"Figure_6_internal_calibration.{ext}"),dpi=600,
                bbox_inches="tight",facecolor="white",**kw)
print("saved hi-res Figure 6")
