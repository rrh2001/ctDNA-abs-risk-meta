# -*- coding: utf-8 -*-
"""
M6: (1) quantitative competing-risk sensitivity (parametric Aalen-Johansen scenarios)
    (2) internal calibration: leave-one-study-out + cluster bootstrap slope/intercept
    (3) external calibration summary (Li/Ogaard): O-E, E/O, calibration-in-the-large
All inputs are EXISTING M1 outputs; no new data collected.
"""
import numpy as np, pandas as pd, os
np.random.seed(20261014)
# Run from anywhere; paths are relative to the repository root.
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
M1 = os.path.join(ROOT, "outputs", "M1")
OUT = os.path.join(ROOT, "outputs", "M6"); os.makedirs(OUT, exist_ok=True)

# ---------- inputs ----------
ipd = pd.read_csv(os.path.join(M1, "M1_ipd_reconstructed.csv"), encoding="utf-8-sig")
ps  = pd.read_csv(os.path.join(M1, "M1_absrisk_perstudy.csv"), encoding="utf-8-sig")
pool= pd.read_csv(os.path.join(M1, "M1_absrisk_pooled.csv"), encoding="utf-8-sig")
# Released S5 is written by M1d when the (non-redistributed) appendix IPD is
# present; otherwise the released S5 table is used.
ext = pd.read_csv(os.path.join(M1, "S5_external_calibration.csv"), encoding="utf-8-sig")

MAIN = [("CRLM","postACT","DFS"),("CRLM","postopMRD","DFS"),("CRLM","preop","DFS"),
        ("LARC","postopMRD","DFS"),("MRD","postACT","DFS"),("MRD","postopMRD","DFS")]
def in_main(r): return (r.scen_code,r.win,r.endp_code) in MAIN
psm = ps[ps.apply(in_main,axis=1)].copy()

# kmid -> cell map (scen_code, win, endp)
kmap = ps[["kmid","scen_code","win","endp_code"]].drop_duplicates("kmid")
ipd = ipd.merge(kmap, on="kmid", how="left", suffixes=("","_m"))
# prefer mapped endp; ipd.endpoint is Chinese/EN; use mapped endp_code where present
ipd["endp"] = ipd["endp_code"]
ipd = ipd[ipd.apply(lambda r:(r.scen_code,r.win,r.endp) in MAIN, axis=1)].copy()

def logit(p): return np.log(p/(1-p))
def invlogit(x): return 1/(1+np.exp(-x))

# ============================================================
# PART 1  competing-risk sensitivity
# ============================================================
def km_steps(t, e):
    """return (event_times, S_before, S_after) KM"""
    t=np.asarray(t,float); e=np.asarray(e,int)
    od=np.argsort(t); t=t[od]; e=e[od]
    ut=np.unique(t)
    S=1.0; rows=[]
    for u in ut:
        idx=t==u
        n=np.sum(t>=u)               # at risk
        d=np.sum(e[idx]==1)
        if d>0:
            Sb=S; Sa=S*(1-d/n)
            rows.append((u,Sb,Sa))
            S=Sa
    return pd.DataFrame(rows,columns=["time","Sb","Sa"])

def cir_km(step,T):
    return 1-step[step.time<=T].Sa.min() if (step.time<=T).any() else 0.0

def cir_aj(step,T,lamc_m):
    """Aalen-Johansen CIR with constant competing hazard (per month):
       CIR_AJ(T)=sum_{u<=T} exp(-lamc*u)*(S(u-)-S(u))"""
    s=step[step.time<=T]
    if len(s)==0: return 0.0
    return float(np.sum(np.exp(-lamc_m*s.time)*(s.Sb-s.Sa)))

rows=[]
LAM_YEAR=[0.005,0.01,0.02,0.03]
for (sc,win,endp) in MAIN:
    d=ipd[(ipd.scen_code==sc)&(ipd.win==win)&(ipd.endp==endp)]
    for g in [0,1]:
        dd=d[d.ctdna==g]
        if len(dd)==0: continue
        st=km_steps(dd.time.values, dd.status.values)
        for T in [12,24,36]:
            if T> dd.time.max(): continue
            ckm=cir_km(st,T)
            r={"scen_code":sc,"win":win,"endp":endp,"ctdna":g,"tp":T,"CIR_KM":ckm}
            for ly in LAM_YEAR:
                r[f"CIR_AJ_{ly}"]=cir_aj(st,T,ly/12.0)
            r["diff_pp_1pct"]=(r["CIR_AJ_0.01"]-ckm)*100
            r["diff_pp_3pct"]=(r["CIR_AJ_0.03"]-ckm)*100
            rows.append(r)
cr=pd.DataFrame(rows)
cr.to_csv(f"{OUT}/M6_competing_risk_sensitivity.csv",index=False,encoding="utf-8-sig")

# ============================================================
# PART 2  LOO internal calibration
# ============================================================
def re_pool(theta,var):
    """DerSimonian-Laird random-effects pooling on logit scale. returns mean, tau2"""
    w=1/var
    m=np.sum(w*theta)/np.sum(w)
    Q=np.sum(w*(theta-m)**2); C=np.sum(w)-np.sum(w**2)/np.sum(w)
    tau2=max(0,(Q-(len(theta)-1))/C) if C>0 else 0
    wr=1/(var+tau2); mm=np.sum(wr*theta)/np.sum(wr)
    return mm,tau2

def build_pairs(df):
    out=[]
    for keys,g in df.groupby(["scen_code","win","endp_code","ctdna","tp"]):
        sc,win,endp,gg,tp=keys
        if len(g)<2: continue
        y=g.CIR.values; se=g.seCIR.values.astype(float); nn=g.n.values.astype(float)
        # continuity correction for zero-event / boundary arms (se reported as 0)
        se=np.where(se<=0, 0.5/nn, se)
        y=np.clip(y,1e-4,1-1e-4)
        th=logit(y); va=(se**2)/((y*(1-y))**2)
        va=np.clip(va,1e-6,None)
        full,_=re_pool(th,va)
        for i in range(len(g)):
            keep=np.ones(len(g),bool); keep[i]=False
            if keep.sum()<1: continue
            mloo,tau=re_pool(th[keep],va[keep])
            out.append(dict(kmid=g.kmid.iloc[i],study_id=g.study_id.iloc[i],
                scen_code=sc,win=win,endp=endp,ctdna=gg,tp=tp,
                n=g.n.iloc[i],events=g.events.iloc[i],
                obs=g.CIR.iloc[i],
                pred_loo=invlogit(mloo),pred_app=invlogit(full)))
    return pd.DataFrame(out)

pairs=build_pairs(psm)
pairs.to_csv(f"{OUT}/M6_internal_calibration_pairs.csv",index=False,encoding="utf-8-sig")

def irls_calib(obs,n,pred):
    """binomial logistic calibration: events~Bin(n,p), logit(p)=a+b*logit(pred)"""
    ev=np.clip(np.round(n*obs),0,n).astype(int)
    q=logit(np.clip(pred,1e-4,1-1e-4))
    X=np.column_stack([np.ones(len(q)),q])
    beta=np.zeros(2)
    for _ in range(60):
        p=invlogit(X@beta); p=np.clip(p,1e-9,1-1e-9)
        W=n*p*(1-p); z=X@beta+(ev-n*p)/(n*p*(1-p))
        XW=X*W[:,None]
        try: nb=np.linalg.solve(XW.T@X, XW.T@z)
        except np.linalg.LinAlgError: return np.array([np.nan,np.nan])
        if np.max(np.abs(nb-beta))<1e-9: beta=nb;break
        beta=nb
    return beta

def calib_metrics(d, predcol, label):
    res=[]
    for name,sub in [("overall",d),("negative",d[d.ctdna==0]),("positive",d[d.ctdna==1])]:
        if len(sub)<8: continue
        b=irls_calib(sub.obs.values,sub.n.values.astype(float),sub[predcol].values)
        res.append(dict(method=label,stratum=name,n_pts=len(sub),
                        intercept=b[0],slope=b[1]))
    return res

met=[]
met+=calib_metrics(pairs,"pred_app","apparent(in-sample)")
met+=calib_metrics(pairs,"pred_loo","LOO cross-validated")

# cluster bootstrap by study_id (study-level LOO pairs, stratum-specific)
B=1000
ci=lambda v:(np.nanpercentile(v,2.5),np.nanpercentile(v,97.5))
for name,mask in [("overall",np.ones(len(pairs),bool)),
                  ("negative",pairs.ctdna.values==0),
                  ("positive",pairs.ctdna.values==1)]:
    d=pairs[mask]
    if len(d)<8: continue
    st=d.study_id.unique()
    bb=[]
    for b in range(B):
        bss=np.random.choice(st,len(st),replace=True)
        db=pd.concat([d[d.study_id==s] for s in bss],ignore_index=True)
        bb.append(irls_calib(db.obs.values,db.n.values.astype(float),db.pred_loo.values))
    bb=np.array(bb)
    bp=irls_calib(d.obs.values,d.n.values.astype(float),d.pred_loo.values)
    ci_a=ci(bb[:,0]); ci_b=ci(bb[:,1])
    met.append(dict(method="LOO bootstrap",stratum=name,n_pts=len(d),
        intercept=bp[0],slope=bp[1],
        intercept_lo=ci_a[0],intercept_hi=ci_a[1],slope_lo=ci_b[0],slope_hi=ci_b[1]))

metr=pd.DataFrame(met)
metr.to_csv(f"{OUT}/M6_internal_calibration_metrics.csv",index=False,encoding="utf-8-sig")

# decile calibration curve points (LOO), overall + by ctdna
dec=[]
for name,sub in [("overall",pairs),("negative",pairs[pairs.ctdna==0]),("positive",pairs[pairs.ctdna==1])]:
    if len(sub)<10: continue
    q=pd.qcut(sub.pred_loo.rank(method="first"),10,labels=False,duplicates="drop")
    for gg,g in sub.groupby(q):
        ev=np.sum(g.n*g.obs); nt=np.sum(g.n)
        dec.append(dict(stratum=name,bin=int(gg),n=nt,
                        pred=np.average(g.pred_loo,weights=g.n),
                        obs=ev/nt))
pd.DataFrame(dec).to_csv(f"{OUT}/M6_calibration_deciles.csv",index=False,encoding="utf-8-sig")

# ============================================================
# PART 3  external calibration
# ============================================================
e=ext[ext.predicted.notna()].copy()
e["OE"]=e.observed-e.predicted
e["EoverO"]=e.predicted/e.observed
e["calib_in_large_logit"]=logit(np.clip(e.observed,1e-4,1-1e-4))-logit(np.clip(e.predicted,1e-4,1-1e-4))
e.to_csv(f"{OUT}/S7_external_calibration_summary.csv",index=False,encoding="utf-8-sig")

# ---------- printed summary ----------
def pp(x):
    if pd.isna(x):return "NA"
    return f"{x:.3f}"
with open(f"{OUT}/M6_summary.txt","w",encoding="utf-8") as f:
    f.write("===== PART1 COMPETING-RISK SENSITIVITY (one-minus-KM vs AJ scenarios) =====\n")
    f.write("Annual non-recurrence-death hazard scenarios: 0.5%,1%,2%,3% per year\n")
    show=cr[cr.tp.isin([12,24,36])].copy()
    for _,r in show.iterrows():
        f.write(f"{r.scen_code:4s} {r.win:9s} ctDNA{r.ctdna} {int(r.tp)}m: "
                f"KM={r.CIR_KM:.3f} | AJ0.5%={r['CIR_AJ_0.005']:.3f} AJ1%={r['CIR_AJ_0.01']:.3f} "
                f"AJ2%={r['CIR_AJ_0.02']:.3f} AJ3%={r['CIR_AJ_0.03']:.3f} | "
                f"d@1%={r.diff_pp_1pct:+.2f}pp d@3%={r.diff_pp_3pct:+.2f}pp\n")
    f.write(f"\nMax |delta| pp at 1%/yr: {cr.diff_pp_1pct.abs().max():.2f} ; at 3%/yr: {cr.diff_pp_3pct.abs().max():.2f}\n")
    f.write("\n===== PART2 INTERNAL CALIBRATION (LOO, cluster bootstrap B=1000) =====\n")
    f.write(f"calibration pairs: {len(pairs)} from {pairs.study_id.nunique()} studies, "
            f"{pairs.kmid.nunique()} comparisons\n")
    for _,r in metr.iterrows():
        ci_s=""
        if "slope_lo" in r and not pd.isna(r.get("slope_lo",np.nan)):
            ci_s=f" slope95CI[{r['slope_lo']:.2f},{r['slope_hi']:.2f}] int95CI[{r['intercept_lo']:.2f},{r['intercept_hi']:.2f}]"
        f.write(f"{r['method']:24s} {r['stratum']:9s} n={r['n_pts']:3d} "
                f"intercept={pp(r['intercept'])} slope={pp(r['slope'])}{ci_s}\n")
    f.write("\n===== PART3 EXTERNAL CALIBRATION =====\n")
    for _,r in e.iterrows():
        f.write(f"{r['dataset']} ctDNA{r.ctdna} {int(r.tp)}m: obs={r.observed:.3f} pred={r.predicted:.3f} "
                f"O-E={r.OE:+.3f} E/O={r.EoverO:.2f} calib-large(logit)={r.calib_in_large_logit:+.2f}\n")
print(open(f"{OUT}/M6_summary.txt",encoding="utf-8").read())
