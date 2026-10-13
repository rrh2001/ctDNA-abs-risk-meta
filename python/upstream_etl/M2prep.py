# -*- coding: utf-8 -*-
# =============================================================================
# UPSTREAM, ONE-TIME EXTRACTION / TRANSFORMATION SCRIPT.
# NOT part of the sealed reproduction chain, which starts from data/ and R/.
#
# It reads the authors' ORIGINAL master extraction workbook and the digitized
# curve folder, which use Chinese worksheet names, column headers and group
# labels and are NOT redistributed. The Chinese string literals below are
# therefore deliberate matchers against those source files, not user-facing
# text. To re-run, set the LOCAL input paths to your own copy (available from
# the corresponding author on reasonable request, subject to publisher terms).
# =============================================================================
import pandas as pd, numpy as np, os, glob
out=open("M2prep_result.txt","w",encoding="utf-8")
# --- EDIT: local input paths (source files not redistributed) ---------------
WORK_DIR = os.environ.get('CTDNA_WORK_DIR', r'<LOCAL_PATH>/work')
DIGITIZED_DIR = os.environ.get('CTDNA_DIGITIZED', r'<LOCAL_PATH>/digitized')
MASTER_XLSX = os.environ.get('CTDNA_MASTER_XLSX', r'<LOCAL_PATH>/master_extraction.xlsx')  # original '全表.xlsx'
B=WORK_DIR; DIG=DIGITIZED_DIR; QB=MASTER_XLSX
os.makedirs(os.path.join(B,"M2_input_bundle"),exist_ok=True)

# 11个可重建多态KM（kmid,digfile,study,author,year,scen,win,endpoint）
KMS=[
 ("3-KM4","3-KM4",3,"Kotani D",2023,"MRD","dynamic","DFS"),
 ("21-KM2c","21-KM2c",21,"Rubio-Alarcon C",2026,"MRD","postACT","TTR"),
 ("18-KM2a","18-KM2a",18,"Raunkilde L",2025,"CRLM","dynamic","PFS"),
 ("18-KM2b","18-KM2b",18,"Raunkilde L",2025,"CRLM","dynamic","OS"),
 ("27-KM2A","27-KM2A",27,"Newhook TE",2023,"CRLM","dynamic","RFS"),
 ("30-KM4c","30-KM4c",30,"Zhou J",2025,"LARC","dynamic","PFS"),
 ("30-KM4c-OS","30-KM4c-OS",30,"Zhou J",2025,"LARC","dynamic","OS"),
 ("30-KM4d","30-KM4d",30,"Zhou J",2025,"LARC","preop","PFS"),
 ("30-KM4d-OS","30-KM4d-OS",30,"Zhou J",2025,"LARC","preop","OS"),
 ("13-KM3A","13-KM3A",13,"Li Y",2024,"CRLM","postACT","RFS"),
 ("4-KM3C","4-KM3C",4,"Mo S",2023,"MRD","dynamic","RFS"),
 ("4-KM3D","4-KM3D",4,"Mo S",2023,"MRD","onTx","RFS"),
 ("5-KM-S6A","5-KM-S6A",5,"Wang Y",2021,"LARC","dynamic","RFS"),
]
# curve_raw -> 标准轨迹
TRAJ={
 "3-KM4":{1:"sust_neg",2:"emergent",3:"clearance",4:"persistent"},
 "21-KM2c":{1:"sust_neg",2:"emergent",3:"clearance",4:"persistent"},
 "18-KM2a":{1:"persistent",2:"sust_neg",3:"clearance",4:"emergent"},
 "18-KM2b":{1:"persistent",2:"sust_neg",3:"clearance",4:"emergent"},
 "27-KM2A":{1:"persistent",2:"clearance",3:"sust_neg"},
 "30-KM4c":{1:"persistent",2:"clearance",3:"sust_neg"},
 "30-KM4c-OS":{1:"persistent",2:"clearance",3:"sust_neg"},
 "30-KM4d":{1:"persistent",2:"clearance",3:"sust_neg"},
 "30-KM4d-OS":{1:"persistent",2:"clearance",3:"sust_neg"},
 "13-KM3A":{1:"clearance",2:"persistent"},
 "4-KM3C":{1:"sust_neg",2:"ever_pos"},
 "4-KM3D":{1:"sust_neg",2:"persistent"},
 "5-KM-S6A":{1:"clearance",2:"persistent"},
}
narlong=pd.read_excel(QB,sheet_name="KM时间点长表")
out.write("NAR long cols: "+repr(list(narlong.columns))+"\n")
# 找时间点与NAR列
tcol=[c for c in narlong.columns if "时间点_月" in str(c)][0]
ncol=[c for c in narlong.columns if "Number" in str(c) or "风险" in str(c) or "at_risk" in str(c).lower()][0]
ccol=[c for c in narlong.columns if "digitized" in str(c).lower()][0]
out.write("use tcol=%s ncol=%s ccol=%s\n"%(tcol,ncol,ccol))

curves=[]; nar=[]; meta=[]
for kmid,digf,sid,au,yr,scen,win,ep in KMS:
    d=pd.read_csv(os.path.join(DIG,digf+"_digitized.csv"))
    d=d.rename(columns={"digitized_curve":"curve"})
    d=d[["KMID","curve","time","St"]].copy()
    d["kmid"]=kmid
    curves.append(d)
    # NAR：主KMID或派生前缀
    nn=narlong[narlong["KMID"].astype(str).str.startswith(kmid)]
    # 若NAR登记在派生KMID（如3-KM4-CN），也按主曲线curves；统一只保留有curve编号且与digfile同轨迹
    if len(nn):
        for col in [ccol]:
            tmp=nn[["KMID",col,tcol,ncol]].copy()
            tmp.columns=["nar_kmid","curve","time","nrisk"]
            tmp=tmp[pd.to_numeric(tmp["curve"],errors="coerce").notna()]
            tmp["curve"]=tmp["curve"].astype(int)
            tmp["time"]=pd.to_numeric(tmp["time"],errors="coerce")
            tmp["nrisk"]=pd.to_numeric(tmp["nrisk"],errors="coerce")
            tmp=tmp.dropna(subset=["time","nrisk"])
            # 多派生KMID会重复curve，按curve+time去重（取最大nrisk，因派生CN/CP/PP可能各是子图）
            tmp["kmid"]=kmid
            nar.append(tmp)
    traj=TRAJ[kmid]
    for c in sorted(d["curve"].unique()):
        meta.append(dict(kmid=kmid,study_id=sid,author=au,year=yr,scen_code=scen,
            win=win,endp_code=ep,curve=int(c),trajectory=traj.get(int(c),"?")))
curves=pd.concat(curves,ignore_index=True)
meta=pd.DataFrame(meta)
if nar:
    nar=pd.concat(nar,ignore_index=True)
    nar=nar.sort_values(["kmid","curve","time"]).drop_duplicates(["kmid","curve","time"],keep="first")
else:
    nar=pd.DataFrame(columns=["kmid","curve","time","nrisk"])

# NAR覆盖核查
summ=[]
for (kmid,cv),g in curves.groupby(["kmid","curve"]):
    nsub=nar[(nar.kmid==kmid)&(nar.curve==cv)]
    summ.append(dict(kmid=kmid,curve=cv,npts=len(g),
        nar_pts=nsub.nrisk.notna().sum(), nar_first=nsub.nrisk.iloc[0] if len(nsub) else np.nan,
        nar_times=";".join(str(round(x,1)) for x in nsub.time.head(12))))
summ=pd.DataFrame(summ).merge(meta,on=["kmid","curve"])
summ.to_csv(os.path.join(B,"M2_input_bundle","M2_traj_coverage.csv"),index=False,encoding="utf-8-sig")
curves.to_csv(os.path.join(B,"M2_input_bundle","curves2.csv"),index=False,encoding="utf-8-sig")
nar.to_csv(os.path.join(B,"M2_input_bundle","nar2.csv"),index=False,encoding="utf-8-sig")
meta.to_csv(os.path.join(B,"M2_input_bundle","traj_meta.csv"),index=False,encoding="utf-8-sig")
out.write("\n=== NAR coverage per trajectory ===\n")
out.write(summ[["kmid","author","scen_code","endp_code","curve","trajectory","npts","nar_pts","nar_first"]].to_string(index=False)+"\n")
out.write("\ncurves rows %d, nar rows %d, meta rows %d\n"%(len(curves),len(nar),len(meta)))
out.close();print("done")
