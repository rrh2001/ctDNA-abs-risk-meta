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
import pandas as pd, os
out=open("m2_hr.txt","w",encoding="utf-8")
# --- EDIT: local input paths (source files not redistributed) ---------------
WORK_DIR = os.environ.get('CTDNA_WORK_DIR', r'<LOCAL_PATH>/work')
MATRIX_XLSX = os.environ.get('CTDNA_MATRIX', r'<LOCAL_PATH>/data_availability_matrix.xlsx')
B=WORK_DIR; MAT=MATRIX_XLSX
s03=pd.read_excel(MAT,sheet_name="03_动态轨迹层(M2)")
# 治愈性场景（2 LARC/3 非转移/4 CRLM），仅HR类，排除5a/5b/6
cure=s03[s03["场景"].astype(str).str[:1].isin(["2","3","4"])].copy()
cure=cure[cure["可用性等级"].astype(str).str.contains("C_仅HR",na=False)]
cure["HR"]=pd.to_numeric(cure["HR"],errors="coerce")
cure=cure.dropna(subset=["HR"])
cols=["场景","StudyID","第一作者","年份","终点","轨迹","比较组名称","HR"]
cure=cure[cure["轨迹"].isin(["持续阳性 persistent","新发/转阳 emergent","清除/下降 clearance","参照-持续阴性"])]
cure.to_csv(os.path.join(B,"M2_out","M2_traj_HR_narrative.csv"),index=False,encoding="utf-8-sig")
out.write("=== 治愈性场景 仅HR 轨迹效应（叙述性，不合并）n=%d ===\n"%len(cure))
out.write(cure[cols].to_string(index=False)+"\n")
out.write("\n按轨迹HR范围:\n")
for t,g in cure.groupby("轨迹"):
    out.write("  %-28s n=%d HR %.2f-%.2f 中位%.2f\n"%(t,len(g),g.HR.min(),g.HR.max(),g.HR.median()))
out.close();print("done")
