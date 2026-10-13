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

import sys, io, os, glob, re
import pandas as pd, numpy as np
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
pd.set_option('display.width',340); pd.set_option('display.max_columns',70); pd.set_option('display.max_colwidth',44)
# --- EDIT: local input paths (source files not redistributed) ---------------
WORK_DIR = os.environ.get('CTDNA_WORK_DIR', r'<LOCAL_PATH>/work')
MASTER_XLSX = os.environ.get('CTDNA_MASTER_XLSX', r'<LOCAL_PATH>/master_extraction.xlsx')  # original '全表.xlsx'
DIGITIZED_DIR = os.environ.get('CTDNA_DIGITIZED', r'<LOCAL_PATH>/digitized')
MATRIX_XLSX = os.environ.get('CTDNA_MATRIX', r'<LOCAL_PATH>/data_availability_matrix.xlsx')
base=WORK_DIR; QB=MASTER_XLSX; DIG=DIGITIZED_DIR; xl=MATRIX_XLSX
OUT=os.path.join(base,"M1_input_bundle"); os.makedirs(OUT,exist_ok=True)

inc=pd.read_csv(os.path.join(base,"M1_A级纳入清单_130KM.csv")); inc['KMID']=inc.KMID.astype(str)
cur=[]
for f in glob.glob(os.path.join(DIG,"*.csv")):
    d=pd.read_csv(f); d.columns=[x.strip() for x in d.columns]
    cur.append(d[['KMID','digitized_curve','time','St']].rename(columns={'digitized_curve':'curve'}))
cur=pd.concat(cur,ignore_index=True); cur['KMID']=cur.KMID.astype(str); cur['curve']=cur.curve.astype(int)
cur['time']=pd.to_numeric(cur.time,errors='coerce'); cur['St']=pd.to_numeric(cur.St,errors='coerce')
cur=cur.dropna(subset=['time','St'])
digids=set(cur.KMID.unique())
nar=pd.read_excel(QB,sheet_name='KM时间点长表')[['KMID','digitized_curve','曲线组','时间点_月','Number_at_risk']]
nar.columns=['KMID','curve','grp','time','nrisk']; nar['KMID']=nar.KMID.astype(str); nar['curve']=nar.curve.astype(int)
nar['time']=pd.to_numeric(nar.time,errors='coerce'); nar['nrisk']=pd.to_numeric(nar.nrisk,errors='coerce')
nar=nar.dropna(subset=['time','nrisk'])
s7=pd.read_excel(xl,sheet_name='07_数字化资产曲线级'); s7['KMID']=s7.KMID.astype(str); s7['curve']=s7.curve.astype(int)
s2=pd.read_excel(xl,sheet_name='02_逐KM主清单(M1)'); s2['KMID']=s2.KMID.astype(str)
cmp=pd.read_excel(QB,sheet_name='KM比较层'); cmp['KMID']=cmp.KMID.astype(str)
cmp=cmp[['KMID','组1事件数','组2事件数','总事件数']]
def nz(s):
    s=str(s)
    if '未报告' in s or s.lower()=='nan': return np.nan
    mm=re.search(r'\d+',s.replace(',','')); return float(mm.group()) if mm else np.nan
cmp['evt']=cmp['总事件数'].map(nz)

def norm(x): return re.sub(r'[\s_()（）]+','',str(x)).lower()
def _cls(t):
    t=norm(t)
    if re.match(r'c[12]$',t):  # 由外层结合02定义处理
        return 'c'
    if re.search(r'not\s*detect|undetected|未检出|不可检测',t): return 0
    if re.search(r'detect|presence|阳性|positive|persistent|emergent|metasta',t): return 1
    if re.search(r'non[\s\-]?pcr|nonpcr|无清除|non[\s_]?clearance|nonclearance',t): return 1
    if re.search(r'pcr|clearance|清除',t): return 0
    if re.search(r'pni[\s\-]*[+＋]',t): return 1
    if re.search(r'pni[\s\-]*[-\u2013\u2212]',t): return 0
    if re.search(r'negative|absence|阴性|sustainedneg|持续阴|低危',t): return 0
    if re.search(r'positive|阳性|高危',t): return 1
    if re.search(r'≥|>0|high|increase|升高',t): return 1
    if re.search(r'≤|<0|low|decrease|降低',t): return 0
    return None
def classify(grp, grp_def='', s2def1='', s2def2=''):
    g=norm(grp)
    # c1/c2：回退02组定义（C1低危 C2高危）
    if re.match(r'c[12]',g):
        blob=norm(s2def1)+norm(s2def2)+norm(grp_def)
        if g.startswith('c1'): return 0 if ('c1低危' in blob or '低危' in blob) else (1 if 'c1高危' in blob else None)
        if g.startswith('c2'): return 1 if ('c2高危' in blob or '高危' in blob) else (0 if 'c2低危' in blob else None)
    v=_cls(g)  # 先只按曲线组名判
    if v is not None: return v
    return _cls(g+norm(grp_def))  # 组名泛化(如ctDNA/采样时点)才用组定义兜底

# 标签源：07曲线级（曲线组+组定义）
lab7={}
for _,r in s7.iterrows():
    lab7[(r.KMID,int(r.curve))]=(str(r['曲线组']), str(r.get('组定义','')))
grp_nar=nar.groupby('KMID').apply(lambda d:{int(c):str(g) for c,g in zip(d.curve,d.grp)},include_groups=False).to_dict()
s2idx=s2.set_index('KMID')
nard=nar.groupby(['KMID','curve']).size().groupby('KMID').min().rename('nar_minpts')

TREAT=['observation','观察','actvs','adjuvantchemo','chemotherapy']
rows=[];excl=[]
for _,r in inc.iterrows():
    k=r.KMID
    if k not in digids:
        excl.append(dict(KMID=k,场景=r['场景'],landmark=r['landmark'],终点=r['终点'],分流='M2动态/多态派生(无独立digitized)')); continue
    curves=sorted(cur[cur.KMID==k].curve.unique())
    d1=s2idx.loc[k,'组1定义'] if k in s2idx.index else ''; d2=s2idx.loc[k,'组2定义'] if k in s2idx.index else ''
    gnames={}; gdefs={}
    for cv in curves:
        tup=lab7.get((k,cv))
        if tup: gnames[cv],gdefs[cv]=tup
        else: gnames[cv]=grp_nar.get(k,{}).get(cv,''); gdefs[cv]=''
    gn={cv:norm(gnames[cv]) for cv in curves}
    # 治疗/干预比较（按曲线组名，不受组定义里“ctDNA阳性人群”限定语影响）
    if any(re.search(r'rlts|重复局部治疗',gn[cv]) for cv in curves):
        excl.append(dict(KMID=k,场景=r['场景'],landmark=r['landmark'],终点=r['终点'],分流='非ctDNA干预比较(有无rLTs)')); continue
    if any(re.search(r'observation|观察|^(act|adjuvant|chemotherapy)',gn[cv]) for cv in curves):
        excl.append(dict(KMID=k,场景=r['场景'],landmark=r['landmark'],终点=r['终点'],分流='P4治疗交互(ACTvs观察)')); continue
    # 动态转换轨迹：≥3态，或 pos-neg/pos-pos 等转换命名，或组定义含转阴/持续阳
    trans_name=any(re.search(r'(pos|neg)[\s\-]?(pos|neg)',gn[cv]) for cv in curves)
    trans_def=any(re.search(r'转阴|转阳|持续阳|持续阴|转换|clearance|persistent',norm(gdefs[cv])) for cv in curves)
    if len(curves)>=3 or trans_name or trans_def:
        excl.append(dict(KMID=k,场景=r['场景'],landmark=r['landmark'],终点=r['终点'],分流='M2多态/转换轨迹')); continue
    if len(curves)<=1:
        excl.append(dict(KMID=k,场景=r['场景'],landmark=r['landmark'],终点=r['终点'],分流='单组,人工复核: '+' '.join(gnames.values()))); continue
    smap={cv:classify(gnames[cv],gdefs[cv],d1,d2) for cv in curves}
    hi=[cv for cv,v in smap.items() if v==1]; lo=[cv for cv,v in smap.items() if v==0]
    if len(hi)==1 and len(lo)==1:
        rows.append(dict(KMID=k,hicurve=hi[0],locurve=lo[0],
                         hi_label=gnames[hi[0]]+' '+gdefs[hi[0]],lo_label=gnames[lo[0]]+' '+gdefs[lo[0]]))
    else:
        excl.append(dict(KMID=k,场景=r['场景'],landmark=r['landmark'],终点=r['终点'],
                         分流='两态映射歧义,人工复核: '+' | '.join(f'{cv}:{gnames[cv]}->{smap[cv]}' for cv in curves)))

two=pd.DataFrame(rows); ex=pd.DataFrame(excl)
cmpu=cmp[['KMID','evt']].drop_duplicates('KMID').rename(columns={'evt':'cmp_evt'})
m1=two.merge(inc,on='KMID',how='left').merge(nard,on='KMID',how='left').merge(cmpu,on='KMID',how='left')
print("A级130 分流：")
print(ex.分流.str.split('[,:(|]').str[0].value_counts().to_string())
print("\n进入M1两态重建:",len(m1),"KM,",m1.StudyID.nunique(),"研究")
print(m1.pivot_table(index='场景',columns='landmark',values='StudyID',aggfunc=lambda s:s.nunique(),fill_value=0).to_string())
amb=ex[ex.分流.str.contains('歧义|单组')]
print("\n仍需人工复核(",len(amb),"):")
print(amb[['KMID','场景','landmark','分流']].to_string(index=False))
cur[cur.KMID.isin(set(m1.KMID))].to_csv(os.path.join(OUT,'curves.csv'),index=False,encoding='utf-8-sig')
nar[nar.KMID.isin(set(m1.KMID))].to_csv(os.path.join(OUT,'nar.csv'),index=False,encoding='utf-8-sig')
mc=['KMID','StudyID','第一作者','年份','场景','landmark','终点','比较组名称','hicurve','locurve','hi_label','lo_label',
    '各组初始N','M4_高危事件估算','M4_事件层','M4_单研究报告资格','nar_minpts','cmp_evt','HR','HR类型']
m1out=m1[mc].copy()
m1out.columns=['kmid','study_id','author','year','scenario','landmark','endpoint','comparison','hicurve','locurve',
               'hi_label','lo_label','group_n','ev_est_m4','ev_layer','single_report','nar_minpts','cmp_evt','hr','hr_type']
def scen(s):
    s=str(s)
    return 'LARC' if s.startswith('2') else 'MRD' if s.startswith('3') else 'CRLM' if s.startswith('4') else 'mCRC' if s.startswith('5a') else 'other'
def win(l):
    l=str(l)
    if l.startswith('L0'): return 'preop'
    if l.startswith('L1') or l.startswith('L2'): return 'postopMRD'
    if l.startswith('L3'): return 'postACT'
    if l.startswith('L4'): return 'postNAC'
    if l.startswith('L5'): return 'dynamic'
    if l.startswith('L6'): return 'onTx'
    return 'other'
def endp(e):
    e=str(e).lower()
    if ('os' in e) or ('总生存' in e) or ('死亡' in e) or ('overall' in e): return 'OS'
    if any(x in e for x in ['dfs','rfs','pfs','无病','无复','复发','recurrence','disease-free','recurrence-free','progression','无进展']): return 'DFS'
    return 'other'
m1out['scen_code']=m1['场景'].map(scen).values
m1out['win']=m1['landmark'].map(win).values
m1out['endp_code']=m1['终点'].map(endp).values
m1out.to_csv(os.path.join(OUT,'km_meta.csv'),index=False,encoding='utf-8-sig')
ex.to_csv(os.path.join(OUT,'M1_exclusion_log.csv'),index=False,encoding='utf-8-sig')
print("\nNAR min点分布:",m1.nar_minpts.value_counts().sort_index().to_dict())
print("有精确总事件数:",m1.cmp_evt.notna().sum(),"/",len(m1))
print("Bundle ->",OUT)
