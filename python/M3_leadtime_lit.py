# -*- coding: utf-8 -*-
# M3: published ctDNA lead-time compilation (publicly reported median/mean;
# narrative, no meta-analytic pooling). Self-contained; writes outputs/M3/.
import os
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "outputs", "M3")
os.makedirs(OUT, exist_ok=True)
rows=[
 ("Reinert T",2019,"Stage I-III CRC, Signatera","mean 8.7 (range 0.8-16.5)","vs radiological imaging","PMC6976395/PMC6512280"),
 ("Tumor-informed WGS cohort",2024,"Stage III CRC","median 8.7","vs radiological recurrence","PubMed 39316995"),
 ("I-IV resected CRC cohort",2024,"Stage I-IV","median 10.1","vs imaging-confirmed recurrence","PMC10822076"),
 ("Tumor-informed MRD cohort",2024,"Resected CRC","median 3.6 (IQR 1.5-5.5)","vs radiologic disease","PMC10805428"),
 ("Early meta-analysis",2019,"CRC, 4 studies","mean 7.1","vs conventional surveillance","meta-analysis (value circulation)"),
 ("Liquid biopsy early-detection",2025,"CRC recurrence","median 8.2; sensitivity 89% vs CEA 54%; NPV 94%","vs radiography/CEA","cjcpt.org/799"),
 ("Prospective multicentre serial",2024,"Stage II/III CRC","mean 5.01","vs imaging","Semantic Scholar 5ddb"),
 ("Scholer LV",2020,"CRC (matrix 04)","median 9.4 (95% CI 7.8-11)","vs CEA","matrix 04 Study49"),
 ("Systematic review/meta (17 studies)",2025,"CRC / LARC / mCRC","CRC 1.0-11.5; LARC 10.2; mCRC 3.1-9.4","ctDNA detection to radiographic recurrence","PMC12010414"),
]
df=pd.DataFrame(rows,columns=["study","year","population","lead_time_months","comparator","source"])
df.to_csv(os.path.join(OUT,"M3_leadtime_literature.csv"),index=False,encoding="utf-8-sig")
print("written",len(df))
