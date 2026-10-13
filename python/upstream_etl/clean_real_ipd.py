# -*- coding: utf-8 -*-
# =============================================================================
# UPSTREAM, ONE-TIME CLEANING OF PUBLISHED APPENDIX IPD.
# NOT part of the sealed reproduction chain, and the cleaned patient-level
# files are NOT redistributed in this repository.
#
# Reads patient-level appendix data published by Li 2022 and Ogaard 2022
# (third-party data, used for academic validation only). The original .tsv
# files are excluded from this repository under the data-sharing policy; set
# REAL_IPD_DIR below to a local copy reconstructed from the published appendices.
#
# Outputs (ipd_clean_*.csv) are consumed only by M4b (empirical reconstruction
# validation) and the external-calibration checks, and are deliberately not
# committed.
# =============================================================================
import os
import sys
import io
import numpy as np
import pandas as pd

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

# --- EDIT: local folder containing the published appendix files 85_Li2022.tsv
# and 24_Ogaard2022.tsv (not redistributed).
REAL_IPD_DIR = os.environ.get("CTDNA_REAL_IPD", r"<LOCAL_PATH>/real_IPD")
OUT_DIR = os.environ.get("CTDNA_REAL_IPD_OUT", ".")


def km_absrisk(t, e, times):
    """Kaplan-Meier absolute (cumulative) risk at each requested time."""
    t = np.asarray(t, float); e = np.asarray(e, int)
    o = np.argsort(t); t = t[o]; e = e[o]; n = len(t); res = {}
    for T in times:
        tt = np.unique(t[t <= T]); s = 1.0; nr = n
        for ut in tt:
            d = int(np.sum((t == ut) & (e == 1))); c = int(np.sum((t == ut) & (e == 0)))
            if nr > 0:
                s *= (1 - d / nr)
            nr -= (d + c)
        res[T] = 1 - s
    return res


# ---------- Li2022: stage III post-surgery, pre/post-chemo two timepoints ----
li = pd.read_csv(os.path.join(REAL_IPD_DIR, "85_Li2022.tsv"), sep="\t")
# patient-level endpoint (last follow-up / recurrence indicator)
pend = li.groupby("patient_id").agg(
    time_month=("time_month", "max"), event=("event", "max"),
    age=("age", "first"), sex=("sex", "first"), stage=("stage", "first"),
    adj_chemo=("adj_chemo", "first"), CEA=("CEA", "first")).reset_index()


def tp(df, label):
    return df[df["ctDNA_timepoint"] == label][["patient_id", "ctDNA_status"]].rename(
        columns={"ctDNA_status": label})


post = tp(li, "post_chemo"); pre = tp(li, "pre_chemo")
L = pend.merge(post, on="patient_id", how="inner").merge(pre, on="patient_id", how="left")
L = L.rename(columns={"post_chemo": "ctDNA"})
L["dataset"] = "Li2022_postchemo_III"
L["landmark"] = "L3_post_ACT"
L["scenario"] = "non_metastatic_postop_MRD"
L.to_csv(os.path.join(OUT_DIR, "ipd_clean_Li2022_postchemo.csv"), index=False, encoding="utf-8-sig")

# four longitudinal states (pre -> post)
D = L.dropna(subset=["pre_chemo"]).copy()


def state(r):
    a, b = int(r["pre_chemo"]), int(r["ctDNA"])
    return {(0, 0): "sustained_neg", (1, 0): "clearance",
            (0, 1): "emergent", (1, 1): "persistent"}[(a, b)]


D["trajectory"] = D.apply(state, axis=1)
D.to_csv(os.path.join(OUT_DIR, "ipd_clean_Li2022_dynamic.csv"), index=False, encoding="utf-8-sig")

# ---------- Ogaard2022: stage IV mCRC, multiple timepoints -------------------
og = pd.read_csv(os.path.join(REAL_IPD_DIR, "24_Ogaard2022.tsv"), sep="\t")
og["ctDNA_timepoint"] = og["ctDNA_timepoint"].replace({"post_ACT": "post-ACT"})
# patient-level endpoint: last follow-up (max time), event = any recurrence
op = og.groupby("patient_id").agg(time_month=("time_month", "max"), event=("event", "max")).reset_index()
op["dataset"] = "Ogaard2022_IV"
op["scenario"] = "mCRC_systemic"


def og_landmark(label, lm):
    z = og[og["ctDNA_timepoint"] == label].sort_values("time_month").groupby("patient_id").tail(1)[
        ["patient_id", "ctDNA_status", "time_month"]].rename(
        columns={"ctDNA_status": "ctDNA", "time_month": "lm_time"})
    m = op.merge(z, on="patient_id", how="inner")
    # follow-up after landmark: exclude recurrence at/before landmark; time from landmark
    m = m[m["time_month"] > m["lm_time"]].copy()
    m["fu_from_lm"] = m["time_month"] - m["lm_time"]
    m["landmark"] = lm
    m.to_csv(os.path.join(OUT_DIR, f"ipd_clean_Ogaard2022_{lm}.csv"), index=False, encoding="utf-8-sig")
    return m


o_postop = og_landmark("post-OP", "L1_postOP")
o_postact = og_landmark("post-ACT", "L3_postACT")
og.sort_values(["patient_id", "time_month"]).to_csv(
    os.path.join(OUT_DIR, "ipd_clean_Ogaard2022_serial.csv"), index=False, encoding="utf-8-sig")


def summ(name, df, tcol="time_month", ecol="event", g="ctDNA", times=(12, 24, 36)):
    print("=" * 60)
    print(name, " n patients", df[g].notna().sum() if g in df else len(df))
    if g in df:
        for v in sorted(df[g].dropna().unique()):
            s = df[df[g] == v]; ar = km_absrisk(s[tcol], s[ecol], times)
            print(f"  ctDNA={int(v)}: n={len(s)}, events={int(s[ecol].sum())}, "
                  f"group fraction={len(s)/len(df):.2f}, absolute recurrence risk@" +
                  ", ".join(f"{T}m={ar[T]:.3f}" for T in times))


summ("Li2022 post-ACT landmark (stage III)", L)
print("\nLi2022 four-state distribution:"); print(D["trajectory"].value_counts().to_string())
for _, r in D.groupby("trajectory"):
    ar = km_absrisk(r["time_month"], r["event"], (12, 24, 36))
    print(f"  {r['trajectory'].iloc[0]}: n={len(r)}, events={int(r['event'].sum())}, " +
          ", ".join(f"{T}m={ar[T]:.3f}" for T in (12, 24, 36)))
summ("Ogaard2022 post-OP landmark (stage IV)", o_postop, tcol="fu_from_lm")
summ("Ogaard2022 post-ACT landmark (stage IV)", o_postact, tcol="fu_from_lm")
print("\nOgaard serial long-table rows", len(og), "patients", og.patient_id.nunique(),
      "timepoint types", og.ctDNA_timepoint.value_counts().to_dict())
