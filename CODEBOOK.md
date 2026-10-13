# CODEBOOK — Data and Output Dictionary

All released delimiter files are UTF-8 with BOM (`utf-8-sig`) and open directly
in R, Python, Excel and Stata. Decimal separator is `.`. Times are in
**months** unless the column name says otherwise.

Module numbering: **M1** two-state absolute-risk reconstruction; **M2**
longitudinal trajectories; **M3** risk timing / lead time; **M4**
reconstruction reliability gate (M4a parametric, M4b empirical real-IPD);
**M5** the Shiny calculator; **M6** competing-risk and calibration.
**S1**–**S7** are supplementary tables (S1 QUIPS, S2 CINeMA, S3 reconstructed-vs
-published HR, S4 heterogeneity/influence, S5 external calibration, S6
competing risk, S7 calibration).

---

## 1. Source data

### `data-raw/study_sources.csv` — the 24 included studies (one row each)
| Column | Meaning |
|---|---|
| study_id | internal study identifier (used in every table) |
| first_author, year | citation author/year |
| title, journal, volume, issue, pages_or_elocator | bibliographic fields |
| doi | verified DOI when available (`doi_status=confirmed`) |
| pmid | PubMed ID |
| other_identifier | PII, article number or trial registration when no DOI |
| local_pdf_filename | archived PDF name for audit (not redistributed) |
| assay_class | tumor-informed / tumor-agnostic / methylation / other |
| population | clinical population covered |
| design | e.g. prospective; multicenter |
| median_followup_months | reported median follow-up |
| doi_status | `confirmed` (DOI verified) or `verify_identifier` (use PII/PMID; do not invent) |

### `data/digitized/{KMID}_digitized.csv` — 164 digitised Kaplan–Meier curves
| Column | Meaning |
|---|---|
| KMID | `{study_id}-KM{...}` curve-set identifier |
| id | raw digitiser record id |
| digitized_curve / curve | curve number within the KM set (maps to a group/state) |
| time | time in months |
| St | digitised survival S(t) on 0–1 |

### `data/meta/` — sealed English bundles consumed by the R pipeline
| File | Rows | Meaning |
|---|---|---|
| km_meta.csv | 96 | one row per KM curve-set: study, scenario, landmark, endpoint, high/low curve numbers, group labels, event estimates, published HR, M4 grade (`single_report`, `ev_layer`) |
| curves.csv | 5,724 | digitised points for the M1 two-state sets (`KMID,curve,time,St`) |
| nar.csv | 1,166 | number-at-risk table for M1 (`KMID,curve,grp,time,nrisk`) |
| curves2.csv / nar2.csv | 960 / 242 | digitised points and number-at-risk for the M2 multi-state sets |
| traj_meta.csv | 39 | maps each M2 curve to a `trajectory` |
| M2_traj_coverage.csv | 39 | digitised-point / NAR coverage audit per trajectory |
| study_moderators.csv | 24 | study-level moderators (informed, metastatic, prospective, multicenter, follow-up) |
| M1_exclusion_log.csv | 34 | why each A-level KM set was routed to M2 or excluded (`scenario,landmark,endpoint,routing`) |

---

## 2. Outputs

### `outputs/M1/` — two-state absolute risk
| File | Meaning |
|---|---|
| M1_ipd_reconstructed.csv | 11,481 reconstructed person-level records (`time,status,kmid,curve,ctdna,study_id,scenario,landmark,endpoint,ev_layer`) used in the main analysis |
| M1_absrisk_perstudy.csv | per-study cumulative incidence CIR=1−S(t) at 12/24/36 m by ctDNA (`n,events,CIR,seCIR`) |
| M1_absrisk_pooled.csv | random-effects pooled CIR with CI (`pool,lo,hi`), prediction interval (`pi_lo,pi_hi`), `i2,tau2`, `k` studies, `insufficient` flag |
| M1_contrast_pooled.csv | pooled risk difference `RD` and risk ratio `RR` with CIs |
| M1_core_results_DFS.csv | headline DFS absolute risks / RD / RR / I² |
| M1c_continuous_curves.csv | one-stage Royston–Parmar continuous CIR curves (**also drives the Shiny app**; `spline_df`, per-month `CIR,CIRlo,CIRhi`) |
| M1c_onestage_vs_twostage.csv | one-stage vs two-stage agreement (`abs_diff`) |
| M1_event_counts_reconstructed.csv / M1_reconstruction_diagnostics.csv | reconstructed N, events, censoring, points, NAR count per curve |
| S3_reconstructed_vs_published_HR.csv | 61 contrasts; reconstructed vs published HR (`hr_reconstructed`,`hr_published`,`within_CI`); 52 matched, 50 within CI |
| S4_leave_one_study_out.csv | leave-one-study-out RD influence (`RD_full`,`RD_loo`,`delta`) |
| S4_egger_small_study.csv | Egger test (`Egger_z`,`Egger_p`) |
| S4_meta_regression.csv / S4_subgroup_RD.csv | moderator meta-regression (`beta,p,R2_percent`) and subgroup RD |
| S5_external_calibration.csv | external calibration vs published appendix IPD |

### `outputs/M2/` — longitudinal trajectories
| File | Meaning |
|---|---|
| M2_ipd_trajectory_reconstructed.csv | reconstructed trajectory person-level records |
| M2_trajectory_absolute_risk.csv | per-trajectory CIR with CIs |
| M2_event_counts.csv | event counts / reconstruction quality per trajectory |
| M2_traj_HR_narrative_EN.csv | 17 published trajectory HRs available as HR only (no poolable curve; narrative, not pooled), with `usability_grade` and handling recommendation |

### `outputs/M3/` — timing and lead time
| File | Meaning |
|---|---|
| M3_risk_timing.csv | CIR at 6/12/18/24/36 m and share of events accrued by 12/24 m (`pct_events_by12/24`) |
| M3_leadtime_literature.csv | 9 published ctDNA lead-time values (narrative, no pooling) |

### `outputs/M4/` — reliability gate
| File | Meaning |
|---|---|
| M4_threshold_gate_dual_evidence.csv | combined M4a/M4b grid (`evidence`, `point_estimate_reconstructible`, MAE, coverage) |
| M4_threshold_summary.csv | parametric grid: cells passed, pass rate, median MAE / coverage by points × risk-table × noise |
| M4_threshold_by_events.csv | parametric gate performance by positive-event stratum |
| M4b_empirical_threshold.csv | real-IPD gate performance by event stratum (points × risk table) |

Tolerances: relative bias ≤10 %, MAE ≤0.05 absolute, 95 %-CI coverage 0.90–0.98.

### `outputs/M6/` — competing risk and calibration (S6/S7)
| File | Meaning |
|---|---|
| M6_competing_risk_sensitivity.csv | KM CIR vs competing-risk-adjusted CIR under non-recurrence death hazards λc = 0.5/1/2/3 % per year (`CIR_AJ_*`, difference in percentage points) |
| M6_internal_calibration_pairs.csv | 129 leave-one-study-out observed vs predicted pairs (19 studies) |
| M6_internal_calibration_metrics.csv | calibration slope/intercept with CIs (overall and by ctDNA) |
| M6_calibration_deciles.csv | calibration decile bins (n, predicted, observed) |
| S7_external_calibration_summary.csv | external O−E, E/O, O/E and logit calibration slope/intercept |
| M6_summary.txt | text summary of the module |

### `outputs/S1_QUIPS/` and `outputs/S2_CINeMA/`
- **S1_QUIPS_24_studies.csv** — six QUIPS domains (D1–D6), overall risk, the
  primary and smallest positive-stratum event counts, and a preliminary-rating
  basis column. The companion `*_README.csv` documents the rating rules.
- **S2_CINeMA_confidence.csv** — five core absolute-risk outcomes with k, events,
  pooled %, I², prediction interval, Egger p, the six CINeMA domains and overall
  confidence plus a rationale; `*_methods_README.csv` documents the method.

---

## 3. Coding dictionary

### Clinical scenario (`scenario` / `scen_code`)
| Code | Meaning |
|---|---|
| non_metastatic_postop_MRD (MRD) | non-metastatic (stage I–III) post-surgery / MRD |
| LARC_neoadjuvant (LARC) | locally advanced rectal cancer, peri-neoadjuvant |
| resectable_CRLM (CRLM) | resectable colorectal cancer liver metastases |
| mCRC_systemic (mCRC) | metastatic CRC on systemic therapy |

### Window / landmark (`win`; original `landmark` codes)
| Code | Landmark |
|---|---|
| preop | L0 pre-operative baseline |
| postopMRD | L1 early post-op (≤1 m) and L2 post-op 1–3 m |
| postACT | L3 after adjuvant / systemic treatment |
| postNAC | L4 after neoadjuvant (chemoradiation) |
| dynamic | L5 longitudinal surveillance |
| onTx | L6 early on-treatment response |
| other | L7 treatment interaction / L8 other non-landmark |

### Endpoint (`endp_code`)
`DFS` (covers DFS/RFS/PFS/TTR/recurrence) · `OS` (overall survival) · `other`.

### ctDNA group
`ctdna = 0` negative/clearance reference; `ctdna = 1` positive/non-clearance.
Text labels (`ctdna_label`, `hi_label`, `lo_label`) are the published group names.

### Longitudinal trajectory (`trajectory`)
| Code | Meaning |
|---|---|
| reference_sustained_negative / sustained_neg | reference: never detected |
| clearance / clearance_decline | positive → negative (clearance) |
| newly_emergent_positive / emergent | negative → positive (new emergence) |
| persistent_positive / persistent | persistently positive |
| ever_pos | ever positive during follow-up (study-specific collapse) |

### Reliability gate (`ev_layer`, `evidence`)
- Event strata: `positive_events_le10`, `positive_events_11to20`,
  `positive_events_21to40`, and the wide simulation stratum
  `events_14to284_not_covering_le10`.
- Evidence: `M4a_parametric_simulation_supporting_optimistic` and
  `M4b_real_IPD`; a curve is used only when **both** support reconstruction.
- `usability_grade = C_HR_only_no_curve`: HR reported but no usable curve → not
  pooled (narrative only).

### External validation datasets (`dataset` / `mapped`)
`Li2022_postchemo_stageIII`; `Ogaard2022_L1_postOP_mCRC`;
`Ogaard2022_L3_postACT_mCRC`. Advanced-stage anchors excluded from the
curative-intent pool are flagged `advanced_anchor_excluded_from_pool`.

### Risk-of-bias / confidence
- QUIPS `overall_risk_of_bias`: Low / Moderate / High (domain-level ratings in
  the six D1–D6 columns).
- CINeMA `overall_confidence`: Low / Moderate / High, with the six domain
  columns and a text `rationale`.

### Metric definitions
- **CIR** cumulative incidence of recurrence = 1 − Kaplan–Meier S(t); `seCIR`
  its Greenwood standard error.
- **pool / lo / hi** random-effects pooled CIR and 95 % CI (logit scale);
  **pi_lo / pi_hi** the 95 % prediction interval; **i2**, **tau2** heterogeneity.
- **RD** risk difference (positive minus negative); **RR** risk ratio.
- **hr_reconstructed / hr_published** hazard ratios derived from reconstructed
  curves versus reported; **within_CI** indicates agreement.
- **median_MAE** median absolute error in reconstructed CIR; **median_coverage**
  empirical 95 %-CI coverage; **pass_rate** fraction of grid cells meeting the
  tolerances.
- **CIR_AJ_λ** competing-risk-adjusted CIR under non-recurrence death hazard λ;
  **diff_pp** difference in percentage points.
- **OE / EoverO** observed-over-expected and expected-over-observed ratios;
  **slope / intercept** calibration slope and intercept (1 and 0 = ideal).
