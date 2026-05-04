# 07 — Stage 3: Target Trial Emulation (TTE)

Two parallel target-trial emulations testing whether **biomarker-guided
deferral of immune-based therapy** yields non-inferior OS compared with
**early combination**.

> **Paper output:** Fig 7 (two-panel KM + ΔRMST forest + IPCW diagnostics).

## Two parallel cohorts

| Cohort | Question | Biomarker rule |
|---|---|---|
| **Cohort A — ICI-only emulation** | HAIC → ICI vs HAIC + ICI early concurrent | **NLR-based** triggers (Rule 1-5 + ExRule 2-3) |
| **Cohort B — ICI + Anti-angio emulation** | HAIC → ICI + Anti-angio vs HAIC + ICI + Anti-angio early concurrent | **PIV-based** triggers |

## Strategy specifications

```
Strategy A (Dynamic / on-demand)
  Add ICI/Anti-angio only if biomarker rule fires:
    Rule 1: Vp3/4 PVTT, distant metastasis, HVTT, IVC-RA thrombus, lymph-node mets
    Rule 2: baseline or pre-HAIC-3 NLR ≥ 2.5
    Rule 3: baseline largest tumor diameter > 13 cm
    Rule 4: baseline PIVKA-II > 12,000 mAU/mL
    Rule 5: baseline AFP < 20 ng/mL
  Exemptions (cycle ≥ 4):
    ExRule 2: NLR persistently < 2.5 → continue observation
    ExRule 3: AFP decline > 50% from baseline → continue observation

Strategy B (Early combination)
  Initiate systemic therapy within 14 days of HAIC initiation, regardless of biomarkers.
```

## Methodology

**Clone-Censor-Weight (CCW)** framework (Hernán & Robins; Cain et al.):

```
For each patient:
  Clone into both arms (Strategy A and Strategy B).
  Apply artificial censoring whenever observed treatment violates the assigned strategy.
  Re-weight by stabilized inverse-probability-of-censoring weights (IPCW)
    estimated from a pooled person-period logistic regression.
  Fit weighted Cox with robust sandwich SE for HR.
  Compute weighted KM and integrate for RMST; bootstrap 500× re-estimating IPCW.
```

**Sensitivity analyses included:**
- IPTW × IPCW double weighting
- Weight truncation at 1st/99th percentiles
- Different grace periods (7d / 14d / 21d / 28d)
- E-value for unmeasured confounding

## Pipeline

```
analysis_ready.csv  +  matched_ids_*.csv
       │
       ▼
tte_piv_R_core_cohort_7group_psm02.R   ── Core CCW + IPCW + weighted Cox + RMST
tte_IT_R_core_cohort_3matched.R        ── 3-cohort matched variant
tte_IT_R_two_cohorts.R                 ── Two-cohort version (ICI / ICI+T)
       │  (writes CSVs: HR, RMST, sensitivity, KM, risk-table)
       ▼
tte_nlr_R_figures.py                   ── Publication figures (NLR rule)
tte_IT_R_figures.py                    ── Single-cohort figures
tte_IT_R_figures_two_cohorts.py        ── Two-cohort figures (ICI vs ICI+T)
       │
       ▼
tte_pathway_visualization.R            ── Strategy pathway diagram
generate_tte_flow_drawio*.R            ── CONSORT-style flow diagrams (drawio XML)
```

## Outputs

| Output | Paper element |
|---|---|
| KM curves of Strategy A vs B (Cohort A + Cohort B) | **Fig 7A, B** |
| ΔRMST forest (sensitivity panel)                   | **Fig 7C** |
| IPCW diagnostics (weight distribution, ESS)        | Supplementary |
| TTE flow diagram (clones, censoring)               | **Fig 7** + Supplements |
| HR / RMST / E-value / sensitivity CSVs             | Supplementary tables |
