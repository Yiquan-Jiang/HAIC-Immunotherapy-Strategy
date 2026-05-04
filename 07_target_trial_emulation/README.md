# 07 — Stage 3: Target Trial Emulation (TTE)

Two parallel target-trial emulations testing whether **biomarker-guided
deferral of immune-based therapy** ("Adaptive On Demand") yields non-inferior
OS compared with **early combination**.

> **Paper output:** Fig 7 (two-panel KM + ΔRMST forest + IPCW diagnostics).

## The single canonical TTE driver

`tte_IT_R_two_cohorts.R` — **IT_RULES_v2** applied to two cohorts in one run.
Trigger logic (per cohort):

| | **Cohort A**: `cohort_3matched` | **Cohort B**: `cohort_7group_psm02` |
|---|---|---|
| Source                 | matched_06 ∪ all HAIC + I + T concurrent | matched_02 ∪ all HAIC + I concurrent |
| Add-on                 | ICI + antiangiogenic                     | ICI only                             |
| AFP drop trigger       | < −32.5 %                                | < −40 %                              |
| PLR trigger            | > 102.4                                  | > 98.7                               |
| SII trigger            | > 390.9                                  | (NA)                                 |
| NLR trigger            | (NA)                                     | > 2.68                               |
| PIVKA Δ trigger        | < −45.6 %                                | < −51.2 %                            |
| Lymph-node metastasis  | active                                   | active                               |

## Strategy specifications

```
Strategy A (Adaptive On Demand)
  Cycle ≥ 3: trigger on any of: AFP/PIVKA insufficient response, PLR/SII/NLR
             elevation, baseline distant metastasis, baseline lymph-node mets
  Untriggered post-HAIC: trigger on AFP > 20 ng/mL OR PIVKA > 40 mAU/mL

Strategy B (Early Combination)
  Initiate systemic therapy within EARLY_GRACE_DAYS = 14 days of HAIC start,
  regardless of biomarkers.

Dynamic grace: 90 days; IPCW truncation: 99th percentile.
```

## Methodology

**Clone-Censor-Weight (CCW)** framework with **stabilized IPCW**:

```
For each patient:
  Clone into both arms (Strategy A and Strategy B).
  Apply artificial censoring whenever observed treatment violates the
    assigned strategy.
  Re-weight by stabilized inverse-probability-of-censoring weights (IPCW)
    estimated from a pooled person-period logistic regression.
  Fit weighted Cox with robust sandwich SE for HR.
  Compute weighted KM and integrate for RMST at τ ∈ {12, 18, 24, 36} months;
    bootstrap 500× re-estimating IPCW each iteration.
```

**Sensitivity analyses included:**
- Different RMST τ (12, 18, 24, 36 months)
- IPCW truncation at 99th percentile
- E-value for unmeasured confounding

## Layout

```
07_target_trial_emulation/
├── tte_IT_R_two_cohorts.R              ── R core: CCW + IPCW + weighted Cox + RMST
│                                          (writes CSVs to output/step3_tte/IT_RULES_R_two_cohorts/)
├── tte_IT_R_figures_two_cohorts.py     ── Python figures from those CSVs
├── tte_pathway_visualization.R         ── Strategy pathway diagram
├── tte_pathway_visualization_alt_samples.R  Alternate-sample variant
└── generate_tte_flow_drawio_two_cohorts.R   CONSORT-style flow (drawio XML)
```

## Run

```bash
# Convenience runner (recommended)
bash utils_runners/run_tte_two_cohorts.sh

# Manual two-step
Rscript 07_target_trial_emulation/tte_IT_R_two_cohorts.R  data/
python  07_target_trial_emulation/tte_IT_R_figures_two_cohorts.py
```

## Outputs

| Output | Paper element |
|---|---|
| KM curves of Strategy A vs B (Cohort A + B)        | **Fig 7A, B** |
| ΔRMST forest (sensitivity panel)                   | **Fig 7C** |
| IPCW diagnostics (weight distribution, ESS)        | Supplementary |
| TTE flow diagram (clones, censoring)               | **Fig 7** + Supplements |
| HR / RMST / E-value / sensitivity CSVs             | Supplementary tables |
