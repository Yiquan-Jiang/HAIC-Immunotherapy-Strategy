# 05 — Stage 2: RCS × treatment interaction (continuous biomarkers)

Restricted cubic spline (RCS) interaction Cox models on each pairwise
comparison of HAIC alone vs. each combination strategy, to identify
**continuous biomarkers that modify the relative efficacy** of adding
systemic therapy to HAIC.

> **Paper output:** Fig 6 (RCS interaction matrix) + supplementary panels.

> **Authoritative script index:** see [`FINAL_SCRIPTS.md`](./FINAL_SCRIPTS.md)
> in this folder. Two parallel methodological routes (PSM 1:1 unweighted Cox
> vs. composite-cohort IPTW-weighted Cox) act as **independent sensitivity
> analyses of one another**.

## Final scripts (per FINAL_SCRIPTS.md)

|  | **Single-indicator combined plots**<br/>(one figure per biomarker) | **8 × 5 matrix overview**<br/>(one large figure: 8 biomarkers × 5 timepoints) |
|---|---|---|
| **Route A — PSM 1:1, unweighted Cox**<br/>6 pair-wise contrasts | `RCS_PSM_dual_timescale.R` | `RCS_PSM_matrix_panel.R` |
| **Route B — composite cohort, IPTW-weighted Cox**<br/>5 arms vs HAIC alone | `afp_pivka_composite/01_rcs_afp_pivka_composite.R` | `afp_pivka_composite/02_rcs_matrix_panel.R` |

## Layout

```
05_rcs_interaction/
├── FINAL_SCRIPTS.md                   ← Authoritative script index
├── build_cohort_psm.py                ← Prerequisite for Route A
├── RCS_PSM_dual_timescale.R           ── Route A: single-indicator
├── RCS_PSM_matrix_panel.R             ── Route A: 8×5 matrix
├── afp_pivka_composite/
│   ├── 00_build_composite_cohorts.py  ← Prerequisite for Route B
│   ├── 00a_extract_pre_it_labs.py     ── Pre-IT lab extraction helper
│   ├── 01_rcs_afp_pivka_composite.R   ── Route B: single-indicator
│   └── 02_rcs_matrix_panel.R          ── Route B: 8×5 matrix
└── publication_figures/               ← Final 5×6 publication panels
    ├── README.md
    ├── make_publication_figures_iptw.R   (Route B → main panel)
    └── make_publication_figures_psm.R    (Route A → sensitivity panel)
```

## Pipeline

```bash
# Build cohorts (run first)
python  build_cohort_psm.py                                # Route A cohorts
python  afp_pivka_composite/00_build_composite_cohorts.py  # Route B cohorts

# Route A: PSM + unweighted Cox
Rscript RCS_PSM_dual_timescale.R                           # single-indicator
Rscript RCS_PSM_matrix_panel.R                             # 8×5 matrix

# Route B: IPTW (CBPS / ridge) weighted Cox
Rscript afp_pivka_composite/01_rcs_afp_pivka_composite.R ALL
Rscript afp_pivka_composite/02_rcs_matrix_panel.R         ALL

# Publication-ready 5×6 panels (assemble Route B for main, Route A for sensitivity)
RMS_RCS_N_BOOT=300 Rscript publication_figures/make_publication_figures_iptw.R
RMS_RCS_N_BOOT=300 Rscript publication_figures/make_publication_figures_psm.R
```

Optional environment variables: `RMS_RCS_NK` (default = 3), `RMS_RCS_N_BOOT` (default = 200).

## Output locations

- Route A → `output/step1_rcs_interaction/psm/pair_XX/{landmark,total_os}/`
- Route B → `output/step1_rcs_interaction/relaxed_dual/<cohort_key>/{landmark,total_os}/`
- Publication panels → `output/step1_rcs_interaction/publication/{iptw,psm}/{landmark,total_os}/`

## Statistical specifications

- **Effect modifier search**: 8 candidate continuous biomarkers
  (AFP, PIVKA-II, PIV, SII, NLR, PLR, MONOCYTE, ALBI) × 5 timepoints
  (Baseline, Pre-HAIC-3, Pre-IT, Pre-HAIC-3 change rate, Pre-IT change rate)
  fitted as `rcs(x, nk = 3)`.
- **Model**: `coxph(Surv(time, event) ~ trt + rcs(x, 3) + trt:rcs(x, 3) + covariates)`.
- **Weighting**: Route A — PSM 1:1 (unweighted); Route B — IPTW via ridge logistic glmnet (α = 0, 5-fold CV), stabilized weights clipped to [0.05, 0.95].
- **Bootstrap**: 200 replicates (300 for publication panels) → 95 % confidence band on HR(x).
- **Dual timescale**: 42-day landmark survival (primary, IT-appropriate) + total OS (sensitivity).
- **Crossing-point detection**: HR(x) crossing 1.0, 0.85, 0.7 → triggers the cut-points used in the on-demand decision rule (Stage 3 TTE).
