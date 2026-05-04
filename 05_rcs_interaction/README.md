# 05 — Stage 2: RCS × treatment interaction (continuous biomarkers)

Restricted cubic spline (RCS) interaction Cox models on each of the six
pairwise comparisons of HAIC alone vs each combination strategy, to identify
**continuous biomarkers that modify the relative efficacy** of adding
systemic therapy to HAIC.

> **Paper output:** Fig 6 (RCS interaction matrix) + supplementary panels.

## Two complementary versions

```
05_rcs_interaction/
├── tte_project_version/      ← Comprehensive: builds all 6 pairs at once,
│                                runs unified analysis, produces matrix panel
├── group7_project_version/   ← Per-pair: one R + one Python per cohort,
│                                useful for debugging individual pairs
├── afp_pivka_composite/      ← Composite biomarker analysis
│                                (AFP, PIVKA-II combined)
└── publication_figures/      ← Final figure assembly (PSM + IPTW)
```

## Pipeline (recommended)

```
build_all_pairs_cohorts.py    ── Build 6 IPTW pairwise cohorts
        ↓
rcs_all_pairs_dual_timescale.R── Main analysis: Surv(time, event) ~ trt × rcs(x, nk=3)
                                  Both 42-day landmark + total OS timescales
        ↓
build_cohort_psm.py           ── Build 6 PSM pairwise cohorts
        ↓
RCS_PSM_dual_timescale.R      ── Sensitivity analysis (PSM, unweighted)
        ↓
RCS_PSM_matrix_panel.R        ── Assemble Fig 6 matrix panel
        ↓
publication_figures/
    make_publication_figures_iptw.R   ── Final IPTW panels
    make_publication_figures_psm.R    ── Final PSM panels
```

## Statistical specifications

- **Effect modifier search**: 14 candidate continuous biomarkers
  (pre-IT NLR, PLR, SII, PIV, ALBI, AFP, PIVKA-II, ALT, AST, total bilirubin, albumin, platelet, lymphocyte, neutrophil)
  fitted as `rcs(x, nk = 3)`.
- **Model**: `coxph(Surv(time, event) ~ trt + rcs(x, 3) + trt:rcs(x, 3) + covariates)`.
- **Weighting**: relaxed-dual IPTW (CBPS) for primary; PSM-matched (unweighted) for sensitivity.
- **Bootstrap**: 200 replicates → 95% confidence band on HR(x).
- **Dual timescale**: 42-day landmark survival (primary) + total OS (sensitivity).
- **Crossing point detection**: HR(x) = 1.0, 0.85, 0.7 → triggers candidate cut-points for the on-demand decision rule used in Stage 3.

## Outputs

| Output | Paper element |
|---|---|
| Per-biomarker HR(x) curves with 95% CI band | **Fig 6** + supplements |
| Interaction P (overall + non-linear, ANOVA Type-II) per biomarker × pair | Supplementary tables |
| Crossing-point CSVs (used to define triggers in Stage 3) | Methods §"Construction of biomarker-guided on-demand decision rules" |
