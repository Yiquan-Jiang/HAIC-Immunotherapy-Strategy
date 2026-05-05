# Publication figures — RCS non-linear interaction

Two independent causal-adjustment methods, generated as parallel figure sets.

| Role            | Method                                             | Script                              | Output dir |
|-----------------|----------------------------------------------------|-------------------------------------|------------|
| **Primary**     | Route B: composite cohort + IPTW-weighted Cox      | `make_publication_figures_iptw.R`   | `iptw/`    |
| **Sensitivity** | Route A: PSM 1:1 match + unweighted Cox            | `make_publication_figures_psm.R`    | `psm/`     |

Re-run:
```bash
RMS_RCS_N_BOOT=300 Rscript make_publication_figures_iptw.R   # primary
RMS_RCS_N_BOOT=300 Rscript make_publication_figures_psm.R    # sensitivity
```

## Comparisons

| ID       | Treatment arm   | Control arm | IPTW cohort file                  | PSM cohort file                               |
|----------|-----------------|-------------|-----------------------------------|-----------------------------------------------|
| THEN_I   | HAIC_then_I     | HAIC_alone  | `composite_THEN_I_cohort.csv`     | `cohort_02_HAIC_alone_vs_HAIC_then_I.csv`     |
| THEN_IT  | HAIC_then_I+T   | HAIC_alone  | `composite_THEN_IT_cohort.csv`    | `cohort_06_HAIC_alone_vs_HAIC_then_IT.csv`    |

## Biomarkers & timepoints

- Biomarkers (rows, 5): AFP, PIVKA, SII, PLR, NLR
- Timepoints (columns, 3): Baseline, Pre-IT, Pre-IT Change Rate (%)

## Files

```
iptw/  or  psm/
├── landmark/         ← 42-day landmark OS (primary, IT-appropriate)
│   Fig_Main_RCS_5indicators.{pdf,png,tiff}        5×6 landscape, 17.5×10 in
│   Fig_Supp_RCS_Baseline.{pdf,png,tiff}           5×2
│   Fig_Supp_RCS_PreIT.{pdf,png,tiff}              5×2
│   Fig_Supp_RCS_PreIT_ChangeRate.{pdf,png,tiff}   5×2
│   (output filenames defined in the R scripts)
└── total_os/         ← Total OS (sensitivity)
```

- `.pdf` = vector (submission master)
- `.png` = 600 DPI raster
- `.tiff` = 600 DPI LZW compressed

## Main figure layout (5 × 6 landscape)

```
            ┌ Panel A ──────────────────────────┐ ┌ Panel B ──────────────────────────┐
            │ HAIC then I  vs  HAIC alone       │ │ HAIC then I+T vs HAIC alone       │
            │ Baseline  |  Pre-IT  |  Change %  │ │ Baseline  |  Pre-IT  |  Change %  │
      AFP   │   [ ]        [ ]        [ ]       │ │   [ ]        [ ]        [ ]       │
      PIVKA │   [ ]        [ ]        [ ]       │ │   [ ]        [ ]        [ ]       │
      SII   │   [ ]        [ ]        [ ]       │ │   [ ]        [ ]        [ ]       │
      PLR   │   [ ]        [ ]        [ ]       │ │   [ ]        [ ]        [ ]       │
      NLR   │   [ ]        [ ]        [ ]       │ │   [ ]        [ ]        [ ]       │
            └───────────────────────────────────┘ └───────────────────────────────────┘
```

## Panel elements

- **Blue curve**: HR point estimate (treatment vs. HAIC alone) across biomarker range
- **Shaded band**: bootstrap 95% CI (B = 300)
- **Reference lines**: HR = 1.0 (solid grey), 0.85 (blue dot-dash), 0.70 (red dot-dash)
- **Labelled dots**: biomarker value where the HR curve crosses each threshold
- **Green / amber tints**: HR < 1 favors treatment; HR > 1 favors HAIC alone
- **Rug**: observed biomarker distribution
- **Subtitle**: `P[int]` = overall interaction p, `P[nl]` = non-linear component p

## Model details

### Route B (IPTW — primary)
- Propensity score: ridge logistic regression (glmnet, α = 0, 5-fold CV) on PS covariates
- Stabilized weights, clipped to [0.05, 0.95]
- Cox: `Surv(time, death) ~ trt * rcs(biomarker, knots = 3)`, `weights = sw`
- PS covariates — static: ALBI, INR, PLT, tumor size/count, PVTT, HVTT, IVC/RA, ascites,
  log(AFP/PIVKA), metastasis, LN status, neutrophils/lymphocytes/monocytes
- PS covariates — dynamic (change-rate panels add): pre-IT AFP/PIVKA/ALBI/neut/AST/ALT,
  AFP/PIVKA change rate

### Route A (PSM — sensitivity)
- 1:1 nearest-neighbor matching on the same PS covariates
- Cox (unweighted): `Surv(time, death) ~ trt * rcs(biomarker, knots = 3)`

Both routes: static biomarkers (AFP/PIVKA/SII) `log1p`-transformed, axis labels show raw values;
change-rate variables untransformed, x-axis in %.
