# 06 — Stage 2: categorical interaction forest plots

Apply the cut-points identified in Stage 2 (RCS interaction, see folder `05_rcs_interaction/`)
to generate **publication-ready interaction forest plots** showing subgroup
heterogeneity of treatment effect, focusing on:

- **PSM02**: HAIC alone vs HAIC → ICI
- **PSM05**: HAIC alone vs HAIC + ICI + Anti-angio concurrent
- **PSM06**: HAIC alone vs HAIC → ICI + Anti-angio (delayed triplet)

> **Paper output:** Fig 5 forest panels + supplementary forests.

## Layout

```
06_categorical_forest_interaction/
└── final_publication/        ← Publication-ready, IPTW + PSM versions
    ├── 01_publication_figures.py                         PSM02
    ├── 02_publication_figures_ids06_IplusT.py            PSM06
    ├── 03_publication_figures_iptw_psm02.py              IPTW+PSM02 combined panel
    ├── 03_publication_figures_ids05_IplusT_concurrent.py PSM05
    ├── 04_publication_figures_iptw_psm06_IplusT.py       IPTW+PSM06
    └── 05_publication_figures_iptw_psm05_IplusT_concurrent.py IPTW+PSM05
```

## Statistical specifications

- For each subgroup, fit a Cox model and compute HR + 95% CI.
- ΔRMST (restricted mean survival time difference) with 300-replicate bootstrap CI.
- **Likelihood-ratio test** of treatment × subgroup interaction for P_interaction.
- **Landmark sensitivity** at 42-day landmark to address immortal-time bias.

## Outputs

| Output | Paper element |
|---|---|
| `Fig_PSM02_forest.{pdf,png}` | **Fig 5** (or related Stage-2 forest) |
| `Fig_PSM06_forest.{pdf,png}` | **Fig 5** (or related Stage-2 forest) |
| `Fig_PSM05_forest.{pdf,png}` | Supplements |
| Subgroup analysis tables (HR, ΔRMST, P_interaction) | Supplementary tables |
