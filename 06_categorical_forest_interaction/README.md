# 06 — Stage 2: categorical interaction forest plots

> **Note**: although the upstream source folder was named `publication_figures/`,
> these scripts specifically draw **immunotherapy-benefit × categorical-variable
> interaction** forest plots (i.e. the categorical companion to the continuous
> RCS interaction analysis in folder `05_rcs_interaction/`).

Apply the cut-points identified in Stage 2 (RCS interaction) to generate
publication-ready interaction forest plots showing subgroup heterogeneity of
treatment effect.

> **Paper output:** Fig 5 / Fig 6 forest panels + supplementary forests.

## Comparisons covered

| Script | Pair (PSM ID) | Treatment vs. Reference |
|---|---|---|
| `01_publication_figures.py` | PSM02 | HAIC → ICI vs HAIC alone |
| `02_publication_figures_ids06_IplusT.py` | PSM06 | HAIC → ICI + Anti-angio vs HAIC alone |
| `03_publication_figures_iptw_psm02.py` | PSM02 (IPTW + PSM) | HAIC → ICI — combined panel |
| `03_publication_figures_ids05_IplusT_concurrent.py` | PSM05 | HAIC + ICI + Anti-angio (concurrent) vs HAIC alone |
| `04_publication_figures_iptw_psm06_IplusT.py` | PSM06 (IPTW + PSM) | HAIC → ICI + Anti-angio — combined panel |
| `05_publication_figures_iptw_psm05_IplusT_concurrent.py` | PSM05 (IPTW + PSM) | concurrent triplet — combined panel |

## Statistical specifications

- For each subgroup, fit a Cox model and compute HR + 95 % CI.
- ΔRMST (restricted-mean-survival-time difference) with 300-replicate bootstrap CI.
- **Likelihood-ratio test** of treatment × subgroup interaction → P_interaction.
- **Landmark sensitivity** at 42-day landmark to address immortal-time bias.

## Outputs

| Output | Paper element |
|---|---|
| `Fig_PSM02_forest.{pdf,png}` | **Fig 5** (or related Stage-2 forest) |
| `Fig_PSM06_forest.{pdf,png}` | **Fig 5** (or related Stage-2 forest) |
| `Fig_PSM05_forest.{pdf,png}` | Supplements |
| Subgroup analysis tables (HR, ΔRMST, P_interaction) | Supplementary tables |
