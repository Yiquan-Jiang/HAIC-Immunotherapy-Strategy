# HAIC-Immunotherapy-Strategy

**Code repository for the manuscript:**

> *Immunotherapy Following Induction Arterial Chemotherapy for Unresectable Hepatocellular Carcinoma: A Biomarker-Guided Strategy*

A causal-inference framework that compares **seven first-line treatment strategies** in 3,885 patients with unresectable HCC treated with induction FOLFOX-HAIC, and develops a **biomarker-guided on-demand decision rule** for adding immunotherapy / antiangiogenic therapy (validated in two parallel **target-trial emulations**).

---

## Highlights

- Decision model to guide sequential immunotherapy after induction HAIC.
- Triggers: AFP / PIVKA-II changes during HAIC plus PLR / NLR / SII inflammatory indices.
- HAIC + ICI may match HAIC + ICI + antiangiogenic triplet survival.
- Delayed sequential: comparable or numerically longer OS than concurrent.

---

## Seven treatment groups

| ID | Group label                | Description                                  |
| -- | -------------------------- | -------------------------------------------- |
| 1  | `HAIC alone`               | Induction HAIC only                          |
| 2  | `HAIC + ICI concurrent`    | HAIC + immunotherapy concurrently            |
| 3  | `HAIC → ICI`               | Delayed sequential immunotherapy             |
| 4  | `HAIC + Anti-angio conc.`  | HAIC + antiangiogenic concurrently           |
| 5  | `HAIC → Anti-angio`        | Delayed sequential antiangiogenic            |
| 6  | `HAIC + ICI + Anti-angio` (triplet) | All three concurrently              |
| 7  | `HAIC → ICI + Anti-angio`  | Delayed sequential triplet                   |

---

## Analytic framework (3 sequential stages of increasing causal rigor)

```
Stage 1  Confounder-adjusted comparisons
         IPTW (CBPS) + 21 pairwise PSM + Overlap Weighting
                       │
                       ▼
Stage 2  Biomarker × treatment interaction
         RCS continuous interaction + categorical forest plots
                       │
                       ▼
Stage 3  Target Trial Emulation (CCW + IPCW)
         Dynamic biomarker-guided strategy vs. early combination
```

---

## Repository layout

```
HAIC-Immunotherapy-Strategy/
├── 00_data_preparation/                   ← Build analysis-ready dataset (7 groups)
│   └── step0_prepare_data.py
│
├── 01_psm_iptw_overall_survival/          ← Stage 1
│   ├── step3_psm_analysis.R                  one-to-one nearest-neighbor PSM
│   ├── step3b_psm_vs_template.R              standardized PSM vs template
│   ├── step4_km_curves.py                    KM curves (overall, unweighted / IPTW)
│   ├── step4b_km_template_matched.py         KM curves of PSM-matched pairs
│   ├── step5_forest_plot.py                  HR / RMST forest summary
│   ├── step5b_forest_vs_IT_concurrent.py     vs HAIC + ICI concurrent
│   ├── step5c_forest_vs_HAIC_alone.py        vs HAIC alone (21 pairwise)
│   └── step6_tables_and_loveplots.R          baseline tables + love plots (balance)
│
├── 02_subgroup_overlap_weighting/         ← Subgroup analysis (overlap weighting)
│   ├── step7_ow_balance_table.R              OW balance tables
│   ├── step7_ow_balance_figure.py            OW love-plot
│   ├── step7_subgroup_ow.R                   Cox HR per subgroup (OW-weighted)
│   ├── step7_subgroup_analysis.py            ΔRMST + interaction tests
│   └── step7_subgroup_plots.py               Subgroup forest panels
│
├── 03_swimmer_plot/                       ← Per-patient timelines
│   ├── swimmer_plot_7groups.R                Swimmer plots for the 7 groups
│   └── plot_haic_then_i_to_target_interval.R HAIC → ICI interval / target-window adherence
│
├── 04_biomarker_dynamics/                 ← Longitudinal AFP/PIVKA-II/inflammatory indices
│   └── psm_afp_pivka_dynamics.py             Trajectories within PSM-matched pairs
│
├── 05_rcs_interaction/                    ← Stage 2 (continuous biomarkers)
│   │   See FINAL_SCRIPTS.md for the authoritative 4-script index
│   │   (Route A: PSM unweighted; Route B: composite-cohort IPTW)
│   ├── FINAL_SCRIPTS.md
│   ├── build_cohort_psm.py                   Prerequisite — Route A cohorts
│   ├── RCS_PSM_dual_timescale.R              Route A: single-indicator RCS
│   ├── RCS_PSM_matrix_panel.R                Route A: 8 × 5 matrix panel
│   ├── afp_pivka_composite/
│   │   ├── 00_build_composite_cohorts.py     Prerequisite — Route B cohorts
│   │   ├── 00a_extract_pre_it_labs.py        Pre-IT lab extraction helper
│   │   ├── 01_rcs_afp_pivka_composite.R      Route B: single-indicator RCS
│   │   └── 02_rcs_matrix_panel.R             Route B: 8 × 5 matrix panel
│   └── publication_figures/                  5 × 6 publication panels
│       ├── README.md
│       ├── make_publication_figures_iptw.R   Route B (primary)
│       └── make_publication_figures_psm.R    Route A (sensitivity)
│
├── 06_categorical_forest_interaction/     ← Stage 2 (categorical biomarkers)
│   │   "Immunotherapy benefit × categorical-variable" interaction forest plots
│   ├── 01_publication_figures.py                  PSM02: HAIC → ICI
│   ├── 02_publication_figures_ids06_IplusT.py     PSM06: HAIC → ICI + Anti-angio
│   ├── 03_publication_figures_iptw_psm02.py       PSM02 IPTW + PSM combined
│   ├── 03_publication_figures_ids05_IplusT_concurrent.py  PSM05 concurrent triplet
│   ├── 04_publication_figures_iptw_psm06_IplusT.py
│   └── 05_publication_figures_iptw_psm05_IplusT_concurrent.py
│
├── 07_target_trial_emulation/             ← Stage 3
│   │   Single canonical driver: tte_IT_R_two_cohorts.R (IT_RULES_v2)
│   │   handles both cohorts in one invocation:
│   │     - Cohort A: cohort_3matched     — ICI + antiangiogenic add-on
│   │     - Cohort B: cohort_7group_psm02 — ICI-only add-on
│   │   Adaptive On Demand vs Early Combination, CCW + stabilized IPCW.
│   ├── tte_IT_R_two_cohorts.R                R core: CCW + IPCW + Cox + RMST
│   ├── tte_IT_R_figures_two_cohorts.py       Figures (KM, ΔRMST, IPCW)
│   ├── tte_pathway_visualization.R           Strategy pathway diagram
│   ├── tte_pathway_visualization_alt_samples.R
│   └── generate_tte_flow_drawio_two_cohorts.R   CONSORT-style flow (drawio XML)
│
└── docs/                                  ← Original project READMEs (preserved)
    ├── README_TTE_project_original.md
    ├── README_group7_project_original.md
    └── ORGANIZATION_REPORT_original.md
```

---

## Quick start

### Environment

```bash
# Python (>=3.11)
pip install pandas numpy matplotlib lifelines scipy statsmodels seaborn

# R (>=4.3)
install.packages(c(
  "dplyr", "tidyr", "stringr",
  "survival", "survey", "rms", "MatchIt", "WeightIt", "cobalt",
  "boot", "glmnet", "ggplot2", "gtsummary", "gridExtra"
))
```

### Pipeline

```bash
# Stage 0: prepare analysis-ready dataset (requires raw inputs)
python  00_data_preparation/step0_prepare_data.py

# Stage 1: PSM / IPTW + KM + tables
Rscript 01_psm_iptw_overall_survival/step3_psm_analysis.R
python  01_psm_iptw_overall_survival/step4_km_curves.py
python  01_psm_iptw_overall_survival/step4b_km_template_matched.py
python  01_psm_iptw_overall_survival/step5_forest_plot.py
Rscript 01_psm_iptw_overall_survival/step6_tables_and_loveplots.R

# Subgroup analysis with overlap weighting
Rscript 02_subgroup_overlap_weighting/step7_ow_balance_table.R
Rscript 02_subgroup_overlap_weighting/step7_subgroup_ow.R
python  02_subgroup_overlap_weighting/step7_subgroup_analysis.py
python  02_subgroup_overlap_weighting/step7_subgroup_plots.py

# Swimmer plots + HAIC → ICI interval visualization
Rscript 03_swimmer_plot/swimmer_plot_7groups.R
Rscript 03_swimmer_plot/plot_haic_then_i_to_target_interval.R

# Longitudinal biomarker dynamics
python  04_biomarker_dynamics/psm_afp_pivka_dynamics.py

# Stage 2a: RCS interaction — see 05_rcs_interaction/FINAL_SCRIPTS.md
python  05_rcs_interaction/build_cohort_psm.py                                # Route A
python  05_rcs_interaction/afp_pivka_composite/00_build_composite_cohorts.py  # Route B
Rscript 05_rcs_interaction/RCS_PSM_dual_timescale.R                           # Route A: single-indicator
Rscript 05_rcs_interaction/RCS_PSM_matrix_panel.R                             # Route A: 8×5 matrix
Rscript 05_rcs_interaction/afp_pivka_composite/01_rcs_afp_pivka_composite.R ALL  # Route B: single-indicator
Rscript 05_rcs_interaction/afp_pivka_composite/02_rcs_matrix_panel.R         ALL  # Route B: 8×5 matrix
RMS_RCS_N_BOOT=300 Rscript 05_rcs_interaction/publication_figures/make_publication_figures_iptw.R
RMS_RCS_N_BOOT=300 Rscript 05_rcs_interaction/publication_figures/make_publication_figures_psm.R

# Stage 2b: categorical interaction forest panels
python  06_categorical_forest_interaction/01_publication_figures.py
python  06_categorical_forest_interaction/02_publication_figures_ids06_IplusT.py

# Stage 3: Target Trial Emulation — IT_RULES_v2 drives both cohorts
Rscript 07_target_trial_emulation/tte_IT_R_two_cohorts.R  data/
python  07_target_trial_emulation/tte_IT_R_figures_two_cohorts.py
```

---

## Statistical methods reference

| Method                          | Implementation                                |
| ------------------------------- | --------------------------------------------- |
| CBPS-IPTW                       | `WeightIt` (`method = "cbps"`) / `glmnet`     |
| 1:1 nearest-neighbor PSM        | `MatchIt` (`method = "nearest"`, `caliper`)   |
| Overlap weighting               | `PSweight` / hand-coded weights               |
| Cox regression (weighted)       | `survival::coxph` + robust sandwich SE        |
| RMST                            | Weighted KM step integration + bootstrap CI   |
| RCS × treatment interaction     | `rms::cph` with `rcs(x, nk = 3)`              |
| Clone-Censor-Weight (CCW)       | Custom implementation                         |
| Stabilized IPCW                 | Person-period pooled logistic regression      |
| Bootstrap CIs                   | 500 (TTE) / 200 (RCS) / 300 (forest)          |
| E-value                         | Sensitivity to unmeasured confounding         |

---

## Data availability

> Patient-level data are **not** included in this repository due to ethical/privacy
> restrictions. De-identified data may be available from the corresponding author
> on reasonable request and approval by the institutional review board.
> The repository contains analytic code only.

### Required input files (placed in `data/`, not tracked by git)

| File                                       | Purpose                                       |
| ------------------------------------------ | --------------------------------------------- |
| `HAIC_NO_TACE_4_TIDY_baseline.csv`         | Baseline characteristics (no TACE patients)   |
| `HAIC_NO_TACE_4_TIDY_baseline_imputed.csv` | Imputed baseline                              |
| `HAIC_NO_TACE_4_TIDY_longitudinal.csv`     | Per-cycle longitudinal labs                   |
| `analysis_ready.csv`                       | 7-group classified analytic dataset           |
| `matched_ids_*.csv`                        | Pre-computed PSM matched IDs (6 pairs)        |

---

## Citation

If you use any portion of this code, please cite the paper (see `CITATION.cff`).

## License

Released under the MIT License — see `LICENSE` file.

## Contact

Corresponding author: **Sun Yat-sen University Cancer Center, Department of Hepatobiliary Oncology**.
For questions about the code, open an issue in this repository.
