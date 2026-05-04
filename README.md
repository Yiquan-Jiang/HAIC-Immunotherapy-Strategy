# HAIC-Immunotherapy-Strategy-JCO

**Code repository for the manuscript:**

> *Immunotherapy Following Induction Arterial Chemotherapy for Unresectable Hepatocellular Carcinoma: A Biomarker-Guided On-Demand Decision Framework*

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
         → Tables 1-5, Fig 2, Fig 3, Fig 4
                       │
                       ▼
Stage 2  Biomarker × treatment interaction
         RCS continuous interaction + categorical forest plots
         → Fig 5 (swimmer), Fig 6 (RCS matrix)
                       │
                       ▼
Stage 3  Target Trial Emulation (CCW + IPCW)
         Dynamic biomarker-guided strategy vs. early combination
         → Fig 7 (two parallel cohorts: ICI-only & ICI + Anti-angio)
```

---

## Repository layout — mapped to the paper

```
HAIC-Immunotherapy-Strategy-JCO/
├── 00_data_preparation/                   ← Build analysis-ready dataset (7 groups)
│
├── 01_psm_iptw_overall_survival/          ← Stage 1
│   │   Output: Tables 1-5, Fig 2, Fig 3
│   ├── step3_psm_analysis.R                  one-to-one nearest-neighbor PSM
│   ├── step3b_psm_vs_template.R              standardized PSM vs template
│   ├── step4_km_curves.py                    KM curves overall (Fig 2A unweighted)
│   ├── step4b_km_template_matched.py         KM curves of PSM-matched pairs (Fig 3)
│   ├── step5_forest_plot.py                  HR/RMST forest summary
│   ├── step5b_forest_vs_IT_concurrent.py     vs HAIC + ICI concurrent
│   ├── step5c_forest_vs_HAIC_alone.py        vs HAIC alone (21 pairwise)
│   └── step6_tables_and_loveplots.R          Tables 1-5 + love plots (balance)
│
├── 02_subgroup_overlap_weighting/         ← Fig 4 high-risk subgroups
│   │   Concurrent vs delayed-sequential under overlap weighting
│   ├── step7_ow_balance_table.R              OW balance tables
│   ├── step7_ow_balance_figure.py            OW love-plot
│   ├── step7_subgroup_ow.R                   Cox HR per subgroup (OW-weighted)
│   ├── step7_subgroup_analysis.py            ΔRMST + interaction tests
│   └── step7_subgroup_plots.py               Forest panels for Fig 4
│
├── 03_swimmer_plot/                       ← Fig 5
│   └── swimmer_plot_7groups.R                Per-patient timelines for 7 groups
│
├── 04_biomarker_dynamics/                 ← Longitudinal AFP/PIVKA-II/inflammatory indices
│   ├── psm_afp_pivka_dynamics.py             Trajectories within PSM-matched pairs
│   └── plot_haic_then_i_to_target_interval.R Time-to-trigger interval visualization
│
├── 05_rcs_interaction/                    ← Stage 2 (continuous) — Fig 6
│   ├── tte_project_version/                  Comprehensive 6-pair RCS × treatment Cox
│   │   ├── build_all_pairs_cohorts.py        Build 6 IPTW pairwise cohorts
│   │   ├── rcs_all_pairs_dual_timescale.R    Main analysis (IPTW)
│   │   ├── build_cohort_psm.py               Build 6 PSM cohorts
│   │   ├── RCS_PSM_dual_timescale.R          Sensitivity analysis (PSM)
│   │   ├── RCS_PSM_matrix_panel.R            Assemble Fig 6 matrix panel
│   │   └── extract_pre_it_for_psm.py         Extract pre-IT lab values
│   ├── afp_pivka_composite/                  AFP-PIVKA composite biomarker analysis
│   │   ├── 00a_extract_pre_it_labs.py
│   │   ├── 01_rcs_afp_pivka_composite.R
│   │   └── 02_rcs_matrix_panel.R
│   └── publication_figures/                  Final figure assembly (PSM + IPTW)
│       ├── make_publication_figures_iptw.R
│       └── make_publication_figures_psm.R
│
├── 06_categorical_forest_interaction/     ← Stage 2 (categorical) — Fig 5 forest
│   └── final_publication/                    Publication-ready interaction forest plots
│       ├── 01_publication_figures.py             PSM02: HAIC alone vs HAIC → ICI
│       ├── 02_publication_figures_ids06_IplusT.py PSM06: vs HAIC → ICI + Anti-angio
│       ├── 03_publication_figures_iptw_psm02.py   IPTW + PSM02 combined panel
│       ├── 03_publication_figures_ids05_IplusT_concurrent.py
│       ├── 04_publication_figures_iptw_psm06_IplusT.py
│       └── 05_publication_figures_iptw_psm05_IplusT_concurrent.py
│
├── 07_target_trial_emulation/             ← Stage 3 — Fig 7
│   │   Clone-Censor-Weight (CCW) + stabilized IPCW + weighted Cox
│   │   Two parallel cohorts:
│   │     - Cohort A (NLR-based): ICI-only emulation
│   │     - Cohort B (PIV-based): ICI + Anti-angio emulation
│   ├── tte_piv_R_core_cohort_7group_psm02.R   TTE core (PIV rule, v3, 7-group PSM02)
│   ├── tte_nlr_R_figures.py                   Publication figures (works for NLR & PIV)
│   ├── tte_pathway_visualization.R            Strategy pathway diagram
│   ├── tte_pathway_visualization_alt_samples.R
│   └── generate_tte_flow_drawio_two_cohorts.R CONSORT-style flow (drawio XML)
│
├── 08_schematic_figures/                  ← Conceptual figures
│   ├── plot_tte_schematic_v3.py              TTE conceptual schematic (Nature-style)
│   ├── plot_tte_nlr_strategy_flowchart.py    NLR-based decision-rule flowchart
│   └── plot_tte_IT_rules_v2_schematic.py     ICI-rule schematic (v2)
│
├── utils_runners/                         ← One-shot runners (.sh)
│   ├── run_tte_cohort_7group_psm02_piv.sh    PIV-rule TTE pipeline
│   ├── run_all_group7.sh                     Group-7 project end-to-end
│   └── run_rcs_interaction_group7.sh         RCS interaction (per-pair) driver
│
└── docs/                                  ← Original project READMEs (preserved)
    ├── README_TTE_project_original.md
    ├── README_group7_project_original.md
    └── ORGANIZATION_REPORT_original.md
```

---

## Script ↔ Figure / Table mapping

| Paper element | Folder | Primary scripts |
|---|---|---|
| **Table 1** Baseline characteristics — overall cohort | `01_psm_iptw_overall_survival/` | `step6_tables_and_loveplots.R` |
| **Table 2** Baseline before PSM (7 groups)            | `01_psm_iptw_overall_survival/` | `step6_tables_and_loveplots.R` |
| **Tables 3-5** PSM balance for selected pairs         | `01_psm_iptw_overall_survival/` | `step3_psm_analysis.R` + `step6_tables_and_loveplots.R` |
| **Fig 1** Study flow                                  | `08_schematic_figures/` | (consort-style schematic) |
| **Fig 2** OS across 7 strategies (KM + IPTW)          | `01_psm_iptw_overall_survival/` | `step4_km_curves.py`, `step5_forest_plot.py` |
| **Fig 3** PSM head-to-head OS                          | `01_psm_iptw_overall_survival/` | `step4b_km_template_matched.py`, `step5b_forest_vs_IT_concurrent.py`, `step5c_forest_vs_HAIC_alone.py` |
| **Fig 4** Concurrent vs delayed-sequential (subgroups, OW) | `02_subgroup_overlap_weighting/` | `step7_subgroup_*` + `step7_ow_balance_*` |
| **Fig 5** Swimmer plots (7 groups)                    | `03_swimmer_plot/` | `swimmer_plot_7groups.R` |
| **Fig 6** RCS interaction matrix                      | `05_rcs_interaction/` | `rcs_all_pairs_dual_timescale.R`, `RCS_PSM_matrix_panel.R`, `make_publication_figures_*.R` |
| **Fig 6 (forest panels)** Categorical interaction     | `06_categorical_forest_interaction/` | `01_publication_figures.py` (PSM02), `02_publication_figures_ids06_IplusT.py` (PSM06) |
| **Fig 7** Target Trial Emulation — two cohorts        | `07_target_trial_emulation/` | `tte_piv_R_core_cohort_7group_psm02.R` (TTE core, v3), `tte_nlr_R_figures.py` (publication figures) |
| Longitudinal biomarker dynamics (Results §Longitudinal tumor biomarker dynamics) | `04_biomarker_dynamics/` | `psm_afp_pivka_dynamics.py` |

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

### Pipeline (paper-order)

```bash
# Stage 0: prepare analysis-ready dataset (requires raw inputs)
python  00_data_preparation/step0_prepare_data.py

# Stage 1: PSM / IPTW + KM + tables
Rscript 01_psm_iptw_overall_survival/step3_psm_analysis.R
python  01_psm_iptw_overall_survival/step4_km_curves.py
python  01_psm_iptw_overall_survival/step4b_km_template_matched.py
python  01_psm_iptw_overall_survival/step5_forest_plot.py
Rscript 01_psm_iptw_overall_survival/step6_tables_and_loveplots.R

# Stage 1b: subgroup OW (Fig 4)
Rscript 02_subgroup_overlap_weighting/step7_ow_balance_table.R
Rscript 02_subgroup_overlap_weighting/step7_subgroup_ow.R
python  02_subgroup_overlap_weighting/step7_subgroup_analysis.py
python  02_subgroup_overlap_weighting/step7_subgroup_plots.py

# Fig 5: swimmer
Rscript 03_swimmer_plot/swimmer_plot_7groups.R

# Stage 2a: RCS interaction (continuous biomarkers)
python  05_rcs_interaction/tte_project_version/build_all_pairs_cohorts.py
Rscript 05_rcs_interaction/tte_project_version/rcs_all_pairs_dual_timescale.R
Rscript 05_rcs_interaction/publication_figures/make_publication_figures_iptw.R
# Or per-pair:
bash    utils_runners/run_rcs_interaction_group7.sh

# Stage 2b: categorical interaction forest (Fig 6 panels)
python  06_categorical_forest_interaction/final_publication/01_publication_figures.py
python  06_categorical_forest_interaction/final_publication/02_publication_figures_ids06_IplusT.py

# Stage 3: Target Trial Emulation (Fig 7)
bash    utils_runners/run_tte_cohort_7group_psm02_piv.sh    # PIV rule (Cohort B: ICI + Anti-angio)
# Manual: TTE core also drives the NLR-rule analysis (Cohort A) when invoked
#        with the NLR cohort IDs as the second argument:
#   Rscript 07_target_trial_emulation/tte_piv_R_core_cohort_7group_psm02.R \
#           data/  data/cohort_ids_nlr.csv
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
