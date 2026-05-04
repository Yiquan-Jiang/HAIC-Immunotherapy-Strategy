# utils_runners — End-to-end pipeline drivers

One-shot bash scripts that chain Stage 0 → Stage N for a specific analysis.

| Runner | What it does | Paper output |
|---|---|---|
| `run_tte_two_cohorts.sh`         | Stage 3 — TTE pipeline (`tte_IT_R_two_cohorts.R` core + `tte_IT_R_figures_two_cohorts.py` figures). Drives **both** Cohort A (cohort_3matched, IT add-on) and Cohort B (cohort_7group_psm02, ICI add-on) in one invocation. | **Fig 7** |
| `run_all_group7.sh`              | Stage 0 → Stage 1 → Stage 2 → swimmer plot end-to-end | Fig 2-5 + Tables 1-5 |
| `run_rcs_interaction_group7.sh`  | Stage 2 RCS interaction (per-pair driver) — see also the canonical FINAL_SCRIPTS.md flow in `05_rcs_interaction/` | **Fig 6** |

Make sure your `data/` directory is populated and R / Python environments
are configured (see root `README.md`) before running these scripts.

```bash
# example
bash utils_runners/run_tte_two_cohorts.sh
```
