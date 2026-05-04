# utils_runners — End-to-end pipeline drivers

One-shot bash scripts that chain Stage 0 → Stage N for a specific analysis.

| Runner | What it does | Paper output |
|---|---|---|
| `run_tte_cohort_7group_psm02.sh`     | Stage 3 — TTE pipeline using **NLR-based** decision rule (Cohort A) | **Fig 7A** |
| `run_tte_cohort_7group_psm02_piv.sh` | Stage 3 — TTE pipeline using **PIV-based** decision rule (Cohort B) | **Fig 7B** |
| `run_all_group7.sh`                  | Group-7 project end-to-end (Stage 0 → Stage 1 → Stage 2 → swimmer plot) | Fig 2-5 + Tables 1-5 |
| `run_rcs_interaction_group7.sh`      | Stage 2 RCS interaction (per-pair version) | **Fig 6** |

> Make sure your `data/` directory is populated and R / Python environments
> are configured (see root `README.md`) before running these scripts.

```bash
# example
bash utils_runners/run_tte_cohort_7group_psm02.sh
bash utils_runners/run_tte_cohort_7group_psm02_piv.sh
```
