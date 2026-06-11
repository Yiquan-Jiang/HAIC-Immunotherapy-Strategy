# 00 — Data preparation

Build the analysis-ready dataset by classifying patients into eight mutually
exclusive treatment groups. Groups 1–7 are the induction-HAIC strategies,
defined by whether and when systemic therapy (immunotherapy and/or
antiangiogenic therapy) was added to HAIC. Group 8 is a systemic-therapy-only
comparator that received **no** induction HAIC, added as an external benchmark
for the Stage-1 overall-survival analysis.

## Scripts

| Script | Description |
|---|---|
| `step0_prepare_data.py` | Reads raw baseline + longitudinal CSVs, classifies patients into the 7 HAIC groups, writes `analysis_ready.csv` |
| `00b_prepare_systemic_it_group.py` | Builds the systemic-only 8th group (`Systemic_I+T`) from the systemic-therapy source, imputing missing covariates with the same MICE recipe (fit on the new group alone, seed 42), and appends it to produce `analysis_ready_8group.csv` |

## Inputs (expected in `data/`, not tracked)

- `HAIC_NO_TACE_4_TIDY_baseline.csv`
- `HAIC_NO_TACE_4_TIDY_longitudinal.csv`
- systemic-therapy source workbook (for `00b`, not tracked)

## Outputs

- `analysis_ready.csv` — 7-group HAIC analytic dataset (~3,885 patients)
- `analysis_ready_8group.csv` — adds the systemic-only 8th group (N = 570 → pooled 4,455 patients)

## Run

```bash
python step0_prepare_data.py              # 7 HAIC groups
python 00b_prepare_systemic_it_group.py   # + systemic-only 8th group
```
