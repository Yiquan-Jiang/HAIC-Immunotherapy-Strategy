# 00 — Data preparation

Build the analysis-ready dataset by classifying patients into seven mutually
exclusive treatment groups based on whether and when systemic therapy
(immunotherapy and/or antiangiogenic therapy) was added to induction HAIC.

## Scripts

| Script | Description |
|---|---|
| `step0_prepare_data.py` | Reads raw baseline + longitudinal CSVs, classifies patients into the 7 groups, writes `analysis_ready.csv` |

## Inputs (expected in `data/`, not tracked)

- `HAIC_NO_TACE_4_TIDY_baseline.csv`
- `HAIC_NO_TACE_4_TIDY_longitudinal.csv`

## Outputs

- `analysis_ready.csv` — 7-group classified analytic dataset (~3,885 patients)

## Run

```bash
python step0_prepare_data.py
```
