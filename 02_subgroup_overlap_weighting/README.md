# 02 — Subgroup analysis with Overlap Weighting

Concurrent vs delayed-sequential systemic therapy in **high-risk subgroups**,
re-weighted by overlap weights to handle small-cell extreme propensity scores
in selected pairwise comparisons.

> **Paper output:** Fig 4.

## Pairwise comparisons addressed

1. `HAIC + ICI concurrent` vs `HAIC → ICI`
2. `HAIC + ICI + Anti-angio concurrent` vs `HAIC → ICI + Anti-angio`

## Subgroups examined

Vascular invasion (Vp3/4 PVTT, HVTT, IVC-RA), distant metastasis,
ALBI grade, baseline AFP / PIVKA-II tertiles, baseline NLR, baseline tumor size, etc.

## Workflow

```
matched_ids_*.csv  +  analysis_ready.csv
       │
       ▼
step7_ow_balance_table.R     ── Overlap-weighting balance table per subgroup
step7_ow_balance_figure.py   ── Love plots after OW
       │
       ▼
step7_subgroup_ow.R          ── OW-weighted Cox HR per subgroup
step7_subgroup_analysis.py   ── ΔRMST + interaction P (LRT)
       │
       ▼
step7_subgroup_plots.py      ── Forest panels for Fig 4
```

## Outputs

| Output                                  | Paper element |
|-----------------------------------------|---------------|
| Subgroup forest panels (concurrent vs sequential) | **Fig 4** |
| Per-subgroup balance tables             | Supplements   |

## Statistical notes

- **Overlap weighting** (Li et al., *J Am Stat Assoc* 2018) gives weight `e(1-e)` per patient, dampening tails and yielding the average treatment effect in the overlap population (ATO).
- Interaction tested by **likelihood-ratio test** of treatment × subgroup-indicator term in a weighted Cox model.
- Subgroup-specific ΔRMST estimated by weighted KM integration with **bootstrap 95% CI** (500 replicates).
