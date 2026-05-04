# 03 — Swimmer plot of individual patient trajectories

Per-patient timelines for up to 50 randomly sampled patients per group,
showing HAIC cycles, immunotherapy / antiangiogenic initiation, response,
TACE / surgery / radiofrequency ablation events, and death / last follow-up.

> **Paper output:** Fig 5.

## Scripts

| Script | Description |
|---|---|
| `swimmer_plot_7groups.R`              | Swimmer plot for all seven treatment groups (50 patients each, sampled with fixed seed) |
| `plot_haic_then_i_to_target_interval.R` | Distribution of the HAIC → ICI initiation interval and target-window adherence (companion timing visualization) |

## Outputs

| Output | Paper element |
|---|---|
| `swimmer_plot_7groups.{pdf,png}` | **Fig 5** |
| `haic_then_i_to_target_interval.{pdf,png}` | Supplementary timing panel |
| Per-patient event-time table     | Supplements |
