# 03 — Swimmer plot of individual patient trajectories

Per-patient timelines for up to 50 randomly sampled patients per group,
showing HAIC cycles, immunotherapy / antiangiogenic initiation, response,
TACE / surgery / radiofrequency ablation events, and death / last follow-up.

## Scripts

| Script | Description |
|---|---|
| `swimmer_plot_7groups.R`              | Swimmer plot for all seven treatment groups (50 patients each, sampled with fixed seed) |
| `plot_haic_then_i_to_target_interval.R` | Distribution of the HAIC → ICI initiation interval and target-window adherence (companion timing visualization) |

## Outputs

- `swimmer_plot_7groups.{pdf,png}` — per-patient timelines for the 7 treatment groups
- `haic_then_i_to_target_interval.{pdf,png}` — HAIC → ICI interval / target-window adherence
- Per-patient event-time table
