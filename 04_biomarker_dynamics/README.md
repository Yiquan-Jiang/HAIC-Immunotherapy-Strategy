# 04 — Longitudinal tumor-biomarker dynamics

Track AFP, PIVKA-II, ALBI score, NLR, SII and PIV across pre-HAIC and the
6-month, 12-month, and 24-month landmarks **within each PSM-matched pair**,
to provide a biological correlate for the equivalent OS observed across
HAIC + systemic-therapy strategies.

> **Paper output:** Results §"Longitudinal tumor biomarker dynamics" + supporting supplementary panels.

## Scripts

| Script | Description |
|---|---|
| `psm_afp_pivka_dynamics.py`            | AFP / PIVKA-II / inflammatory-index trajectories within the 21 PSM pairs |
| `plot_haic_then_i_to_target_interval.R`| Distribution of HAIC → ICI initiation interval, target-window adherence |

## Outputs

- Median trajectory plots (AFP, PIVKA-II, ALBI, NLR, SII, PIV) per PSM pair
- "Time to systemic-therapy add-on" histograms used to motivate the dynamic-decision triggers in Stage 3
