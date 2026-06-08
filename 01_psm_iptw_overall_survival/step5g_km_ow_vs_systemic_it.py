#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""step5g — ATO overlap-weighted KM curves: each HAIC strategy vs Systemic I+T.

Plots the SAME balanced (ATO) basis as the step5f forest, so the survival curves match the
overlap-weighted hazard ratios. Reads ow_km_data.csv + ow_forest_data.csv (from step5e).
One panel per contrast (HAIC group vs Systemic I+T), HR(95% CI) and post-OW max|SMD| annotated.
A single 8-arm weighting cannot balance the no-HAIC arm against all 7 others at once, so the
balanced comparison is necessarily pairwise.
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

BASE_DIR = "/Users/yqj/Nutstore Files/我的坚果云/Liver_tumor_big_data/FIRST_LINE_HAIC_2025-12-30/HAIC_NO_TACE_4_TIDY/update_group_7"
RES_DIR = os.path.join(BASE_DIR, "results", "ow_vs_systemic_it_8group")
OUT_DIR = os.path.join(BASE_DIR, "figures", "ow_vs_systemic_it_8group")
os.makedirs(OUT_DIR, exist_ok=True)

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 9, "axes.labelsize": 9, "axes.titlesize": 9.5,
    "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8,
    "axes.linewidth": 0.7, "axes.spines.top": False, "axes.spines.right": False,
    "savefig.dpi": 300, "savefig.bbox": "tight", "pdf.fonttype": 42, "ps.fonttype": 42,
})

LABELS = {
    "HAIC_alone": "HAIC alone", "HAIC+I_concurrent": "HAIC + I (conc.)",
    "HAIC_then_I": "HAIC → I (seq.)", "HAIC+T_concurrent": "HAIC + T (conc.)",
    "HAIC_then_T": "HAIC → T (seq.)", "HAIC+I+T_concurrent": "HAIC + I + T (conc.)",
    "HAIC_then_I+T": "HAIC → I + T (seq.)",
}
ORDER = list(LABELS.keys())
COL_HAIC = "#D55E00"
COL_SYS = "#009E73"

km = pd.read_csv(os.path.join(RES_DIR, "ow_km_data.csv"))
fd = pd.read_csv(os.path.join(RES_DIR, "ow_forest_data.csv")).set_index("group")

fig, axes = plt.subplots(2, 4, figsize=(15, 7.2), sharex=True, sharey=True)
axes = axes.ravel()

for i, g in enumerate(ORDER):
    ax = axes[i]
    sub = km[km["group"] == g]
    for arm, col in [(g, COL_HAIC), ("Systemic_I+T", COL_SYS)]:
        a = sub[sub["arm"] == arm].sort_values("time")
        if a.empty:
            continue
        t = np.concatenate([[0], a["time"].values])
        s = np.concatenate([[1], a["surv"].values])
        lbl = LABELS[g] if arm == g else "Systemic I+T"
        ax.step(t, s, where="post", color=col, lw=1.7, label=lbl)
    r = fd.loc[g]
    p = "<0.001" if r["p"] < 0.001 else f"{r['p']:.3f}"
    ax.text(0.96, 0.96, f"HR {r['HR']:.2f} ({r['CI_lower']:.2f}–{r['CI_upper']:.2f})\n"
                        f"p {p}   max|SMD| {r['max_smd_adj']:.2f}",
            transform=ax.transAxes, ha="right", va="top", fontsize=7.5,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#cccccc", lw=0.5))
    ax.set_title(f"{LABELS[g]}  vs  Systemic I+T", fontsize=9)
    ax.set_xlim(0, 60)
    ax.set_ylim(0, 1.0)
    ax.legend(loc="lower left", frameon=False, fontsize=7.5)
    ax.grid(alpha=0.25, lw=0.4)
    if i % 4 == 0:
        ax.set_ylabel("Overall survival")
    if i >= 4:
        ax.set_xlabel("Months")

axes[3].set_xlabel("Months")
axes[-1].axis("off")
fig.suptitle("ATO overlap-weighted KM: adding HAIC vs systemic therapy alone "
             "(balanced, overlap population)", y=1.01, fontsize=11)
fig.tight_layout()
base = os.path.join(OUT_DIR, "KM_ow_vs_systemic_it_panels")
fig.savefig(f"{base}.pdf", bbox_inches="tight", pad_inches=0.08)
fig.savefig(f"{base}.png", dpi=300, bbox_inches="tight", pad_inches=0.08)
print("Saved:", base + ".pdf/.png")
