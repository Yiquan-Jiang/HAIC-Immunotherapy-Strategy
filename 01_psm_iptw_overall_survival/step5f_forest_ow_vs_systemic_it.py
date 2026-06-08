#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""step5f — Forest plot: each HAIC strategy vs Systemic I+T, overlap-weighted (ATO).

Reads results/ow_vs_systemic_it_8group/ow_forest_data.csv (produced by step5e) and renders
a publication forest. HR = hazard of the HAIC group relative to Systemic I+T; HR<1 (left of
the line) => the HAIC strategy prolongs OS vs systemic therapy alone, in the overlap
(equipoise) population with measured confounders balanced (max post-OW |SMD| annotated).
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
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 9, "axes.labelsize": 10, "axes.titlesize": 11,
    "xtick.labelsize": 8.5, "ytick.labelsize": 9,
    "axes.linewidth": 0.7, "axes.spines.top": False, "axes.spines.right": False,
    "savefig.dpi": 300, "savefig.bbox": "tight", "pdf.fonttype": 42, "ps.fonttype": 42,
})
COL_SIG = "#E64B35"
COL_NS = "#3C5488"
COL_RULE = "#333333"

LABELS = {
    "HAIC_alone": "HAIC alone",
    "HAIC+I_concurrent": "HAIC + Immunotherapy (conc.)",
    "HAIC_then_I": "HAIC → Immunotherapy (seq.)",
    "HAIC+T_concurrent": "HAIC + Targeted therapy (conc.)",
    "HAIC_then_T": "HAIC → Targeted therapy (seq.)",
    "HAIC+I+T_concurrent": "HAIC + I + T (concurrent)",
    "HAIC_then_I+T": "HAIC → I + T (sequential)",
}
ORDER = list(LABELS.keys())

df = pd.read_csv(os.path.join(RES_DIR, "ow_forest_data.csv"))
df = df.set_index("group").loc[ORDER].reset_index()

n = len(df)
ypos = np.arange(n)[::-1]
fig, ax = plt.subplots(figsize=(8.2, 0.62 * n + 1.7))

for y, (_, r) in zip(ypos, df.iterrows()):
    sig = r["CI_upper"] < 1 or r["CI_lower"] > 1
    c = COL_SIG if sig else COL_NS
    ax.plot([r["CI_lower"], r["CI_upper"]], [y, y], color=c, lw=1.6, zorder=2)
    ax.plot(r["HR"], y, "s", color=c, ms=7, zorder=3)

ax.axvline(1.0, color=COL_RULE, lw=0.9, ls="--", zorder=1)
ax.set_yticks(ypos)
ax.set_yticklabels([LABELS[g] for g in df["group"]])
ax.set_xscale("log")
ax.set_xticks([0.25, 0.5, 1.0, 2.0])
ax.set_xticklabels(["0.25", "0.5", "1.0", "2.0"])
ax.set_xlim(0.2, 2.2)
ax.set_xlabel("Hazard ratio for death (HAIC strategy vs Systemic I+T), ATO-weighted")
ax.set_ylim(-0.8, n - 0.2)

# right-hand annotation: HR (95% CI), p, balance
xann = 2.45
ax.text(xann, n - 0.2, "HR (95% CI)        p        max|SMD|  ESS",
        fontsize=8, color="#222", va="bottom", fontweight="bold")
for y, (_, r) in zip(ypos, df.iterrows()):
    p = r["p"]
    ptxt = "<0.001" if p < 0.001 else f"{p:.3f}"
    ax.text(xann, y, f"{r['HR']:.2f} ({r['CI_lower']:.2f}–{r['CI_upper']:.2f})"
                     f"   {ptxt:>7}   {r['max_smd_adj']:.3f}    {int(r['ess_sys'])}",
            fontsize=8, va="center", color="#222")

ax.annotate("HAIC better", xy=(0.22, -0.7), xytext=(0.45, -0.7), fontsize=8, color=COL_SIG,
            va="center", ha="left", arrowprops=dict(arrowstyle="->", color=COL_SIG, lw=1))
ax.annotate("Systemic I+T better", xy=(2.1, -0.7), xytext=(1.25, -0.7), fontsize=8,
            color="#777", va="center", ha="left",
            arrowprops=dict(arrowstyle="->", color="#999", lw=1))
ax.set_title("Adding HAIC vs systemic therapy alone (overlap-weighted, ATO population)",
             loc="left", pad=10)

base = os.path.join(OUT_DIR, "HR_forest_ow_vs_systemic_it")
fig.savefig(f"{base}.pdf", bbox_inches="tight", pad_inches=0.06)
fig.savefig(f"{base}.png", dpi=300, bbox_inches="tight", pad_inches=0.06)
print("Saved:", base + ".pdf/.png")
print(df[["group", "HR", "CI_lower", "CI_upper", "p", "max_smd_adj", "ess_sys"]].to_string(index=False))
