#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Publication-quality forest plot: each treatment group vs HAIC+I+T (concurrent)
HR = h(row group) / h(reference).  Data from PSM survival_analysis_final.csv.
Style: Nature / Cell — table-style layout with aligned columns.
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import numpy as np
import os
import pandas as pd

BASE_DIR = (
    "/Users/yqj/Nutstore Files/我的坚果云/Liver_tumor_big_data/"
    "FIRST_LINE_HAIC_2025-12-30/HAIC_NO_TACE_4_TIDY/update_group_7"
)
RES_DIR = os.path.join(BASE_DIR, "results", "psm_balance_tables_complete")
OUTPUT_DIR = os.path.join(BASE_DIR, "figures", "psm_pub_quality")
os.makedirs(OUTPUT_DIR, exist_ok=True)

REF_GROUP = "HAIC+I+T_concurrent"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8,
    "axes.labelsize": 9,
    "axes.titlesize": 10,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "legend.fontsize": 7,
    "axes.linewidth": 0.6,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 3.0,
    "ytick.major.size": 3.0,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.05,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "axes.grid": False,
})

# ── NPG palette ──
COL_NS = "#3C5488"
COL_SIG = "#E64B35"
COL_REF = "#00A087"
COL_GRID = "#EDEDED"
COL_RULE = "#333333"
COL_TXT = "#222222"
COL_TXT2 = "#555555"

GROUP_ORDER = [
    "HAIC_alone",
    "HAIC+I_concurrent",
    "HAIC_then_I",
    "HAIC+T_concurrent",
    "HAIC_then_T",
    "HAIC+I+T_concurrent",
    "HAIC_then_I+T",
]

GROUP_LABELS = {
    "HAIC_alone":            "HAIC alone",
    "HAIC+I_concurrent":     "HAIC + Immunotherapy (conc.)",
    "HAIC_then_I":           "HAIC → Immunotherapy (seq.)",
    "HAIC+T_concurrent":     "HAIC + Targeted therapy (conc.)",
    "HAIC_then_T":           "HAIC → Targeted therapy (seq.)",
    "HAIC+I+T_concurrent":   "HAIC + I + T (concurrent)",
    "HAIC_then_I+T":         "HAIC → I + T (sequential)",
}


def invert_hr_ci(hr, lo, hi):
    return 1.0 / hr, 1.0 / hi, 1.0 / lo


def row_vs_reference(df, other):
    m1 = (df["Group1"] == other) & (df["Group2"] == REF_GROUP)
    m2 = (df["Group1"] == REF_GROUP) & (df["Group2"] == other)
    if m1.any():
        r = df.loc[m1].iloc[0]
        hr, lo, hi = float(r["HR"]), float(r["CI_lower"]), float(r["CI_upper"])
        return invert_hr_ci(hr, lo, hi), int(r["N1_after"]), int(r["N2_after"]), float(r["P_value"])
    if m2.any():
        r = df.loc[m2].iloc[0]
        return (
            (float(r["HR"]), float(r["CI_lower"]), float(r["CI_upper"])),
            int(r["N2_after"]), int(r["N1_after"]), float(r["P_value"]),
        )
    raise ValueError(f"未找到 {other} vs {REF_GROUP}")


def main():
    surv = pd.read_csv(os.path.join(RES_DIR, "survival_analysis_final.csv"))

    others = [g for g in GROUP_ORDER if g != REF_GROUP]
    rows = []
    for g in others:
        (hr, lo, hi), n_other, n_ref, p = row_vs_reference(surv, g)
        rows.append(dict(
            group=g, label=GROUP_LABELS[g],
            HR=hr, CI_lower=lo, CI_upper=hi, P_value=p,
            N_other=n_other, N_ref=n_ref, is_ref=False,
        ))

    ref_idx = 5
    rows.insert(ref_idx, dict(
        group=REF_GROUP, label=GROUP_LABELS[REF_GROUP],
        HR=1.0, CI_lower=np.nan, CI_upper=np.nan, P_value=np.nan,
        N_other=np.nan, N_ref=np.nan, is_ref=True,
    ))

    pd.DataFrame(rows).to_csv(
        os.path.join(OUTPUT_DIR, "forest_vs_IT_concurrent_psm_after.csv"), index=False)

    n = len(rows)

    # ── Figure geometry ──
    fig_w = 7.2
    fig_h = 0.42 * n + 2.0
    fig = plt.figure(figsize=(fig_w, fig_h))

    left_frac = 0.38
    mid_frac = 0.30
    right_frac = 0.32

    ax_forest = fig.add_axes([left_frac, 0.16, mid_frac, 0.74])
    ax_left = fig.add_axes([0.0, 0.16, left_frac, 0.74])
    ax_right = fig.add_axes([left_frac + mid_frac, 0.16, right_frac, 0.74])

    for a in (ax_left, ax_right):
        a.set_xlim(0, 1)
        a.set_ylim(-0.5, n - 0.5)
        a.axis("off")

    ax_forest.set_ylim(-0.5, n - 0.5)
    X_LO, X_HI = 0.3, 3.5
    ax_forest.set_xlim(X_LO, X_HI)
    ax_forest.set_xscale("log")
    ax_forest.spines["left"].set_visible(False)
    ax_forest.spines["top"].set_visible(False)
    ax_forest.spines["right"].set_visible(False)
    ax_forest.set_yticks([])

    ax_forest.axvline(x=1.0, color=COL_RULE, linestyle="--", linewidth=0.7,
                      alpha=0.6, zorder=1)

    for i in range(n):
        y = n - 1 - i
        if i % 2 == 0:
            for a in (ax_left, ax_forest, ax_right):
                a.axhspan(y - 0.5, y + 0.5, color=COL_GRID, zorder=0, linewidth=0)

    # ── Header ──
    hdr_y = n - 0.5 + 0.35
    for a in (ax_left, ax_forest, ax_right):
        a.set_ylim(-0.5, n - 0.5 + 0.7)

    ax_left.text(0.04, hdr_y, "Treatment group", ha="left", va="center",
                 fontsize=8, fontweight="bold", color=COL_TXT)
    ax_left.text(0.88, hdr_y, "n", ha="center", va="center",
                 fontsize=8, fontweight="bold", color=COL_TXT)

    ax_right.text(0.04, hdr_y, "HR (95% CI)", ha="left", va="center",
                  fontsize=8, fontweight="bold", color=COL_TXT)
    ax_right.text(0.82, hdr_y, "P value", ha="center", va="center",
                  fontsize=8, fontweight="bold", color=COL_TXT)

    rule_kw = dict(color=COL_RULE, linewidth=0.8, clip_on=False)
    for a in (ax_left, ax_forest, ax_right):
        a.axhline(y=n - 0.5 + 0.08, **rule_kw)
        a.axhline(y=-0.5, **rule_kw)

    # ── Rows ──
    for i, rp in enumerate(rows):
        y = n - 1 - i
        is_ref = rp["is_ref"]
        hr = rp["HR"]
        lo, hi = rp["CI_lower"], rp["CI_upper"]
        p = rp["P_value"]

        label_style = dict(fontsize=7.5, va="center", color=COL_TXT)
        if is_ref:
            label_style["fontweight"] = "bold"
            label_style["color"] = COL_REF

        ax_left.text(0.04, y, rp["label"], ha="left", **label_style)

        if is_ref:
            ax_left.text(0.88, y, "—", ha="center", va="center",
                         fontsize=7.5, color=COL_REF, fontweight="bold")
            ax_right.text(0.04, y, "1.00 (reference)", ha="left", va="center",
                          fontsize=7.5, color=COL_REF, fontweight="bold")
            ax_right.text(0.82, y, "—", ha="center", va="center",
                          fontsize=7.5, color=COL_REF, fontweight="bold")
            ax_forest.scatter([1.0], [y], marker="D", s=55,
                              color=COL_REF, zorder=5, linewidths=0.4,
                              edgecolors="white")
            continue

        sig = p < 0.05
        col = COL_SIG if sig else COL_NS

        n_str = f"{rp['N_other']}"
        ax_left.text(0.88, y, n_str, ha="center", va="center",
                     fontsize=7.5, color=COL_TXT)

        ci_lo_c = max(lo, X_LO * 1.02)
        ci_hi_c = min(hi, X_HI * 0.98)

        ax_forest.plot([ci_lo_c, ci_hi_c], [y, y],
                       color=col, linewidth=1.6, solid_capstyle="round", zorder=3)
        for xc in (ci_lo_c, ci_hi_c):
            ax_forest.plot([xc, xc], [y - 0.12, y + 0.12],
                           color=col, linewidth=0.9, zorder=3)

        weight = min(max(rp["N_other"] / 500, 0.5), 1.8)
        diamond_s = 50 * weight
        ax_forest.scatter([hr], [y], marker="D", s=diamond_s,
                          color=col, zorder=5, linewidths=0.4,
                          edgecolors="white")

        hr_str = f"{hr:.2f} ({lo:.2f}\u2013{hi:.2f})"
        p_str = "< 0.001" if p < 0.001 else f"{p:.3f}"
        txt_col = "#B03000" if sig else COL_TXT
        ax_right.text(0.04, y, hr_str, ha="left", va="center",
                      fontsize=7.5, color=txt_col,
                      fontweight="bold" if sig else "normal")
        ax_right.text(0.82, y, p_str, ha="center", va="center",
                      fontsize=7.5, color=txt_col,
                      fontweight="bold" if sig else "normal")

    # ── X-axis ──
    xticks = [0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0]
    ax_forest.set_xticks(xticks)
    ax_forest.set_xticklabels([str(x) for x in xticks], fontsize=7)
    ax_forest.set_xlabel("Hazard ratio", fontsize=8.5, labelpad=14)

    ax_forest.annotate(
        "", xy=(X_LO, -0.06), xytext=(0.85, -0.06),
        xycoords=("data", "axes fraction"), textcoords=("data", "axes fraction"),
        arrowprops=dict(arrowstyle="-|>", color=COL_TXT2, lw=0.6),
    )
    ax_forest.annotate(
        "", xy=(X_HI, -0.06), xytext=(1.15, -0.06),
        xycoords=("data", "axes fraction"), textcoords=("data", "axes fraction"),
        arrowprops=dict(arrowstyle="-|>", color=COL_TXT2, lw=0.6),
    )
    ax_forest.text(0.18, -0.12, "Favours\nrow group",
                   ha="center", va="top", fontsize=6, color=COL_TXT2,
                   transform=ax_forest.transAxes, linespacing=0.9)
    ax_forest.text(0.82, -0.12, "Favours\nreference",
                   ha="center", va="top", fontsize=6, color=COL_TXT2,
                   transform=ax_forest.transAxes, linespacing=0.9)

    # ── Legend (bottom) ──
    leg_handles = [
        mlines.Line2D([], [], color=COL_NS, marker="D", markersize=5,
                      markerfacecolor=COL_NS, markeredgecolor="white",
                      markeredgewidth=0.4, linewidth=1.4, label="P \u2265 0.05"),
        mlines.Line2D([], [], color=COL_SIG, marker="D", markersize=5,
                      markerfacecolor=COL_SIG, markeredgecolor="white",
                      markeredgewidth=0.4, linewidth=1.4, label="P < 0.05"),
        mlines.Line2D([], [], color=COL_REF, marker="D", markersize=5,
                      markerfacecolor=COL_REF, markeredgecolor="white",
                      markeredgewidth=0.4, linewidth=0, label="Reference"),
    ]
    ax_forest.legend(handles=leg_handles, loc="upper center",
                     bbox_to_anchor=(0.5, -0.22), ncol=3,
                     fontsize=7, frameon=False, handlelength=1.8,
                     columnspacing=1.5)

    # ── Title ──
    fig.text(0.02, 0.96,
             "Overall survival: hazard ratio vs HAIC + I + T (concurrent)",
             fontsize=10, fontweight="bold", color=COL_TXT, va="top")
    fig.text(0.02, 0.925,
             "Propensity-score matched cohorts  |  1 : 1 nearest-neighbour matching",
             fontsize=7.5, color=COL_TXT2, va="top")

    base = os.path.join(OUTPUT_DIR, "HR_forest_vs_IT_concurrent_psm_after")
    fig.savefig(f"{base}.pdf", bbox_inches="tight", pad_inches=0.08)
    fig.savefig(f"{base}.png", dpi=300, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    print(f"已保存: {base}.pdf / .png")


if __name__ == "__main__":
    main()
