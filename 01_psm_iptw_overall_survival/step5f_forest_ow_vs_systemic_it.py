#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Publication-quality forest plot: each HAIC strategy vs Systemic I+T (reference).
Two panels side by side: Left = Before weighting (unadjusted Cox), Right = After overlap
weighting (ATO). HR = h(HAIC strategy) / h(Systemic I+T). HR < 1 -> HAIC prolongs OS.
Same Nature/Cell table-style layout as step5c_forest_vs_HAIC_alone.py.

Each row is a SEPARATE focused pairwise overlap-weighting (ATO) model (varices dropped),
so the no-HAIC arm is balanced per contrast (max|SMD| <=0.15) instead of via the joint
8-group IPTW (which left it at |SMD| 0.335). After-OW HR/CI/p come from the R analysis
(step5e: ow_forest_data.csv); the Before column is the unadjusted Cox on the same pairwise
complete-case set (ow_weights.csv).
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
import os
import pandas as pd
from lifelines import CoxPHFitter

BASE_DIR = (
    "/Users/yqj/Nutstore Files/我的坚果云/Liver_tumor_big_data/"
    "FIRST_LINE_HAIC_2025-12-30/HAIC_NO_TACE_4_TIDY/update_group_7"
)
RES_DIR = os.path.join(BASE_DIR, "results", "ow_vs_systemic_it_8group")
OUTPUT_DIR = os.path.join(BASE_DIR, "figures", "ow_vs_systemic_it_8group")
os.makedirs(OUTPUT_DIR, exist_ok=True)

REF_GROUP = "Systemic_I+T"
GROUP_ORDER = ["HAIC_alone", "HAIC+I_concurrent", "HAIC_then_I", "HAIC+T_concurrent",
               "HAIC_then_T", "HAIC+I+T_concurrent", "HAIC_then_I+T"]
GROUP_LABELS = {
    "Systemic_I+T": "Systemic I+T",
    "HAIC_alone": "HAIC alone",
    "HAIC+I_concurrent": "HAIC + I (concurrent)",
    "HAIC_then_I": "HAIC → I (deferred)",
    "HAIC+T_concurrent": "HAIC + T (concurrent)",
    "HAIC_then_T": "HAIC → T (deferred)",
    "HAIC+I+T_concurrent": "HAIC + I + T (concurrent)",
    "HAIC_then_I+T": "HAIC → I + T (deferred)",
}

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 9, "axes.linewidth": 0.7, "savefig.dpi": 300, "savefig.bbox": "tight",
    "pdf.fonttype": 42, "ps.fonttype": 42,
})
COL_NS = "#3C5488"
COL_SIG = "#E64B35"
COL_REF = "#00A087"
COL_GRID = "#EDEDED"
COL_RULE = "#333333"
COL_TXT = "#222222"
COL_TXT2 = "#555555"

# ════════════════════════════════════════════════════════════════════
# 1. Compute rows: Before = unadjusted Cox; After = OW (ATO) from R
# ════════════════════════════════════════════════════════════════════
W = pd.read_csv(os.path.join(RES_DIR, "ow_weights.csv"))
fd = pd.read_csv(os.path.join(RES_DIR, "ow_forest_data.csv")).set_index("group")

before_rows, after_rows = [], []
for g in GROUP_ORDER:
    sub = W[W["group"] == g].copy()
    sub["E"] = sub["death"].astype(int)
    # Before: unadjusted Cox on the same pairwise complete-case set (treat=1 HAIC vs 0 ref)
    cph = CoxPHFitter().fit(sub[["os_months", "E", "treat"]],
                            duration_col="os_months", event_col="E")
    hr = float(np.exp(cph.params_["treat"]))
    lo, hi = np.exp(cph.confidence_intervals_.loc["treat"].values)
    pv = float(cph.summary.loc["treat", "p"])
    n_haic = int((sub["treat"] == 1).sum())
    before_rows.append(dict(group=g, label=GROUP_LABELS[g], HR=hr, CI_lower=lo,
                            CI_upper=hi, P_value=pv, N=n_haic, is_ref=False))
    # After: OW (ATO) from R
    r = fd.loc[g]
    after_rows.append(dict(group=g, label=GROUP_LABELS[g], HR=float(r["HR"]),
                           CI_lower=float(r["CI_lower"]), CI_upper=float(r["CI_upper"]),
                           P_value=float(r["p"]), N=int(r["n_haic"]), is_ref=False))
smd_min, smd_max = fd["max_smd_adj"].min(), fd["max_smd_adj"].max()


N_REF = int(fd["n_sys"].iloc[0])   # Systemic I+T n (constant = 570 across all contrasts)


def build_plot_rows(data_rows):
    rows = [dict(group=REF_GROUP, label=GROUP_LABELS[REF_GROUP], HR=1.0, CI_lower=np.nan,
                 CI_upper=np.nan, P_value=np.nan, N=N_REF, is_ref=True)]
    rows.extend(data_rows)
    return rows


def draw_forest_panel(fig, rect, rows, panel_title, subtitle):
    n = len(rows)
    left_frac, mid_frac, right_frac = 0.40, 0.28, 0.32
    ax_left = fig.add_axes([rect[0], rect[1], rect[2] * left_frac, rect[3]])
    ax_forest = fig.add_axes([rect[0] + rect[2] * left_frac, rect[1],
                              rect[2] * mid_frac, rect[3]])
    ax_right = fig.add_axes([rect[0] + rect[2] * (left_frac + mid_frac), rect[1],
                             rect[2] * right_frac, rect[3]])
    for a in (ax_left, ax_right):
        a.set_xlim(0, 1)
        a.set_ylim(-0.5, n - 0.5)
        a.axis("off")
    X_LO, X_HI = 0.25, 4.0
    ax_forest.set_ylim(-0.5, n - 0.5)
    ax_forest.set_xlim(X_LO, X_HI)
    ax_forest.set_xscale("log")
    for sp in ("left", "top", "right"):
        ax_forest.spines[sp].set_visible(False)
    ax_forest.set_yticks([])
    ax_forest.axvline(x=1.0, color=COL_RULE, linestyle="--", linewidth=0.7, alpha=0.6, zorder=1)
    for i in range(n):
        y = n - 1 - i
        if i % 2 == 0:
            for a in (ax_left, ax_forest, ax_right):
                a.axhspan(y - 0.5, y + 0.5, color=COL_GRID, zorder=0, linewidth=0)
    hdr_y = n - 0.5 + 0.35
    for a in (ax_left, ax_forest, ax_right):
        a.set_ylim(-0.5, n - 0.5 + 0.7)
    ax_left.text(0.04, hdr_y, "Treatment group", ha="left", va="center", fontsize=9,
                 fontweight="bold", color=COL_TXT)
    ax_left.text(0.92, hdr_y, "n", ha="center", va="center", fontsize=9,
                 fontweight="bold", color=COL_TXT)
    ax_right.text(0.04, hdr_y, "HR (95% CI)", ha="left", va="center", fontsize=9,
                  fontweight="bold", color=COL_TXT)
    ax_right.text(0.85, hdr_y, "P value", ha="center", va="center", fontsize=9,
                  fontweight="bold", color=COL_TXT)
    rule_kw = dict(color=COL_RULE, linewidth=0.8, clip_on=False)
    for a in (ax_left, ax_forest, ax_right):
        a.axhline(y=n - 0.5 + 0.08, **rule_kw)
        a.axhline(y=-0.5, **rule_kw)

    for i, rp in enumerate(rows):
        y = n - 1 - i
        is_ref = rp["is_ref"]
        hr, lo, hi, p = rp["HR"], rp["CI_lower"], rp["CI_upper"], rp["P_value"]
        label_style = dict(fontsize=8.5, va="center", color=COL_TXT)
        if is_ref:
            label_style["fontweight"] = "bold"
            label_style["color"] = COL_REF
        ax_left.text(0.04, y, rp["label"], ha="left", **label_style)
        if is_ref:
            ax_left.text(0.92, y, f"{int(rp['N'])}", ha="center", va="center", fontsize=8.5,
                         color=COL_REF, fontweight="bold")
            ax_right.text(0.04, y, "1.00 (reference)", ha="left", va="center",
                          fontsize=8.5, color=COL_REF, fontweight="bold")
            ax_right.text(0.85, y, "—", ha="center", va="center", fontsize=8.5,
                          color=COL_REF, fontweight="bold")
            ax_forest.scatter([1.0], [y], marker="D", s=65, color=COL_REF, zorder=5,
                              linewidths=0.4, edgecolors="white")
            continue
        sig = p < 0.05 if not np.isnan(p) else False
        col = COL_SIG if sig else COL_NS
        ax_left.text(0.92, y, f"{int(rp['N'])}", ha="center", va="center",
                     fontsize=8.5, color=COL_TXT)
        ci_lo_c = max(lo, X_LO * 1.02)
        ci_hi_c = min(hi, X_HI * 0.98)
        ax_forest.plot([ci_lo_c, ci_hi_c], [y, y], color=col, linewidth=1.8,
                       solid_capstyle="round", zorder=3)
        for xc in (ci_lo_c, ci_hi_c):
            ax_forest.plot([xc, xc], [y - 0.14, y + 0.14], color=col, linewidth=1.0, zorder=3)
        weight = min(max(rp["N"] / 500, 0.5), 1.8)
        ax_forest.scatter([hr], [y], marker="D", s=60 * weight, color=col, zorder=5,
                          linewidths=0.4, edgecolors="white")
        hr_str = f"{hr:.2f} ({lo:.2f}–{hi:.2f})"
        p_str = "< 0.001" if p < 0.001 else f"{p:.3f}"
        txt_col = "#B03000" if sig else COL_TXT
        ax_right.text(0.04, y, hr_str, ha="left", va="center", fontsize=8.5,
                      color=txt_col, fontweight="bold" if sig else "normal")
        ax_right.text(0.85, y, p_str, ha="center", va="center", fontsize=8.5,
                      color=txt_col, fontweight="bold" if sig else "normal")

    xticks = [0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0]
    ax_forest.set_xticks(xticks)
    ax_forest.set_xticklabels([str(x) for x in xticks], fontsize=8)
    ax_forest.set_xlabel("Hazard ratio", fontsize=9.5, labelpad=10)
    ax_forest.text(0.18, -0.10, "Favours\nHAIC strategy", ha="center", va="top",
                   fontsize=7, color=COL_TXT2, transform=ax_forest.transAxes, linespacing=0.9)
    ax_forest.text(0.82, -0.10, "Favours\nSystemic I+T", ha="center", va="top",
                   fontsize=7, color=COL_TXT2, transform=ax_forest.transAxes, linespacing=0.9)
    ax_forest.set_title(panel_title, fontsize=11, fontweight="bold", loc="center", pad=18)
    ax_forest.text(0.5, 1.06, subtitle, fontsize=7.5, color=COL_TXT2, ha="center",
                   va="bottom", transform=ax_forest.transAxes)
    return ax_forest


rows_before = build_plot_rows(before_rows)
rows_after = build_plot_rows(after_rows)

fig_w = 10.0
fig_h = 0.50 * len(rows_before) + 2.2
fig = plt.figure(figsize=(fig_w, fig_h))
panel_gap = 0.03
panel_w = (1.0 - panel_gap) / 2.0
bot, top_h = 0.14, 0.68

draw_forest_panel(fig, [0.0, bot, panel_w, top_h], rows_before,
                  "A  Before weighting", "Unadjusted Cox regression")
draw_forest_panel(fig, [panel_w + panel_gap, bot, panel_w, top_h], rows_after,
                  "B  After overlap weighting",
                  f"Pairwise ATO; max|SMD| {smd_min:.2f}–{smd_max:.2f} (robust Cox)")

leg_handles = [
    mlines.Line2D([], [], color=COL_NS, marker="D", markersize=6, markerfacecolor=COL_NS,
                  markeredgecolor="white", markeredgewidth=0.4, linewidth=1.6, label="P ≥ 0.05"),
    mlines.Line2D([], [], color=COL_SIG, marker="D", markersize=6, markerfacecolor=COL_SIG,
                  markeredgecolor="white", markeredgewidth=0.4, linewidth=1.6, label="P < 0.05"),
    mlines.Line2D([], [], color=COL_REF, marker="D", markersize=6, markerfacecolor=COL_REF,
                  markeredgecolor="white", markeredgewidth=0.4, linewidth=0,
                  label="Reference (Systemic I+T)"),
]
fig.legend(handles=leg_handles, loc="lower center", bbox_to_anchor=(0.5, 0.01), ncol=3,
           fontsize=8.5, frameon=False, handlelength=1.8, columnspacing=1.5)
fig.text(0.02, 0.985, "Overall survival: hazard ratio vs Systemic I+T (reference)",
         fontsize=12, fontweight="bold", color=COL_TXT, va="top")
fig.text(0.02, 0.955,
         "HR = h(HAIC strategy) / h(Systemic I+T)  |  HR < 1 favours the HAIC strategy  "
         "|  estimand: ATO (overlap population)",
         fontsize=8.5, color=COL_TXT2, va="top")

base = os.path.join(OUTPUT_DIR, "HR_forest_ow_vs_systemic_it")
fig.savefig(f"{base}.pdf", bbox_inches="tight", pad_inches=0.08)
fig.savefig(f"{base}.png", dpi=300, bbox_inches="tight", pad_inches=0.08)
plt.close(fig)
print(f"已保存: {base}.pdf / .png")
