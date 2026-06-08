#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Publication-quality forest plot: each HAIC-containing group vs Systemic I+T (reference).
Two panels side by side: Left = Before IPTW, Right = After CBPS-IPTW.
HR = h(row group) / h(Systemic_I+T).  HR < 1 → row group better than Systemic I+T.
IPTW via WeightIt CBPS + weight truncation (P1/P99) + full-model Cox with robust SE.
Style: Nature / Cell — table-style layout with aligned columns.

This script is ALWAYS 8-group (Systemic I+T is the reference arm):
reads analysis_ready_8group.csv, reads weighted-Cox input from
results/psm_vs_template_8group/survival_gps_final.csv, writes figures to
figures/psm_pub_quality_8group/.
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
SFX = "_8group"
DATA_CSV = "analysis_ready_8group.csv"
DATA_DIR   = os.path.join(BASE_DIR, "data")
OUTPUT_DIR = os.path.join(BASE_DIR, "figures", "psm_pub_quality" + SFX)
os.makedirs(OUTPUT_DIR, exist_ok=True)

REF_GROUP = "Systemic_I+T"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 9,
    "axes.labelsize": 10,
    "axes.titlesize": 11,
    "xtick.labelsize": 8.5,
    "ytick.labelsize": 8.5,
    "legend.fontsize": 8,
    "axes.linewidth": 0.7,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "xtick.major.width": 0.7,
    "ytick.major.width": 0.7,
    "xtick.major.size": 3.5,
    "ytick.major.size": 3.5,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.05,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "axes.grid": False,
})

COL_NS   = "#3C5488"
COL_SIG  = "#E64B35"
COL_REF  = "#00A087"
COL_GRID = "#EDEDED"
COL_RULE = "#333333"
COL_TXT  = "#222222"
COL_TXT2 = "#555555"

GROUP_ORDER = [
    "HAIC_alone",
    "HAIC+I_concurrent",
    "HAIC_then_I",
    "HAIC+T_concurrent",
    "HAIC_then_T",
    "HAIC+I+T_concurrent",
    "HAIC_then_I+T",
    "Systemic_I+T",
]

GROUP_COLORS = {
    "HAIC_alone":            "#00A087",
    "HAIC+I_concurrent":     "#3C5488",
    "HAIC_then_I":           "#4DBBD5",
    "HAIC+T_concurrent":     "#E64B35",
    "HAIC_then_T":           "#F39B7F",
    "HAIC+I+T_concurrent":   "#8491B4",
    "HAIC_then_I+T":         "#91D1C2",
    "Systemic_I+T":          "#009E73",
}

GROUP_LABELS = {
    "HAIC_alone":            "HAIC alone",
    "HAIC+I_concurrent":     "HAIC + Immunotherapy (conc.)",
    "HAIC_then_I":           "HAIC → Immunotherapy (seq.)",
    "HAIC+T_concurrent":     "HAIC + Targeted therapy (conc.)",
    "HAIC_then_T":           "HAIC → Targeted therapy (seq.)",
    "HAIC+I+T_concurrent":   "HAIC + I + T (concurrent)",
    "HAIC_then_I+T":         "HAIC → I + T (sequential)",
    "Systemic_I+T":          "Systemic I+T",
}


# ════════════════════════════════════════════════════════════════════
# 1. Compute Before-IPTW HR (unadjusted Cox, P from Cox Wald test)
# ════════════════════════════════════════════════════════════════════
print("1. 计算 IPTW 前 HR（各组 vs Systemic_I+T, 未调整 Cox）...")

df = pd.read_csv(os.path.join(DATA_DIR, DATA_CSV))
df = df[df["os_months"] >= 0].copy()
df["group"] = df["main_group"]
df["death_status"] = df["death_status"].map(
    {"Yes": 1, "No": 0, "1": 1, "0": 0, 1: 1, 0: 0}
).fillna(0).astype(int)

ref_data = df[df["group"] == REF_GROUP]
others = [g for g in GROUP_ORDER if g != REF_GROUP]

before_rows = []
for g in others:
    sub = df[df["group"] == g]
    tmp = pd.concat([ref_data, sub])[["os_months", "death_status"]].copy()
    tmp["treat"] = [0] * len(ref_data) + [1] * len(sub)
    cph = CoxPHFitter()
    cph.fit(tmp, duration_col="os_months", event_col="death_status")
    hr    = float(np.exp(cph.params_["treat"]))
    ci_lo = float(np.exp(cph.confidence_intervals_["95% lower-bound"]["treat"]))
    ci_hi = float(np.exp(cph.confidence_intervals_["95% upper-bound"]["treat"]))
    p = float(cph.summary["p"]["treat"])
    before_rows.append(dict(
        group=g, label=GROUP_LABELS[g],
        HR=hr, CI_lower=ci_lo, CI_upper=ci_hi, P_value=p,
        N=len(sub), N_ref=len(ref_data), is_ref=False, phase="before",
    ))
    print(f"  Before: {g} | n={len(sub)} | HR={hr:.2f} ({ci_lo:.2f}-{ci_hi:.2f}) P={p:.4f}")

# ════════════════════════════════════════════════════════════════════
# 2. Read After-IPTW HR from survival_gps_final.csv and re-level to Systemic_I+T
# ════════════════════════════════════════════════════════════════════
# survival_gps_final.csv stores HR(g vs HAIC_alone) for the weighted full-model
# Cox.  To express each HAIC group vs Systemic_I+T as reference we re-level the
# SAME weighted Cox via the ratio identity
#   HR(g vs Systemic_I+T) = HR(g vs HAIC_alone) / HR(Systemic_I+T vs HAIC_alone)
# applied to the point estimate and both CI bounds (multiplicative).  No new
# statistics are fitted; this is the re-leveling of the same weighted Cox.
print("\n2. 读取 CBPS-IPTW 加权后 HR (全模型 Cox, robust SE), 重新设定 Systemic_I+T 为参照...")

GPS_RES_DIR = os.path.join(BASE_DIR, "results", "psm_vs_template" + SFX)
surv = pd.read_csv(os.path.join(GPS_RES_DIR, "survival_gps_final.csv"))

ref_mask = surv["Group"] == REF_GROUP
if not ref_mask.any():
    raise ValueError(f"survival_gps_final.csv 中未找到参照组 {REF_GROUP} (vs HAIC_alone) 的加权 Cox 结果")
ref_r = surv.loc[ref_mask].iloc[0]
ref_hr = float(ref_r["HR"])              # HR(Systemic_I+T vs HAIC_alone)
n_ref_after = int(ref_r["N"])            # weighted N of Systemic_I+T arm

after_rows = []
for g in others:
    if g == "HAIC_alone":
        # HAIC_alone is the original CSV reference (HR == 1 vs itself); re-level
        # to Systemic_I+T: HR(HAIC_alone vs Systemic_I+T) = 1 / HR(Sys vs HAIC_alone)
        hr_src, lo_src, hi_src = 1.0, 1.0, 1.0
        p = float(ref_r["P_value"])
        p_holm = float(ref_r["P_holm"]) if "P_holm" in ref_r.index else p
        n_other = int(ref_r["N_ref"])    # N of HAIC_alone arm
    else:
        mask = surv["Group"] == g
        if not mask.any():
            print(f"  ⚠ 未找到 {g} 的 IPTW 结果")
            continue
        r = surv.loc[mask].iloc[0]
        hr_src = float(r["HR"])
        lo_src = float(r["CI_lower"])
        hi_src = float(r["CI_upper"])
        p  = float(r["P_value"])
        p_holm = float(r["P_holm"]) if "P_holm" in r.index else p
        n_other = int(r["N"])

    # Re-level vs Systemic_I+T (divide by HR(Systemic_I+T vs HAIC_alone))
    hr = hr_src / ref_hr
    lo = lo_src / ref_hr
    hi = hi_src / ref_hr

    after_rows.append(dict(
        group=g, label=GROUP_LABELS[g],
        HR=hr, CI_lower=lo, CI_upper=hi, P_value=p, P_holm=p_holm,
        N=n_other, N_ref=n_ref_after, is_ref=False, phase="after",
    ))
    p_fmt = "< 0.0001" if p < 0.0001 else f"{p:.4f}"
    p_holm_fmt = "< 0.0001" if p_holm < 0.0001 else f"{p_holm:.4f}"
    print(f"  After:  {g} | n={n_other} | HR={hr:.2f} ({lo:.2f}-{hi:.2f}) P={p_fmt} P_holm={p_holm_fmt}")

# Save CSV
all_rows = before_rows + after_rows
pd.DataFrame(all_rows).to_csv(
    os.path.join(OUTPUT_DIR, "forest_vs_systemic_it.csv"), index=False)

# ════════════════════════════════════════════════════════════════════
# 3. Draw side-by-side forest plots
# ════════════════════════════════════════════════════════════════════
print("\n3. 绘制 Forest Plot...")


def build_plot_rows(data_rows):
    rows = [dict(
        group=REF_GROUP, label=GROUP_LABELS[REF_GROUP],
        HR=1.0, CI_lower=np.nan, CI_upper=np.nan, P_value=np.nan,
        N=np.nan, N_ref=np.nan, is_ref=True,
    )]
    rows.extend(data_rows)
    return rows


def draw_forest_panel(fig, rect, rows, panel_title, subtitle):
    n = len(rows)

    left_frac  = 0.40
    mid_frac   = 0.28
    right_frac = 0.32

    ax_left   = fig.add_axes([rect[0],
                               rect[1],
                               rect[2] * left_frac,
                               rect[3]])
    ax_forest = fig.add_axes([rect[0] + rect[2] * left_frac,
                               rect[1],
                               rect[2] * mid_frac,
                               rect[3]])
    ax_right  = fig.add_axes([rect[0] + rect[2] * (left_frac + mid_frac),
                               rect[1],
                               rect[2] * right_frac,
                               rect[3]])

    for a in (ax_left, ax_right):
        a.set_xlim(0, 1)
        a.set_ylim(-0.5, n - 0.5)
        a.axis("off")

    X_LO, X_HI = 0.25, 4.0
    ax_forest.set_ylim(-0.5, n - 0.5)
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

    hdr_y = n - 0.5 + 0.35
    for a in (ax_left, ax_forest, ax_right):
        a.set_ylim(-0.5, n - 0.5 + 0.7)

    ax_left.text(0.04, hdr_y, "Treatment group", ha="left", va="center",
                 fontsize=9, fontweight="bold", color=COL_TXT)
    ax_left.text(0.92, hdr_y, "n", ha="center", va="center",
                 fontsize=9, fontweight="bold", color=COL_TXT)
    ax_right.text(0.04, hdr_y, "HR (95% CI)", ha="left", va="center",
                  fontsize=9, fontweight="bold", color=COL_TXT)
    ax_right.text(0.85, hdr_y, "P value", ha="center", va="center",
                  fontsize=9, fontweight="bold", color=COL_TXT)

    rule_kw = dict(color=COL_RULE, linewidth=0.8, clip_on=False)
    for a in (ax_left, ax_forest, ax_right):
        a.axhline(y=n - 0.5 + 0.08, **rule_kw)
        a.axhline(y=-0.5, **rule_kw)

    for i, rp in enumerate(rows):
        y = n - 1 - i
        is_ref = rp["is_ref"]
        hr = rp["HR"]
        lo, hi = rp["CI_lower"], rp["CI_upper"]
        p = rp["P_value"]

        label_style = dict(fontsize=8.5, va="center", color=COL_TXT)
        if is_ref:
            label_style["fontweight"] = "bold"
            label_style["color"] = COL_REF

        ax_left.text(0.04, y, rp["label"], ha="left", **label_style)

        if is_ref:
            ax_left.text(0.92, y, "—", ha="center", va="center",
                         fontsize=8.5, color=COL_REF, fontweight="bold")
            ax_right.text(0.04, y, "1.00 (reference)", ha="left", va="center",
                          fontsize=8.5, color=COL_REF, fontweight="bold")
            ax_right.text(0.85, y, "—", ha="center", va="center",
                          fontsize=8.5, color=COL_REF, fontweight="bold")
            ax_forest.scatter([1.0], [y], marker="D", s=65,
                              color=COL_REF, zorder=5, linewidths=0.4,
                              edgecolors="white")
            continue

        sig = p < 0.05 if not np.isnan(p) else False
        col = COL_SIG if sig else COL_NS

        n_str = f"{int(rp['N'])}"
        ax_left.text(0.92, y, n_str, ha="center", va="center",
                     fontsize=8.5, color=COL_TXT)

        ci_lo_c = max(lo, X_LO * 1.02)
        ci_hi_c = min(hi, X_HI * 0.98)
        ax_forest.plot([ci_lo_c, ci_hi_c], [y, y],
                       color=col, linewidth=1.8, solid_capstyle="round", zorder=3)
        for xc in (ci_lo_c, ci_hi_c):
            ax_forest.plot([xc, xc], [y - 0.14, y + 0.14],
                           color=col, linewidth=1.0, zorder=3)

        weight = min(max(rp["N"] / 500, 0.5), 1.8)
        diamond_s = 60 * weight
        ax_forest.scatter([hr], [y], marker="D", s=diamond_s,
                          color=col, zorder=5, linewidths=0.4,
                          edgecolors="white")

        hr_str = f"{hr:.2f} ({lo:.2f}–{hi:.2f})"
        p_str  = "< 0.001" if p < 0.001 else f"{p:.3f}"
        txt_col = "#B03000" if sig else COL_TXT
        ax_right.text(0.04, y, hr_str, ha="left", va="center",
                      fontsize=8.5, color=txt_col,
                      fontweight="bold" if sig else "normal")
        ax_right.text(0.85, y, p_str, ha="center", va="center",
                      fontsize=8.5, color=txt_col,
                      fontweight="bold" if sig else "normal")

    xticks = [0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0]
    ax_forest.set_xticks(xticks)
    ax_forest.set_xticklabels([str(x) for x in xticks], fontsize=8)
    ax_forest.set_xlabel("Hazard ratio", fontsize=9.5, labelpad=10)

    ax_forest.text(0.18, -0.10, "Favours\nrow group",
                   ha="center", va="top", fontsize=7, color=COL_TXT2,
                   transform=ax_forest.transAxes, linespacing=0.9)
    ax_forest.text(0.82, -0.10, "Favours\nSystemic I+T",
                   ha="center", va="top", fontsize=7, color=COL_TXT2,
                   transform=ax_forest.transAxes, linespacing=0.9)

    ax_forest.set_title(panel_title, fontsize=11, fontweight="bold",
                        loc="center", pad=18)
    ax_forest.text(0.5, 1.06, subtitle, fontsize=7.5, color=COL_TXT2,
                   ha="center", va="bottom", transform=ax_forest.transAxes)

    return ax_forest


rows_before = build_plot_rows(before_rows)
rows_after  = build_plot_rows(after_rows)

fig_w = 10.0
fig_h = 0.48 * len(rows_before) + 1.8
fig = plt.figure(figsize=(fig_w, fig_h))

panel_gap = 0.03
panel_w = (1.0 - panel_gap) / 2.0
bot = 0.15
top_h = 0.73

ax_l = draw_forest_panel(
    fig, [0.0, bot, panel_w, top_h],
    rows_before,
    "A  Before IPTW",
    "Unadjusted Cox regression"
)
ax_r = draw_forest_panel(
    fig, [panel_w + panel_gap, bot, panel_w, top_h],
    rows_after,
    "B  After CBPS-IPTW",
    "CBPS + weight truncation (P1/P99) + full-model Cox (robust SE)"
)

leg_handles = [
    mlines.Line2D([], [], color=COL_NS, marker="D", markersize=6,
                  markerfacecolor=COL_NS, markeredgecolor="white",
                  markeredgewidth=0.4, linewidth=1.6, label="P ≥ 0.05"),
    mlines.Line2D([], [], color=COL_SIG, marker="D", markersize=6,
                  markerfacecolor=COL_SIG, markeredgecolor="white",
                  markeredgewidth=0.4, linewidth=1.6, label="P < 0.05"),
    mlines.Line2D([], [], color=COL_REF, marker="D", markersize=6,
                  markerfacecolor=COL_REF, markeredgecolor="white",
                  markeredgewidth=0.4, linewidth=0, label="Reference (Systemic I+T)"),
]
fig.legend(handles=leg_handles, loc="lower center",
           bbox_to_anchor=(0.5, 0.01), ncol=3,
           fontsize=8.5, frameon=False, handlelength=1.8,
           columnspacing=1.5)

fig.text(0.02, 0.96,
         "HR vs Systemic I+T (does adding HAIC to systemic I+T improve OS?)",
         fontsize=12, fontweight="bold", color=COL_TXT, va="top")
fig.text(0.02, 0.93,
         "HR = h(row group) / h(Systemic I+T)  |  HR < 1 favours row group",
         fontsize=8.5, color=COL_TXT2, va="top")

base = os.path.join(OUTPUT_DIR, "HR_forest_vs_systemic_it")
fig.savefig(f"{base}.pdf", bbox_inches="tight", pad_inches=0.08)
fig.savefig(f"{base}.png", dpi=300, bbox_inches="tight", pad_inches=0.08)
plt.close(fig)
print(f"\n已保存: {base}.pdf / .png")
