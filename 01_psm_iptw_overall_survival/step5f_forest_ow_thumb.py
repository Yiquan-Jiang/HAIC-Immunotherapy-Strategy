#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Thumbnail forest plot for the graphical abstract — 8-group, vs HAIC alone (OW).

Companion to step5f_forest_ow.py. Reads the existing pairwise overlap-weighting
result `ow_vs_haic_alone_8group/ow_forest_data.csv` (no model re-fitting) and
renders the single After-OW (ATO) panel as a compact, legible-when-shrunk forest:
  * 8 rows (HAIC alone reference + 7 comparison groups incl. Systemic I+T)
  * rows ordered by ascending OW HR — same order as the full step5f figure
  * larger, bolder fonts and thicker CI bars so it stays legible at small sizes
  * minimal text columns: group label + HR (95% CI) — no n, no P column
  * significance colour (p < 0.05 → red, else gray), matching step5f's raw-p logic
    (the vs-HAIC-alone OW track has no Holm correction)
  * reference row marked with a green diamond at HR = 1

This is the 8-group analogue of step5c_forest_vs_HAIC_alone_thumb.py, pointed at
the overlap-weighting track (the preferred estimate for the Systemic I+T contrast,
whose joint-IPTW result is overlap-deficient).

Outputs into the same ow_vs_haic_alone_8group/ folder:
  HR_forest_ow_vs_haic_alone_thumb.{png,pdf}
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

BASE_DIR = (
    "/Users/yqj/Nutstore Files/我的坚果云/Liver_tumor_big_data/"
    "FIRST_LINE_HAIC_2025-12-30/HAIC_NO_TACE_4_TIDY/update_group_7"
)
RES_DIR = os.path.join(BASE_DIR, "results", "ow_vs_haic_alone_8group")
FIG_DIR = os.path.join(BASE_DIR, "figures", "ow_vs_haic_alone_8group")
CSV_IN  = os.path.join(RES_DIR, "ow_forest_data.csv")

COL_NS  = "#3C5488"
COL_SIG = "#E64B35"
COL_REF = "#00A087"
COL_TXT = "#1F1F1F"

REF_GROUP = "HAIC_alone"

SHORT_LABELS = {
    "HAIC_alone":          "HAIC alone",
    "HAIC+I_concurrent":   "HAIC + I (conc.)",
    "HAIC_then_I":         "HAIC → I",
    "HAIC+T_concurrent":   "HAIC + T (conc.)",
    "HAIC_then_T":         "HAIC → T",
    "HAIC+I+T_concurrent": "HAIC + I + T (conc.)",
    "HAIC_then_I+T":       "HAIC → I + T",
    "Systemic_I+T":        "Systemic I + T",
}

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


def main():
    fd = pd.read_csv(CSV_IN)

    # Comparison rows ordered by ascending OW HR (same as the full step5f figure).
    comp = fd[fd["group"] != REF_GROUP].sort_values("HR")

    rows = [dict(
        label=SHORT_LABELS[REF_GROUP], HR=1.0, lo=np.nan, hi=np.nan,
        sig=False, is_ref=True,
    )]
    for _, r in comp.iterrows():
        rows.append(dict(
            label=SHORT_LABELS[r["group"]],
            HR=float(r["HR"]), lo=float(r["CI_lower"]), hi=float(r["CI_upper"]),
            sig=bool(r["p"] < 0.05),
            is_ref=False,
        ))

    n = len(rows)

    # Figure sized so 8 rows stay tall enough to keep fonts/lines big when shrunk.
    fig_w, fig_h = 4.6, 4.0
    fig = plt.figure(figsize=(fig_w, fig_h))

    # Layout: [label column | forest | HR text]
    label_frac = 0.30
    hr_frac    = 0.27
    forest_frac = 1.0 - label_frac - hr_frac

    bot = 0.10
    top = 0.92
    height = top - bot

    ax_l = fig.add_axes([0.00, bot, label_frac, height])
    ax_f = fig.add_axes([label_frac, bot, forest_frac, height])
    ax_r = fig.add_axes([label_frac + forest_frac, bot, hr_frac, height])

    for a in (ax_l, ax_r):
        a.set_xlim(0, 1)
        a.set_ylim(-0.5, n - 0.5)
        a.axis("off")

    # Range widened vs the 7-group thumb so Systemic I+T (HR>1) sits inside the axis.
    X_LO, X_HI = 0.35, 1.7
    ax_f.set_xlim(X_LO, X_HI)
    ax_f.set_ylim(-0.5, n - 0.5)
    ax_f.set_xscale("log")
    for s in ("top", "right", "left"):
        ax_f.spines[s].set_visible(False)
    ax_f.spines["bottom"].set_linewidth(1.0)
    ax_f.spines["bottom"].set_color(COL_TXT)
    ax_f.set_yticks([])

    # Reference line at HR = 1
    ax_f.axvline(x=1.0, color="#777777", linestyle=(0, (4, 2)),
                 linewidth=1.0, alpha=0.85, zorder=1)

    for i, rp in enumerate(rows):
        y = n - 1 - i
        is_ref = rp["is_ref"]
        col = COL_REF if is_ref else (COL_SIG if rp["sig"] else COL_NS)

        # Label
        ax_l.text(0.97, y, rp["label"], ha="right", va="center",
                  fontsize=11,
                  fontweight="bold" if is_ref else "semibold",
                  color=col if is_ref else COL_TXT)

        if is_ref:
            ax_f.scatter([1.0], [y], marker="D", s=130,
                         color=col, zorder=5, linewidths=0.6,
                         edgecolors="white")
            ax_r.text(0.03, y, "1.00 (ref)", ha="left", va="center",
                      fontsize=10.5, fontweight="bold", color=col)
            continue

        hr, lo, hi = rp["HR"], rp["lo"], rp["hi"]
        # clip CI to axis range for drawing
        lo_c = max(lo, X_LO * 1.02)
        hi_c = min(hi, X_HI * 0.98)
        ax_f.plot([lo_c, hi_c], [y, y],
                  color=col, linewidth=3.2, solid_capstyle="round", zorder=3)
        # CI end caps
        for xc in (lo_c, hi_c):
            ax_f.plot([xc, xc], [y - 0.18, y + 0.18],
                      color=col, linewidth=2.0, zorder=3)
        # diamond
        ax_f.scatter([hr], [y], marker="D", s=150,
                     color=col, zorder=5, linewidths=0.7,
                     edgecolors="white")

        hr_str = f"{hr:.2f} ({lo:.2f}–{hi:.2f})"
        ax_r.text(0.03, y, hr_str, ha="left", va="center",
                  fontsize=10.5,
                  fontweight="bold" if rp["sig"] else "normal",
                  color=col if rp["sig"] else COL_TXT)

    # X axis ticks — sparse and big; pin both major + minor on the log axis
    # otherwise matplotlib silently overlays scientific-notation minor labels.
    xticks = [0.5, 1.0, 1.5]
    ax_f.xaxis.set_major_locator(mticker.FixedLocator(xticks))
    ax_f.xaxis.set_major_formatter(mticker.FixedFormatter(["0.5", "1", "1.5"]))
    ax_f.xaxis.set_minor_locator(mticker.NullLocator())
    ax_f.xaxis.set_minor_formatter(mticker.NullFormatter())
    ax_f.tick_params(axis="x", which="major", length=4, width=1.0,
                     color=COL_TXT, labelsize=10.5, labelcolor=COL_TXT)
    ax_f.set_xlabel("Hazard ratio (vs HAIC alone)",
                    fontsize=11, color=COL_TXT, labelpad=40)

    # Favours arrows sit just below the tick labels, above the x-label.
    arrow_y = -0.13
    ax_f.annotate("", xy=(0.30, arrow_y), xytext=(0.50, arrow_y),
                  xycoords="axes fraction",
                  arrowprops=dict(arrowstyle="-|>", lw=1.2,
                                  color="#444444"))
    ax_f.annotate("", xy=(0.70, arrow_y), xytext=(0.50, arrow_y),
                  xycoords="axes fraction",
                  arrowprops=dict(arrowstyle="-|>", lw=1.2,
                                  color="#444444"))
    ax_f.text(0.28, arrow_y, "row better",
              transform=ax_f.transAxes,
              ha="right", va="center", fontsize=9, color="#444444",
              fontstyle="italic")
    ax_f.text(0.72, arrow_y, "HAIC alone better",
              transform=ax_f.transAxes,
              ha="left", va="center", fontsize=9, color="#444444",
              fontstyle="italic")

    out_base = os.path.join(FIG_DIR, "HR_forest_ow_vs_haic_alone_thumb")
    for ext in (".png", ".pdf"):
        fig.savefig(out_base + ext,
                    dpi=600 if ext == ".png" else 300,
                    bbox_inches="tight", pad_inches=0.08)
        print(f"Saved: {out_base + ext}")
    plt.close(fig)


if __name__ == "__main__":
    main()
