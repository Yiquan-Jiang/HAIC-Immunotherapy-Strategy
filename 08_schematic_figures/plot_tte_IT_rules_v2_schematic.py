"""
Schematic of the IT-Rules v2 trigger logic (two cohorts).

Layer 1 (pre-HAIC, evaluated per cycle from cycle 3):
  - AFP / PLR / SII / NLR / PIVKA / Distant-Meta / Lymph-Node
Layer 2 (post-HAIC, only if untriggered in Layer 1):
  - AFP > 20 ng/mL

Produces a single figure that:
  (A) Lays out the per-cohort rule panels (A: I+T cohort_3matched, B: I cohort_7group_psm02)
  (B) Shows the timeline / cycle structure
  (C) Shows the Clone-Censor-Weight skeleton with both arms

Output:
  output/schematic_figures/IT_rules_v2_schematic.{pdf,png}
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle
from matplotlib.lines import Line2D

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.normpath(os.path.join(SCRIPT_DIR, "..", ".."))
OUT_DIR = os.path.join(PROJECT_ROOT, "output", "schematic_figures")
os.makedirs(OUT_DIR, exist_ok=True)

# ---- color palette ----------------------------------------------------------
C_HAIC      = "#2C3E50"
C_TIMELINE  = "#34495E"
C_LAYER1    = "#E64B35"
C_LAYER2    = "#F39C12"
C_NEVER     = "#7F8C8D"
C_DYN       = "#E64B35"
C_EARLY     = "#3C5488"
C_BG_A      = "#FFF6F2"
C_BG_B      = "#F2F6FF"
C_TEXT      = "#1B1B1B"
C_GRAY      = "#888888"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 9,
    "axes.linewidth": 0.6,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


def fancy_box(ax, x, y, w, h, text, fc, ec=None, fontsize=8.5, fontweight="normal",
              text_color=C_TEXT, pad=0.35, align="center"):
    if ec is None:
        ec = fc
    box = FancyBboxPatch((x, y), w, h,
                         boxstyle=f"round,pad={pad},rounding_size=0.18",
                         linewidth=0.8, facecolor=fc, edgecolor=ec)
    ax.add_patch(box)
    ha = "center" if align == "center" else "left"
    tx = x + w / 2 if align == "center" else x + 0.15
    ax.text(tx, y + h / 2, text,
            ha=ha, va="center", fontsize=fontsize,
            color=text_color, fontweight=fontweight)
    return box


def arrow(ax, x0, y0, x1, y1, color="#444444", lw=1.0, style="-|>",
          mutation_scale=10):
    a = FancyArrowPatch((x0, y0), (x1, y1), arrowstyle=style, color=color,
                        lw=lw, mutation_scale=mutation_scale)
    ax.add_patch(a)


def cohort_panel(ax, cohort_label, sub_title,
                 layer1_lines, layer2_line, never_line,
                 cohort_n, triggered_n, never_n, hr_text,
                 bg_color):
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 12)
    ax.axis("off")
    ax.add_patch(Rectangle((0, 0), 10, 12, facecolor=bg_color,
                           edgecolor="none", zorder=0))

    # Title bar
    ax.text(0.2, 11.55, cohort_label, fontsize=12.5, fontweight="bold",
            color=C_TEXT, va="center")
    ax.text(0.2, 11.00, sub_title, fontsize=9, color=C_GRAY, va="center")
    ax.plot([0.2, 9.8], [10.65, 10.65], color="#999", lw=0.6)

    # Cycle timeline
    cyc_y = 9.7
    cyc_x = [0.6 + i * 1.45 for i in range(6)]
    ax.plot([cyc_x[0] - 0.25, cyc_x[-1] + 0.4], [cyc_y, cyc_y],
            color=C_TIMELINE, lw=1.1)
    for i, x in enumerate(cyc_x, start=1):
        ax.plot(x, cyc_y, marker="o",
                color=C_HAIC, markersize=7, zorder=3)
        ax.text(x, cyc_y - 0.42, f"HAIC\nC{i}", ha="center", va="top",
                fontsize=7.5, color=C_TEXT)
    ax.text(cyc_x[-1] + 0.95, cyc_y - 0.05, "post-HAIC",
            ha="left", va="center", fontsize=8, color=C_LAYER2,
            fontweight="bold")
    arrow(ax, cyc_x[-1] + 0.45, cyc_y, cyc_x[-1] + 0.95, cyc_y,
          color=C_LAYER2, lw=1.0)

    # Layer 1 box
    L1_top = 8.7
    fancy_box(ax, 0.3, L1_top - 0.6, 9.4, 0.6,
              "Layer 1  (per-cycle, evaluated from C3)",
              fc=C_LAYER1, text_color="white", fontweight="bold", fontsize=9)
    # Each rule
    rule_y = L1_top - 1.2
    for line in layer1_lines:
        fancy_box(ax, 0.55, rule_y - 0.55, 8.9, 0.55, line,
                  fc="white", ec=C_LAYER1, fontsize=8.3, align="left")
        rule_y -= 0.65

    # Layer 2 box
    L2_top = rule_y - 0.15
    fancy_box(ax, 0.3, L2_top - 0.55, 9.4, 0.55,
              "Layer 2  (post-HAIC, only if Layer 1 didn't fire)",
              fc=C_LAYER2, text_color="white", fontweight="bold", fontsize=9)
    fancy_box(ax, 0.55, L2_top - 1.15, 8.9, 0.55, layer2_line,
              fc="white", ec=C_LAYER2, fontsize=8.3, align="left")

    # Never triggered
    N_top = L2_top - 1.45
    fancy_box(ax, 0.3, N_top - 0.55, 9.4, 0.55,
              "Never triggered (continue HAIC alone)",
              fc=C_NEVER, text_color="white", fontweight="bold", fontsize=9)
    fancy_box(ax, 0.55, N_top - 1.15, 8.9, 0.55, never_line,
              fc="white", ec=C_NEVER, fontsize=8.3, align="left")

    # bottom stat strip
    stat_y = 0.55
    fancy_box(ax, 0.3, stat_y - 0.45, 9.4, 0.85,
              f"N eligible = {cohort_n}    |    triggered = {triggered_n}    |    "
              f"never = {never_n}    |    {hr_text}",
              fc="white", ec="#999", fontsize=8.5, fontweight="bold")


def ccw_panel(ax):
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6)
    ax.axis("off")

    # Eligible -> Clone
    fancy_box(ax, 4.0, 5.05, 2.0, 0.7, "Eligible patient",
              fc="white", ec="#555", fontweight="bold", fontsize=9)
    arrow(ax, 5.0, 5.05, 5.0, 4.55, lw=1.2)
    fancy_box(ax, 3.6, 3.85, 2.8, 0.7, "Clone (×2)",
              fc="#222", text_color="white", fontweight="bold", fontsize=9)

    # Two arms
    arrow(ax, 4.5, 3.85, 2.4, 3.05, lw=1.2)
    arrow(ax, 5.5, 3.85, 7.6, 3.05, lw=1.2)

    fancy_box(ax, 0.4, 2.35, 4.0, 0.7, "Arm A: Adaptive On Demand",
              fc=C_DYN, text_color="white", fontweight="bold", fontsize=9.5)
    fancy_box(ax, 5.6, 2.35, 4.0, 0.7, "Arm B: Early Combination (≤14 d)",
              fc=C_EARLY, text_color="white", fontweight="bold", fontsize=9.5)

    arm_a_lines = [
        "Censor if add-on > trigger_day + 90 d",
        "Censor if add-on given before trigger",
        "Censor untriggered patients who got add-on",
    ]
    y = 1.8
    for ln in arm_a_lines:
        fancy_box(ax, 0.4, y - 0.4, 4.0, 0.4, "• " + ln,
                  fc="white", ec=C_DYN, fontsize=7.8, align="left")
        y -= 0.45

    arm_b_lines = [
        "Censor at day 14 if add-on > 14 d",
        "Otherwise follow OS to event/censoring",
    ]
    y = 1.8
    for ln in arm_b_lines:
        fancy_box(ax, 5.6, y - 0.4, 4.0, 0.4, "• " + ln,
                  fc="white", ec=C_EARLY, fontsize=7.8, align="left")
        y -= 0.45

    # IPCW + Cox
    arrow(ax, 2.4, 0.45, 4.7, 0.05, lw=1.2)
    arrow(ax, 7.6, 0.45, 5.3, 0.05, lw=1.2)
    fancy_box(ax, 3.0, -0.55, 4.0, 0.55,
              "Stabilized IPCW  →  weighted Cox / RMST",
              fc="#222", text_color="white", fontweight="bold", fontsize=9)


def main():
    fig = plt.figure(figsize=(15, 12))

    # 2x2 grid: top row = two cohort rule panels; bottom = CCW panel spanning both
    gs = fig.add_gridspec(
        nrows=2, ncols=2,
        height_ratios=[1.65, 1.0], width_ratios=[1, 1],
        hspace=0.18, wspace=0.10,
        left=0.03, right=0.985, top=0.965, bottom=0.05,
    )

    axA = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[0, 1])
    axC = fig.add_subplot(gs[1, :])

    layer1_A = [
        "1. AFP drop < 32.5%   (AFP %Δ > −32.5)",
        "2. PLR > 102.4",
        "3. SII (PLT × NLR) > 390.9",
        "4. PIVKA drop < 45.6%   (PIVKA %Δ > −45.6)",
        "5. Distant metastasis = 1   (fires at C3)",
        "6. Lymph-node metastasis = 1   (fires at C3)",
    ]
    cohort_panel(
        axA,
        cohort_label="Cohort A — HAIC then I+T on demand",
        sub_title="cohort_3matched   (matched IDs 06  +HAIC+I+T concurrent)",
        layer1_lines=layer1_A,
        layer2_line="AFP > 20 ng/mL  OR  PIVKA > 40 mAU/mL after HAIC  →  trigger I+T add-on",
        never_line="No add-on; continue HAIC monotherapy and follow OS",
        cohort_n=1938,
        triggered_n=1328,
        never_n=610,
        hr_text="HR (AoD vs Early) = 0.68 (0.55–0.83), P < 0.001",
        bg_color=C_BG_A,
    )

    layer1_B = [
        "1. AFP drop < 40%   (AFP %Δ > −40)",
        "2. PLR > 98.7",
        "3. NLR > 2.68",
        "4. PIVKA drop < 51.2%   (PIVKA %Δ > −51.2)",
        "5. Distant metastasis = 1   (fires at C3)",
        "6. Lymph-node metastasis = 1   (fires at C3)",
    ]
    cohort_panel(
        axB,
        cohort_label="Cohort B — HAIC then I on demand",
        sub_title="cohort_7group_psm02   (matched IDs 02  +HAIC+I concurrent)",
        layer1_lines=layer1_B,
        layer2_line="AFP > 20 ng/mL  OR  PIVKA > 40 mAU/mL after HAIC  →  trigger I add-on",
        never_line="No add-on; continue HAIC monotherapy and follow OS",
        cohort_n=574,
        triggered_n=377,
        never_n=197,
        hr_text="HR (AoD vs Early) = 0.68 (0.49–0.94), P = 0.020",
        bg_color=C_BG_B,
    )

    ccw_panel(axC)

    fig.suptitle(
        "Target Trial Emulation — IT-Rules v2 trigger logic and "
        "Clone-Censor-Weight design",
        fontsize=13.5, fontweight="bold", color=C_TEXT, y=0.995,
    )

    out_pdf = os.path.join(OUT_DIR, "IT_rules_v2_schematic.pdf")
    out_png = os.path.join(OUT_DIR, "IT_rules_v2_schematic.png")
    fig.savefig(out_pdf, dpi=300, bbox_inches="tight", pad_inches=0.08)
    fig.savefig(out_png, dpi=300, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    print(f"Saved: {out_pdf}")
    print(f"Saved: {out_png}")


if __name__ == "__main__":
    main()
