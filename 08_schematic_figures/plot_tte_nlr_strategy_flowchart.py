#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TTE strategy rules — clean two-arm specification only.
Matches tte_nlr_R_core.R (NLR_BASED_RULES_v2).

  python3 script/plot_tte_nlr_strategy_flowchart.py
"""
from __future__ import annotations
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

SCRIPT_DIR  = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent
DEFAULT_OUT  = PROJECT_ROOT / "output" / "schematic_figures"

CA  = "#3C5488"
CB  = "#E64B35"
CT  = "#00A087"
CCY = "#4DBBD5"
CK  = "#222222"
CG  = "#888888"
CL  = "#C0C0C0"

FW, FH = 7.0, 4.2


def _rc():
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7, "axes.linewidth": 0,
        "figure.dpi": 150, "savefig.dpi": 300,
        "savefig.bbox": "tight", "savefig.pad_inches": 0.03,
        "pdf.fonttype": 42, "ps.fonttype": 42,
    })


def box(ax, x, y, w, h, txt, *,
        fc="white", ec=CK, lw=0.5, fs=6, fw="normal", tc=CK, ls=1.1):
    ax.add_patch(FancyBboxPatch(
        (x - w/2, y - h/2), w, h, boxstyle="round,pad=0.04",
        lw=lw, ec=ec, fc=fc, zorder=2, clip_on=False))
    ax.text(x, y, txt, ha="center", va="center", fontsize=fs,
            fontweight=fw, color=tc, linespacing=ls, zorder=5, clip_on=False)


def dia(ax, cx, cy, w, h, txt, *, fc="#FFF9E6", ec=CK, fs=5.5):
    hw, hh = w/2, h/2
    ax.fill([cx, cx+hw, cx, cx-hw, cx],
            [cy+hh, cy, cy-hh, cy, cy+hh],
            fc=fc, ec=ec, lw=0.5, zorder=2, clip_on=False)
    ax.text(cx, cy, txt, ha="center", va="center", fontsize=fs,
            fontweight="bold", color=CK, zorder=5, clip_on=False)


def arr(ax, x0, y0, x1, y1, *, c=CK, lw=0.5, ms=5.5):
    ax.add_patch(FancyArrowPatch(
        (x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=ms,
        lw=lw, color=c, clip_on=False, zorder=4, shrinkA=0, shrinkB=0))


def ln(ax, x0, y0, x1, y1, *, c=CK, lw=0.4, ls="-"):
    ax.plot([x0, x1], [y0, y1], c=c, lw=lw, ls=ls, clip_on=False, zorder=3)


def lbl(ax, x, y, txt, *, fs=5.5, c=CG, ha="center", fw="normal"):
    ax.text(x, y, txt, ha=ha, va="center", fontsize=fs, fontweight=fw,
            color=c, zorder=5, clip_on=False)


def draw(out_dir: Path) -> None:
    _rc()
    fig, ax = plt.subplots(figsize=(FW, FH))
    ax.set_xlim(0, FW)
    ax.set_ylim(0, FH)
    ax.axis("off")
    fig.patch.set_facecolor("white")

    G = 0.10
    TOP = FH - 0.12

    # ── Title ────────────────────────────────────────────────────────────
    ax.text(FW/2, TOP + 0.05,
            "Target Trial Emulation — Strategy Specification",
            ha="center", va="top", fontsize=9, fontweight="bold", color=CK)

    # ── Divider ──────────────────────────────────────────────────────────
    div_x = 4.55
    ln(ax, div_x, 0.15, div_x, TOP - 0.28, c=CL, lw=0.35, ls=(0, (4, 3)))

    # ══════════════════════════════════════════════════════════════════════
    #  STRATEGY A — Dynamic
    # ══════════════════════════════════════════════════════════════════════
    Ax = 2.0

    # Header
    hdr_y = TOP - 0.40
    hdr_h = 0.26
    ax.add_patch(FancyBboxPatch(
        (0.15, hdr_y - hdr_h/2), div_x - 0.35, hdr_h,
        boxstyle="round,pad=0.04", lw=0, fc=CA, zorder=2, clip_on=False))
    ax.text(div_x/2 - 0.1, hdr_y,
            "Strategy A — Dynamic (NLR-based)",
            ha="center", va="center", fontsize=7, fontweight="bold",
            color="white", zorder=5)

    # Evaluate
    y1 = hdr_y - 0.34
    bh = 0.22
    arr(ax, Ax, hdr_y - hdr_h/2, Ax, y1 + bh/2 + G, c=CA)
    box(ax, Ax, y1, 3.4, bh,
        "Evaluate at pre-HAIC-3 timepoint",
        fs=6, fw="bold", ec=CT, fc="#F0FAF8")

    # Trigger rules
    tr_top = y1 - bh/2 - 0.06
    rule_gap = 0.17
    rules = [
        ("1", "Vascular invasion (PVTT Vp3/4 · HVTT · IVC/RA) or extrahepatic disease"),
        ("2", "NLR ≥ 2.5 (baseline or pre-HAIC-3)"),
        ("3", "Max tumour > 13 cm"),
        ("4", "PIVKA-II > 12,000 mAU/mL"),
        ("5", "AFP < 20 ng/mL"),
    ]
    title_h = 0.16
    block_h = title_h + len(rules) * rule_gap + 0.08
    tr_bot = tr_top - block_h
    tr_w = 4.1

    arr(ax, Ax, y1 - bh/2, Ax, tr_top + G, c=CT)
    ax.add_patch(FancyBboxPatch(
        (Ax - tr_w/2, tr_bot), tr_w, block_h,
        boxstyle="round,pad=0.04", lw=0.5, ec=CT, fc="white",
        zorder=2, clip_on=False))
    ax.text(Ax, tr_top - 0.06,
            "Trigger rules (if any met → add immunotherapy)",
            ha="center", va="top", fontsize=5.5, fontweight="bold",
            color=CT, zorder=5)

    ry = tr_top - title_h
    rx = Ax - tr_w/2 + 0.12
    for num, desc in rules:
        ax.text(rx, ry, num, ha="center", va="center",
                fontsize=5.5, fontweight="bold", color="white", zorder=6,
                bbox=dict(boxstyle="circle,pad=0.06", fc=CT, ec=CT, lw=0))
        ax.text(rx + 0.20, ry, desc, ha="left", va="center",
                fontsize=5.5, color=CK, zorder=5)
        ry -= rule_gap

    # Diamond: any rule met?
    d_w, d_h = 1.1, 0.36
    yd = tr_bot - d_h/2 - 0.12
    arr(ax, Ax, tr_bot, Ax, yd + d_h/2 + G, c=CK)
    dia(ax, Ax, yd, d_w, d_h, "Any rule met?", fs=5.5)

    # YES → right
    xr = 4.05
    ln(ax, Ax + d_w/2, yd, xr - 0.20, yd)
    arr(ax, xr - 0.20, yd, xr, yd, c=CT)
    lbl(ax, Ax + d_w/2 + 0.12, yd + 0.08, "Yes", fs=5, fw="bold", c=CT)
    box(ax, xr + 0.30, yd, 0.55, 0.26,
        "Add\nimmuno.", fs=5, fw="bold", ec=CT, fc="#F0FAF8")

    # NO → down
    lbl(ax, Ax + 0.12, yd - d_h/2 - 0.06, "No", fs=5, fw="bold", c=CA)

    # ── Exemption diamond (two OR conditions, same level) ────────────────
    lbl(ax, Ax + 0.12, yd - d_h/2 - 0.06, "No", fs=5, fw="bold", c=CA)

    ye = yd - 0.60
    arr(ax, Ax, yd - d_h/2, Ax, ye + d_h*0.7 + G, c=CA)

    # Taller diamond to fit two-line text
    de_w, de_h = 2.4, 0.52
    ax.fill([Ax, Ax+de_w/2, Ax, Ax-de_w/2, Ax],
            [ye+de_h/2, ye, ye-de_h/2, ye, ye+de_h/2],
            fc="#EBF7FB", ec=CCY, lw=0.5, zorder=2, clip_on=False)
    ax.text(Ax, ye + 0.07,
            "Cycle ≥ 4: NLR < 1.8?",
            ha="center", va="center", fontsize=5.5, fontweight="bold",
            color=CK, zorder=5, clip_on=False)
    ax.text(Ax, ye - 0.10,
            "OR  AFP drop > 90% AND AFP < 20 ng/mL?",
            ha="center", va="center", fontsize=5.0, fontweight="bold",
            color=CK, zorder=5, clip_on=False)

    # YES → exempt
    ln(ax, Ax + de_w/2, ye, xr - 0.20, ye)
    arr(ax, xr - 0.20, ye, xr, ye, c=CCY)
    lbl(ax, Ax + de_w/2 + 0.13, ye + 0.08, "Yes", fs=5, fw="bold", c=CCY)
    box(ax, xr + 0.30, ye, 0.55, 0.22,
        "Exempt", fs=5.5, ec=CCY, fc="#EBF7FB")

    # NO → fallback
    lbl(ax, Ax + 0.12, ye - de_h/2 - 0.06, "No", fs=5, fw="bold", c=CA)
    y_fb = ye - 0.52
    arr(ax, Ax, ye - de_h/2, Ax, y_fb + 0.10 + G, c=CA)
    box(ax, Ax, y_fb, 2.8, 0.20,
        "Fallback: add immunotherapy at cycle k",
        fs=5.5, fw="bold", ec=CA, fc="#EEF2F9")

    # ══════════════════════════════════════════════════════════════════════
    #  STRATEGY B — Early Combination
    # ══════════════════════════════════════════════════════════════════════
    Bx = 5.70

    ax.add_patch(FancyBboxPatch(
        (div_x + 0.15, hdr_y - hdr_h/2), FW - div_x - 0.35, hdr_h,
        boxstyle="round,pad=0.04", lw=0, fc=CB, zorder=2, clip_on=False))
    ax.text((div_x + FW) / 2, hdr_y,
            "Strategy B — Early Combination",
            ha="center", va="center", fontsize=7, fontweight="bold",
            color="white", zorder=5)

    yb1 = hdr_y - 0.40
    arr(ax, Bx, hdr_y - hdr_h/2, Bx, yb1 + 0.15 + G, c=CB)
    box(ax, Bx, yb1, 1.9, 0.28,
        "Add immunotherapy\nwithin 14 d of first HAIC",
        fs=6, fw="bold", ec=CB, fc="#FFF5F0")

    ybd = yb1 - 0.42
    arr(ax, Bx, yb1 - 0.14, Bx, ybd + d_h/2 + G, c=CK)
    dia(ax, Bx, ybd, 1.1, d_h, "Added ≤ 14 d?", fs=5.5)

    # YES
    ybf = ybd - 0.42
    arr(ax, Bx, ybd - d_h/2, Bx, ybf + 0.09 + G, c=CT)
    lbl(ax, Bx + 0.10, ybd - d_h/2 - 0.06, "Yes", fs=5, fw="bold", c=CT)
    box(ax, Bx, ybf, 1.4, 0.18,
        "Follow strategy", fs=5.5, ec=CT, fc="#F0FAF8")

    # NO
    xno = 6.60
    ln(ax, Bx + 0.55, ybd, xno - 0.12, ybd)
    lbl(ax, Bx + 0.65, ybd + 0.08, "No", fs=5, fw="bold", c=CB)
    box(ax, xno + 0.10, ybd, 0.55, 0.22,
        "Censor\nd 14", fs=5, fw="bold", ec=CB, fc="#FFF5F0")

    # ── Footnote ─────────────────────────────────────────────────────────
    ax.text(0.05, 0.08,
            "PVTT = portal vein tumour thrombus; HVTT = hepatic vein tumour thrombus; "
            "IVC/RA = inferior vena cava / right atrium.",
            ha="left", va="bottom", fontsize=4.5, color=CG, fontstyle="italic")

    # ── Save ─────────────────────────────────────────────────────────────
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = out_dir / "tte_nlr_strategy_rules"
    fig.savefig(f"{stem}.pdf", facecolor="white")
    fig.savefig(f"{stem}.png", dpi=300, facecolor="white")
    plt.close(fig)
    print(f"Saved: {stem}.pdf / .png")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", type=str, default=str(DEFAULT_OUT))
    args = ap.parse_args()
    draw(Path(args.out_dir))


if __name__ == "__main__":
    main()
