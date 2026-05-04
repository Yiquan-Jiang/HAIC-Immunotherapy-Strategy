#!/usr/bin/env python3
"""
TTE Strategy Schematic v3 — Nature-level minimalist
Only 2 accent colours (navy + red-orange) + gray.
Strict grid, zero overlap, no decorative elements.
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Polygon
import numpy as np, os

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 7,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.02,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
})

# Only 4 colours
NAVY  = '#3C5488'
RED   = '#E64B35'
GRAY  = '#666666'
LTGR  = '#CCCCCC'
BG    = '#F8F8F8'
WHITE = '#FFFFFF'
BLACK = '#222222'
MIDGR = '#999999'


def box(ax, x, y, w, h, txt, fc=WHITE, ec=LTGR, fs=6.5, fw='normal',
        tc=BLACK, lw=0.6, zorder=4, pad=0.08):
    p = FancyBboxPatch((x-w/2, y-h/2), w, h,
        boxstyle=f'round,pad={pad}', fc=fc, ec=ec, lw=lw, zorder=zorder)
    ax.add_patch(p)
    ax.text(x, y, txt, ha='center', va='center', fontsize=fs,
            fontweight=fw, color=tc, zorder=zorder+1, linespacing=1.28)


def diam(ax, x, y, w, h, txt, ec=GRAY, fs=6, tc=BLACK):
    v = np.array([[x,y+h/2],[x+w/2,y],[x,y-h/2],[x-w/2,y],[x,y+h/2]])
    ax.add_patch(Polygon(v, closed=True, fc=WHITE, ec=ec, lw=0.6, zorder=4))
    ax.text(x, y, txt, ha='center', va='center', fontsize=fs,
            fontweight='bold', color=tc, zorder=5, linespacing=1.2)


def ar(ax, x1, y1, x2, y2, lbl='', c=GRAY, lw=0.6, fs=5.5,
       lo=(0,0), lc=None, lfw='normal'):
    lc = lc or MIDGR
    ax.annotate('', xy=(x2,y2), xytext=(x1,y1),
        arrowprops=dict(arrowstyle='->', color=c, lw=lw,
                        shrinkA=1, shrinkB=1), zorder=3)
    if lbl:
        ax.text((x1+x2)/2+lo[0], (y1+y2)/2+lo[1], lbl, ha='center',
                va='bottom', fontsize=fs, color=lc, fontweight=lfw, zorder=4)


# ═══════════════════════════════════════════════════════════════════════════
# Canvas: Nature double-col 7.2 in, height ~8 in
# Coordinate system: x 0-20, y 0-28, show y 7-27.5
# ═══════════════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(7.2, 8.0))
ax.set_xlim(0.5, 19.5)
ax.set_ylim(12.7, 27.2)
ax.axis('off')

M  = 10.0   # midpoint
LX = 4.8    # left col
RX = 14.2   # right col

# ═══════════════════════════════════════════════════════════════════════════
# TITLE
# ═══════════════════════════════════════════════════════════════════════════
ax.text(M, 27.0, 'Target Trial Emulation: Strategy Schematic',
        ha='center', fontsize=10, fontweight='bold', color=BLACK)
ax.text(M, 26.55, 'Clone-Censor-Weight framework  |  Outcome: overall survival',
        ha='center', fontsize=6.5, color=MIDGR, fontstyle='italic')

# ═══════════════════════════════════════════════════════════════════════════
# TIME ZERO
# ═══════════════════════════════════════════════════════════════════════════
box(ax, M, 25.9, 5.0, 0.55, 'Time zero: first HAIC date',
    fc=NAVY, ec=NAVY, fs=7.5, fw='bold', tc=WHITE, lw=0.8)

ar(ax, M, 25.62, M, 25.22, c=GRAY, lw=0.7)

box(ax, M, 24.95, 5.8, 0.42,
    'Each eligible patient cloned into both arms',
    fc=BG, ec=LTGR, fs=6.5, fw='bold', tc=NAVY)

# Split
ar(ax, M-1.5, 24.73, LX+1.0, 24.22, c=NAVY, lw=0.8)
ar(ax, M+1.5, 24.73, RX-1.0, 24.22, c=RED, lw=0.8)

# ═══════════════════════════════════════════════════════════════════════════
# ARM HEADERS
# ═══════════════════════════════════════════════════════════════════════════
box(ax, LX, 23.95, 7.0, 0.45, 'Strategy B: Early Combination',
    fc=NAVY, ec=NAVY, fs=7.5, fw='bold', tc=WHITE, lw=0.8)
box(ax, RX, 23.95, 7.0, 0.45, 'Strategy A: Dynamic',
    fc=RED, ec=RED, fs=7.5, fw='bold', tc=WHITE, lw=0.8)

# ═══════════════════════════════════════════════════════════════════════════
# LEFT ARM
# ═══════════════════════════════════════════════════════════════════════════
YL = {
    'rule': 23.05,
    'dec':  22.15,
    'yes':  21.50,
    'no':   21.50,
    'grace':21.05,
}

box(ax, LX, YL['rule'], 6.0, 0.55,
    'Add immunotherapy within 14 days\nof first HAIC',
    fs=6.5)

ar(ax, LX, YL['rule']-0.29, LX, YL['dec']+0.28, c=GRAY, lw=0.6)

diam(ax, LX, YL['dec'], 3.6, 0.60,
     'Immune added\n<  14 days?', ec=NAVY, fs=6)

# Yes
ar(ax, LX-1.8, YL['dec'], LX-3.3, YL['dec'],
   lbl='Yes', c=NAVY, lw=0.6, lc=NAVY, lfw='bold', lo=(0, 0.10))
box(ax, LX-3.3, YL['yes'], 2.0, 0.40,
    'Follow-up', fs=6, ec=NAVY)

# No
ar(ax, LX+1.8, YL['dec'], LX+3.3, YL['dec'],
   lbl='No', c=RED, lw=0.6, lc=RED, lfw='bold', lo=(0, 0.10))
box(ax, LX+3.3, YL['no'], 2.0, 0.40,
    'Censor at day 14', fs=6, ec=RED)

ax.text(LX, YL['grace'], 'Grace period = 14 days',
        ha='center', fontsize=6, color=NAVY, fontweight='bold',
        fontstyle='italic')

# ═══════════════════════════════════════════════════════════════════════════
# RIGHT ARM
# ═══════════════════════════════════════════════════════════════════════════
YR = {
    'trig_h':  23.30,
    'trig_sub':23.05,
    'trig':    21.90,
    'dec1':    20.55,
    'yes1':    20.55,
    'ex_h':    19.75,
    'ex_sub':  19.50,
    'ex':      18.95,
    'dec2':    18.10,
    'reeval':  17.65,
    'fb_h':    17.10,
    'fb':      16.60,
    'grace':   16.05,
}

# ── Trigger Rules ────────────────────────────────────────────────────────
ax.text(RX, YR['trig_h'], 'TRIGGER RULES',
        ha='center', fontsize=7, fontweight='bold', color=RED)
ax.text(RX, YR['trig_sub'],
        'Evaluated at pre-HAIC-3  ·  Absolute (exemptions cannot override)',
        ha='center', fontsize=5, color=MIDGR, fontstyle='italic')

trig = (
    'R1   PVTT Vp3/4 | distant mets | HVTT | IVC/RA | LN mets\n'
    'R2   Baseline NLR >= 2.5   or   pre-HAIC-3 NLR >= 2.5\n'
    'R3   Baseline max tumor diameter > 13 cm\n'
    'R4   Baseline PIVKA > 12 000 mAU/mL\n'
    'R5   Baseline AFP < 20 ng/mL'
)
box(ax, RX, YR['trig'], 8.2, 1.75, trig, fs=6, ec=GRAY, lw=0.5)

ar(ax, RX, YR['trig']-0.89, RX, YR['dec1']+0.27, c=GRAY, lw=0.6)

diam(ax, RX, YR['dec1'], 3.2, 0.55,
     'Any rule\ntriggered?', ec=RED, fs=6)

# Yes — box above the arrow endpoint for clarity
ar(ax, RX-1.6, YR['dec1'], RX-3.2, YR['dec1'],
   lbl='Yes', c=RED, lw=0.6, lc=RED, lfw='bold', lo=(0, 0.12))
box(ax, RX-3.2, YR['yes1']+0.30, 1.8, 0.38,
    'Add immune\nat pre-HAIC-3', fs=5.5, ec=RED)

# No
ar(ax, RX, YR['dec1']-0.29, RX, YR['ex_h']+0.10,
   lbl='No', c=GRAY, lw=0.6, lo=(0.22, 0.03), lfw='bold')

# ── Exemption Rules ──────────────────────────────────────────────────────
ax.text(RX, YR['ex_h'], 'EXEMPTION RULES',
        ha='center', fontsize=6.5, fontweight='bold', color=NAVY)
ax.text(RX, YR['ex_sub'],
        'Per-cycle evaluation from cycle >= 4',
        ha='center', fontsize=5, color=MIDGR, fontstyle='italic')

box(ax, RX, YR['ex'], 6.8, 0.55,
    'Ex1  AFP drop >= 90% vs baseline\n'
    'Ex2  Current NLR < 1.8',
    fs=6, ec=GRAY, lw=0.5)

ar(ax, RX, YR['ex']-0.29, RX, YR['dec2']+0.24, c=GRAY, lw=0.6)

diam(ax, RX, YR['dec2'], 3.0, 0.50,
     'Exempt at\ncycle k?', ec=NAVY, fs=6)

# Yes → right
ar(ax, RX+1.5, YR['dec2'], RX+2.8, YR['dec2'],
   lbl='Yes', c=NAVY, lw=0.6, lc=NAVY, lfw='bold', lo=(0, 0.10))
box(ax, RX+2.8, YR['reeval'], 1.6, 0.35,
    'Re-evaluate\nnext cycle', fs=5.5, ec=NAVY)

# Loop back
ax.annotate('', xy=(RX+2.8, YR['ex']+0.28),
            xytext=(RX+2.8, YR['reeval']+0.18),
            arrowprops=dict(arrowstyle='->', color=NAVY, lw=0.5,
                            connectionstyle='arc3,rad=-0.2',
                            shrinkA=1, shrinkB=1), zorder=2)

# Always exempt note (positioned below re-evaluate, inside bounds)
ax.text(RX+2.8, YR['reeval']-0.38,
        'always exempt\n= never triggered',
        ha='center', fontsize=4.5, color=NAVY, fontstyle='italic',
        linespacing=1.15)

# No → fallback
ar(ax, RX, YR['dec2']-0.27, RX, YR['fb_h']+0.10,
   lbl='No', c=GRAY, lw=0.6, lo=(0.22, 0.03), lfw='bold')

# ── Fallback ─────────────────────────────────────────────────────────────
ax.text(RX, YR['fb_h'], 'FALLBACK',
        ha='center', fontsize=6.5, fontweight='bold', color=RED)

box(ax, RX, YR['fb'], 5.6, 0.38,
    'Not exempt at cycle k  \u2192  add immune at cycle k',
    fs=6, ec=GRAY, lw=0.5)

ar(ax, RX, YR['fb']-0.21, RX, YR['grace']+0.10, c=GRAY, lw=0.6)

ax.text(RX, YR['grace'], 'Grace period = 90 days from trigger day',
        ha='center', fontsize=6, color=RED, fontweight='bold',
        fontstyle='italic')

# ═══════════════════════════════════════════════════════════════════════════
# THIN VERTICAL SEPARATOR
# ═══════════════════════════════════════════════════════════════════════════
ax.plot([M, M], [YR['grace']-0.2, 23.72], color=LTGR, lw=0.6, zorder=1)

# ═══════════════════════════════════════════════════════════════════════════
# BOTTOM — Censoring + Outcome
# ═══════════════════════════════════════════════════════════════════════════

# Light horizontal rule to separate
ax.plot([1.0, 19.0], [15.40, 15.40], color=LTGR, lw=0.4, zorder=1)

# Arrows down
ax.plot([LX, LX], [YL['grace']-0.18, 15.40], color=LTGR, lw=0.4,
        ls=':', zorder=1)
ax.plot([RX, RX], [YR['grace']-0.18, 15.40], color=LTGR, lw=0.4,
        ls=':', zorder=1)

YB_cen = 14.55
YB_out = 13.30

censor_txt = (
    'Artificial censoring (Dynamic arm)\n'
    'Case 1   Triggered at pre-HAIC-3, no immune within 90 d  \u2192  censor at trigger + 90 d\n'
    'Case 2   Fallback triggered, no immune within 90 d  \u2192  censor at fallback + 90 d\n'
    'Case 3   Never triggered but received immune  \u2192  censor at actual immune start'
)
box(ax, M, YB_cen, 16.5, 1.25, censor_txt,
    fs=5.8, ec=GRAY, lw=0.5, pad=0.12)

ar(ax, M, YB_cen-0.65, M, YB_out+0.25, c=GRAY, lw=0.7)

box(ax, M, YB_out, 9.5, 0.50,
    'Stabilized IPCW  \u2192  Weighted Cox PH  \u2192  HR + RMST',
    fc=NAVY, ec=NAVY, fs=7.5, fw='bold', tc=WHITE, lw=0.8)

# ═══════════════════════════════════════════════════════════════════════════
# TIMELINE (minimal, far left)
# ═══════════════════════════════════════════════════════════════════════════
tl = 0.9
ax.plot([tl, tl], [YR['fb'], 25.9], color=LTGR, lw=0.7, zorder=1)
ax.annotate('', xy=(tl, YR['fb']), xytext=(tl, YR['fb']+0.4),
            arrowprops=dict(arrowstyle='->', color=LTGR, lw=0.5), zorder=1)
for yy, lab in [(25.9, 'HAIC #1'), (YL['dec'], 'pre-HAIC-3'),
                (YR['ex'], 'Cycle >= 4'), (YR['fb'], 'Fallback')]:
    ax.plot(tl, yy, 'o', color=NAVY, ms=2, zorder=3)
    ax.text(tl-0.08, yy, lab, ha='right', va='center', fontsize=4.5,
            color=NAVY, fontweight='bold')

# ═══════════════════════════════════════════════════════════════════════════
# SAVE
# ═══════════════════════════════════════════════════════════════════════════
_script_dir = os.path.dirname(os.path.abspath(__file__))
out = os.path.normpath(os.path.join(_script_dir, '..', '..', 'output', 'schematic_figures'))
os.makedirs(out, exist_ok=True)
s = os.path.join(out, 'tte_strategy_schematic_v3')
fig.savefig(f'{s}.pdf', bbox_inches='tight', pad_inches=0.02)
fig.savefig(f'{s}.png', dpi=300, bbox_inches='tight', pad_inches=0.02)
plt.close(fig)
print(f"Saved -> {s}.pdf / .png")
