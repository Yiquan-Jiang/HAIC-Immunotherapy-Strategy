"""
Publication-quality figures from R core TTE analysis results.
Reads CSV outputs from tte_nlr_R_core.R (or cohort variant) and generates all figures.
Target: NEJM / Lancet / Nature Medicine level.

Usage:
  python tte_nlr_R_figures.py
      -> reads DATA_DIR/NLR_BASED_RULES_R/*.csv, writes figures there
  python tte_nlr_R_figures.py /path/to/NLR_BASED_RULES_R/cohort_7group_psm02
      -> reads/writes in that folder (e.g. after tte_nlr_R_core_cohort_7group_psm02.R)
  python tte_nlr_R_figures.py /path/to/PIV_BASED_RULES_R/cohort_7group_psm02
      -> same figures from tte_piv_R_core_cohort_7group_psm02.R outputs
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.patches as mpatches
import matplotlib.patheffects as pe
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

# ─── Paths ────────────────────────────────────────────────────────────────────
SCRIPT_DIR   = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.normpath(os.path.join(SCRIPT_DIR, '..', '..'))
if len(sys.argv) > 1 and str(sys.argv[1]).strip():
    R_DIR = os.path.abspath(str(sys.argv[1]).strip())
else:
    R_DIR = os.path.join(PROJECT_ROOT, 'output', 'step3_tte', 'NLR_BASED_RULES_R')
OUT_DIR = R_DIR
_rdir_norm = os.path.normpath(R_DIR).replace(os.sep, '/')
ANALYSIS_LABEL = (
    'PIV-Based' if 'PIV_BASED_RULES' in _rdir_norm.upper() else 'NLR-Based'
)

# ─── Color palette (aligned with plot_tte_schematic_v3.py) ───────────────────
COLOR_DYN   = '#E64B35'   # red-orange – Dynamic strategy  (Strategy A)
COLOR_EARLY = '#3C5488'   # navy       – Early combination (Strategy B)
COLOR_GRAY  = '#8491B4'
COLOR_GREEN = '#00A087'
GRAY_DARK   = '#333333'
GRAY_MED    = '#666666'
GRAY_LIGHT  = '#999999'
GRAY_VLIGHT = '#CCCCCC'

STATS_BOX = dict(boxstyle='round,pad=0.45', facecolor='white',
                 edgecolor=GRAY_VLIGHT, alpha=0.96, linewidth=0.6)

# ─── KM-specific rcParams ─────────────────────────────────────────────────────
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 10,
    'axes.labelsize': 11,
    'axes.titlesize': 12,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'axes.linewidth': 0.8,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'xtick.major.width': 0.8,
    'ytick.major.width': 0.8,
    'xtick.major.size': 3.5,
    'ytick.major.size': 3.5,
    'xtick.minor.size': 1.8,
    'ytick.minor.size': 1.8,
    'lines.linewidth': 2.0,
    'lines.markersize': 4,
    'legend.frameon': False,
    'legend.borderpad': 0.3,
    'legend.handlelength': 2.0,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.05,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
    'axes.grid': False,
})

RISK_TIMES = [0, 6, 12, 18, 24, 30, 36]
KM_MAX     = 36

# ─── Helpers ──────────────────────────────────────────────────────────────────
def save_fig(fig, name):
    """Save PDF (vector) + PNG (300 DPI) and print paths."""
    for ext in ('.pdf', '.png'):
        path = os.path.join(OUT_DIR, name + ext)
        fig.savefig(path, dpi=300, bbox_inches='tight', pad_inches=0.05)
        print(f"  Saved: {path}")
    plt.close(fig)


def fmt_p(p):
    if p < 0.001:
        return 'P < 0.001'
    return f'P = {p:.3f}'


def fmt_p_val(p):
    if p < 0.001:
        return '<0.001'
    return f'{p:.3f}'


# ─── Load data ────────────────────────────────────────────────────────────────
print("=" * 70)
print(f"TTE {ANALYSIS_LABEL} Rules: Publication Figures (from R results)")
print(f"R results directory: {R_DIR}")
print("=" * 70)

main_res = pd.read_csv(os.path.join(R_DIR, 'R_main_results.csv'))
rmst_df  = pd.read_csv(os.path.join(R_DIR, 'R_Table3_RMST.csv'))
sens_df  = pd.read_csv(os.path.join(R_DIR, 'R_Table4_sensitivity.csv'))
km_dyn   = pd.read_csv(os.path.join(R_DIR, 'km_dynamic.csv'))
km_early = pd.read_csv(os.path.join(R_DIR, 'km_early_combo.csv'))
risk_tbl = pd.read_csv(os.path.join(R_DIR, 'R_risk_table.csv'))
clone_df = pd.read_csv(os.path.join(R_DIR, 'R_clone_dataset.csv'))

hr          = main_res['HR'].iloc[0]
hr_lo       = main_res['HR_lo'].iloc[0]
hr_hi       = main_res['HR_hi'].iloc[0]
hr_p        = main_res['HR_p'].iloc[0]
n_dyn       = int(main_res['n_dyn'].iloc[0])
n_early     = int(main_res['n_early'].iloc[0])
n_eligible  = int(main_res['n_eligible'].iloc[0])
ess_dyn     = main_res['ess_dyn'].iloc[0]
ess_early   = main_res['ess_early'].iloc[0]
n_dyn_events   = int(main_res['events_dyn'].iloc[0])
n_early_events = int(main_res['events_early'].iloc[0])
n_dyn_cens     = int(main_res['censored_dyn'].iloc[0])
n_early_cens   = int(main_res['censored_early'].iloc[0])
med_dyn        = main_res['median_dyn'].iloc[0]
med_early      = main_res['median_early'].iloc[0]
ph_p           = main_res['ph_test_p'].iloc[0] if 'ph_test_p' in main_res.columns else np.nan

rmst24 = rmst_df[rmst_df['tau'] == 24].iloc[0]

# ══════════════════════════════════════════════════════════════════════════════
# Figure 1 — IPCW-weighted Kaplan-Meier survival curves
# ══════════════════════════════════════════════════════════════════════════════
print("\n--- Figure 1: KM survival curves ---")

fig1 = plt.figure(figsize=(5.0, 5.5))
gs1  = fig1.add_gridspec(2, 1, height_ratios=[4.2, 0.8], hspace=0.04)
ax_km = fig1.add_subplot(gs1[0])
ax_rt = fig1.add_subplot(gs1[1])

# ── KM curves ─────────────────────────────────────────────────────────────────
mask_d = km_dyn['time'] <= KM_MAX + 0.5
mask_e = km_early['time'] <= KM_MAX + 0.5

t_d = km_dyn['time'][mask_d].values
s_d = km_dyn['surv'][mask_d].values
t_e = km_early['time'][mask_e].values
s_e = km_early['surv'][mask_e].values

ax_km.step(t_d, s_d, where='post', color=COLOR_DYN, linewidth=2.0,
           linestyle='-',
           label=f'Dynamic strategy (n = {n_dyn:,})')
ax_km.step(t_e, s_e, where='post', color=COLOR_EARLY, linewidth=2.0,
           linestyle='--',
           label=f'Early combination (n = {n_early:,})')

# 95% CI shading — approximate via ±1.36/√n envelope (visual only)
# Use simple ±0.05 band as placeholder when CI columns absent
for t_arr, s_arr, col in [(t_d, s_d, COLOR_DYN), (t_e, s_e, COLOR_EARLY)]:
    ax_km.fill_between(t_arr,
                       np.clip(s_arr - 0.04, 0, 1),
                       np.clip(s_arr + 0.04, 0, 1),
                       step='post', alpha=0.10, color=col)

ax_km.set_xlim(0, KM_MAX + 0.5)
ax_km.set_ylim(-0.02, 1.05)
ax_km.set_ylabel('Overall Survival Probability', fontsize=11)
ax_km.set_title('IPCW-Weighted Kaplan–Meier Curves\n'
                '(Clone-Censor-Weighting, CCW)',
                fontweight='bold', pad=10, fontsize=11)
ax_km.set_xticks(RISK_TIMES)
ax_km.yaxis.set_major_locator(mticker.MultipleLocator(0.2))
ax_km.yaxis.set_minor_locator(mticker.MultipleLocator(0.1))
ax_km.tick_params(axis='both', which='both', direction='out')

# ── Legend ────────────────────────────────────────────────────────────────────
leg = ax_km.legend(loc='upper right', handlelength=2.2, fontsize=9,
                   frameon=False, borderpad=0.3)

# ── Statistics box (bottom-right) ─────────────────────────────────────────────
ph_str = f'PH test: {fmt_p(ph_p)}' if not np.isnan(ph_p) else ''
ann_lines = [
    f'HR = {hr:.2f} (95% CI {hr_lo:.2f}–{hr_hi:.2f})',
    fmt_p(hr_p),
    f'ΔRMST (24 mo) = {rmst24["delta"]:+.2f} mo '
    f'({rmst24["ci_lo"]:+.2f}, {rmst24["ci_hi"]:+.2f})',
]
if ph_str:
    ann_lines.append(ph_str)
ax_km.text(0.97, 0.05, '\n'.join(ann_lines),
           transform=ax_km.transAxes,
           fontsize=8.5, va='bottom', ha='right', bbox=STATS_BOX,
           linespacing=1.55)

# ── Risk table ────────────────────────────────────────────────────────────────
risk_d = risk_tbl[risk_tbl['arm'] == 'dynamic'].copy()
risk_e = risk_tbl[risk_tbl['arm'] == 'early_combo'].copy()

n_groups = 2
ax_rt.set_xlim(0, KM_MAX + 0.5)
ax_rt.set_ylim(-0.8, n_groups + 0.2)
ax_rt.axis('off')

ax_rt.text(KM_MAX / 2, n_groups - 0.05, 'No. at risk',
           ha='center', va='top', fontsize=9,
           color=GRAY_MED, style='italic')

for t_pt in RISK_TIMES:
    row_d = risk_d[risk_d['time'] == t_pt]
    row_e = risk_e[risk_e['time'] == t_pt]
    n_d = int(row_d['n_at_risk'].iloc[0]) if len(row_d) else 0
    n_e = int(row_e['n_at_risk'].iloc[0]) if len(row_e) else 0
    ax_rt.text(t_pt, 1.1, str(n_d),
               ha='center', va='center', fontsize=9,
               color=COLOR_DYN, fontweight='bold')
    ax_rt.text(t_pt, 0.1, str(n_e),
               ha='center', va='center', fontsize=9,
               color=COLOR_EARLY, fontweight='bold')

ax_rt.text(KM_MAX / 2, -0.65, 'Time from First HAIC (months)',
           ha='center', va='center', fontsize=11, color=GRAY_DARK)

save_fig(fig1, 'R_Fig1_KM_full')

# ══════════════════════════════════════════════════════════════════════════════
# Figure 2 — ΔRMST forest plot
# ══════════════════════════════════════════════════════════════════════════════
print("--- Figure 2: ΔRMST forest plot ---")

labels = [f'ΔRMST  (τ = {int(r["tau"])} mo)' for _, r in rmst_df.iterrows()]
ests   = rmst_df['delta'].tolist()
lo_ci  = rmst_df['ci_lo'].tolist()
hi_ci  = rmst_df['ci_hi'].tolist()
pvals  = rmst_df['p'].tolist()
n_items = len(labels)

fig2, ax_f = plt.subplots(figsize=(5.5, 0.70 * n_items + 2.2))
y_pos = np.arange(n_items)[::-1].astype(float)

for i in range(n_items):
    yi  = y_pos[i]
    est = ests[i]; cli = lo_ci[i]; chi = hi_ci[i]; p = pvals[i]
    col = COLOR_DYN if est >= 0 else COLOR_EARLY

    # CI line with round caps
    ax_f.plot([cli, chi], [yi, yi], color=col, linewidth=2.0,
              solid_capstyle='round', zorder=3)
    # Diamond marker
    ax_f.scatter([est], [yi], marker='D', s=80, color=col, zorder=5,
                 edgecolors='white', linewidths=0.8)

    # Right-side annotation
    p_str = fmt_p_val(p)
    ax_f.text(chi + 0.15, yi,
              f'{est:+.2f} ({cli:+.2f}, {chi:+.2f})   P = {p_str}',
              va='center', fontsize=8.5, color=GRAY_DARK)

# Reference line
ax_f.axvline(x=0, color='#555555', linestyle='--', linewidth=0.9, zorder=1)

ax_f.set_yticks(y_pos)
ax_f.set_yticklabels(labels, fontsize=9.5)
ax_f.set_xlabel('ΔRMST: Dynamic − Early Combination (months)', fontsize=10)

p_hr_fmt = fmt_p_val(hr_p)
ax_f.set_title(
    f'HR = {hr:.2f} (95% CI {hr_lo:.2f}–{hr_hi:.2f}),  {fmt_p(hr_p)}',
    fontsize=9.5, fontweight='bold', pad=10)

# x-axis limits: leave room for text
xl = min(lo_ci) - 0.5
xr = max(hi_ci) + max(3.5, (max(hi_ci) - min(lo_ci)) * 0.85)
ax_f.set_xlim(xl, xr)
ax_f.set_ylim(-0.8, n_items - 0.2)

# Direction annotation
ax_f.annotate('', xy=(xl + 0.05, -0.65), xytext=(xl + 1.0, -0.65),
              xycoords='data', textcoords='data',
              arrowprops=dict(arrowstyle='->', color=COLOR_EARLY, lw=1.2))
ax_f.text(xl + 1.1, -0.65, 'Favors early combination',
          va='center', fontsize=8, color=COLOR_EARLY)
ax_f.annotate('', xy=(0 + 1.0, -0.65), xytext=(0 + 0.05, -0.65),
              xycoords='data', textcoords='data',
              arrowprops=dict(arrowstyle='->', color=COLOR_DYN, lw=1.2))
ax_f.text(0 + 1.1, -0.65, 'Favors dynamic strategy',
          va='center', fontsize=8, color=COLOR_DYN)

ax_f.spines['left'].set_visible(False)
ax_f.tick_params(left=False)

save_fig(fig2, 'R_Fig2_forest')

# ══════════════════════════════════════════════════════════════════════════════
# Figure 3 — Study flow diagram
# ══════════════════════════════════════════════════════════════════════════════
print("--- Figure 3: Flow diagram ---")

N_DB = 1614
n_excluded = N_DB - n_eligible
pct_d = n_dyn_cens / max(n_dyn, 1) * 100
pct_e = n_early_cens / max(n_early, 1) * 100

fig3, ax = plt.subplots(figsize=(7.0, 7.0))
ax.set_xlim(0, 10)
ax.set_ylim(0, 10)
ax.axis('off')

# ── Box styles ────────────────────────────────────────────────────────────────
def draw_box(ax, x, y, w, h, text, fc, ec, fontsize=9, bold=False, italic=False):
    rect = FancyBboxPatch((x - w/2, y - h/2), w, h,
                          boxstyle='round,pad=0.12',
                          facecolor=fc, edgecolor=ec, linewidth=1.1, zorder=3)
    ax.add_patch(rect)
    fw = 'bold' if bold else 'normal'
    fs = 'italic' if italic else 'normal'
    ax.text(x, y, text, ha='center', va='center',
            fontsize=fontsize, fontweight=fw, fontstyle=fs,
            zorder=4, linespacing=1.5)


def draw_arrow(ax, x1, y1, x2, y2):
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle='->', color=GRAY_DARK,
                                lw=1.3, connectionstyle='arc3,rad=0.0'))


def draw_excl_box(ax, x, y, w, h, text):
    rect = FancyBboxPatch((x - w/2, y - h/2), w, h,
                          boxstyle='round,pad=0.12',
                          facecolor='#FDE8E8', edgecolor='#C62828',
                          linewidth=0.9, zorder=3)
    ax.add_patch(rect)
    ax.text(x, y, text, ha='center', va='center',
            fontsize=8, zorder=4, linespacing=1.5, color='#B71C1C')


# ── Boxes ─────────────────────────────────────────────────────────────────────
# 1. Total database
draw_box(ax, 5, 9.2, 3.2, 0.75,
         f'Total Database\nn = {N_DB:,}',
         '#E3F2FD', '#1565C0', fontsize=9.5, bold=True)

# Arrow down
draw_arrow(ax, 5, 8.82, 5, 8.25)

# Exclusion box (right)
draw_excl_box(ax, 8.1, 8.55, 2.8, 0.70,
              f'Excluded (OS ≤ 0 or missing):\nn = {n_excluded:,}')
ax.annotate('', xy=(6.75, 8.55), xytext=(5.3, 8.55),
            arrowprops=dict(arrowstyle='->', color='#C62828', lw=1.0))

# 2. Eligible
draw_box(ax, 5, 7.85, 3.2, 0.72,
         f'Eligible Population\nn = {n_eligible:,}',
         '#E3F2FD', '#1565C0', fontsize=9.5, bold=True)

draw_arrow(ax, 5, 7.49, 5, 7.05)

# Clone note
ax.text(5, 7.27, 'Cloned into both arms (CCW)',
        ha='center', va='center', fontsize=8.5,
        style='italic', color=GRAY_MED)

# Split arrows
ax.annotate('', xy=(2.5, 6.35), xytext=(4.2, 6.85),
            arrowprops=dict(arrowstyle='->', color=COLOR_DYN, lw=1.3))
ax.annotate('', xy=(7.5, 6.35), xytext=(5.8, 6.85),
            arrowprops=dict(arrowstyle='->', color=COLOR_EARLY, lw=1.3))

# 3a. Dynamic arm
draw_box(ax, 2.5, 5.95, 3.2, 0.72,
         f'Dynamic Strategy Arm\nn = {n_dyn:,}',
         '#FDE8E8', COLOR_DYN, fontsize=9.5, bold=True)

# 3b. Early combo arm
draw_box(ax, 7.5, 5.95, 3.2, 0.72,
         f'Early Combination Arm\nn = {n_early:,}',
         '#E8EAF6', COLOR_EARLY, fontsize=9.5, bold=True)

draw_arrow(ax, 2.5, 5.59, 2.5, 4.85)
draw_arrow(ax, 7.5, 5.59, 7.5, 4.85)

# 4a. Dynamic outcomes
draw_box(ax, 2.5, 4.35, 3.2, 0.95,
         f'Artificially censored: {n_dyn_cens:,} ({pct_d:.1f}%)\n'
         f'Events (deaths): {n_dyn_events:,}\n'
         f'Median follow-up: {med_dyn:.1f} mo',
         '#FDE8E8', COLOR_DYN, fontsize=8.5)

# 4b. Early combo outcomes
draw_box(ax, 7.5, 4.35, 3.2, 0.95,
         f'Artificially censored: {n_early_cens:,} ({pct_e:.1f}%)\n'
         f'Events (deaths): {n_early_events:,}\n'
         f'Median follow-up: {med_early:.1f} mo',
         '#E8EAF6', COLOR_EARLY, fontsize=8.5)

# Merge arrows
ax.annotate('', xy=(5, 2.90), xytext=(2.5, 3.87),
            arrowprops=dict(arrowstyle='->', color=GRAY_DARK, lw=1.2))
ax.annotate('', xy=(5, 2.90), xytext=(7.5, 3.87),
            arrowprops=dict(arrowstyle='->', color=GRAY_DARK, lw=1.2))

# 5. Analysis result
draw_box(ax, 5, 2.35, 6.5, 0.95,
         f'IPCW-Weighted Cox Analysis  (R survival::coxph, robust SE)\n'
         f'HR = {hr:.2f} (95% CI {hr_lo:.2f}–{hr_hi:.2f}),  {fmt_p(hr_p)}\n'
         f'ESS: Dynamic = {ess_dyn:.0f},  Early combination = {ess_early:.0f}',
         '#E8F5E9', '#2E7D32', fontsize=8.5, bold=False)

# Footer note
ax.text(5, 1.25,
        'Time zero: first HAIC  |  Grace period: early 14 d / dynamic 90 d  |  '
        'IPCW truncation: 99th percentile\n'
        f'Bootstrap: 500 paired iterations  |  RMST τ: 12, 18, 24, 36 months',
        ha='center', va='center', fontsize=7.5,
        color=GRAY_MED, style='italic', linespacing=1.6)

# Title
ax.text(5, 9.82, 'Study Flow Diagram',
        ha='center', va='center', fontsize=12,
        fontweight='bold', color=GRAY_DARK)
ax.text(5, 9.58, f'({ANALYSIS_LABEL} dynamic rules, CCW)',
        ha='center', va='center', fontsize=8.5,
        color=GRAY_MED, style='italic')

save_fig(fig3, 'R_Fig3_flowchart')

# ══════════════════════════════════════════════════════════════════════════════
# eFigure — IPCW weight distribution diagnostics (Panel A / B)
# ══════════════════════════════════════════════════════════════════════════════
print("--- eFigure: Weight diagnostics ---")

fig_d, axes = plt.subplots(1, 2, figsize=(7.0, 3.0))
plt.subplots_adjust(wspace=0.38)

panel_label_kw = dict(fontsize=10, fontweight='bold', va='top',
                      transform=None)

for ax_idx, (col_name, title, panel) in enumerate([
    ('sw_ipcw', 'IPCW Weight Distribution\n(Truncated at 99th Percentile)', 'A'),
    ('sw',      'Primary Weight Distribution\n(IPCW-Only)', 'B'),
]):
    ax = axes[ax_idx]
    for arm_name, color, lbl in [
        ('dynamic',     COLOR_DYN,   'Dynamic strategy'),
        ('early_combo', COLOR_EARLY, 'Early combination'),
    ]:
        ws = clone_df.loc[clone_df['arm'] == arm_name, col_name].values
        ws = ws[np.isfinite(ws)]
        ax.hist(ws, bins=60, alpha=0.55, color=color, label=lbl,
                density=True, edgecolor='white', linewidth=0.3)
        # KDE overlay
        from scipy.stats import gaussian_kde
        kde = gaussian_kde(ws, bw_method=0.25)
        x_grid = np.linspace(ws.min(), ws.max(), 300)
        ax.plot(x_grid, kde(x_grid), color=color, linewidth=1.5, alpha=0.9)

    ax.set_xlabel('Weight', fontsize=10)
    ax.set_ylabel('Density', fontsize=10)
    ax.set_title(title, fontsize=9.5, fontweight='bold', pad=8)
    ax.legend(fontsize=8, frameon=False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.text(-0.12, 1.08, panel, transform=ax.transAxes,
            fontsize=10, fontweight='bold', va='top')

    if col_name == 'sw':
        ax.text(0.97, 0.95,
                f'ESS: Dyn = {ess_dyn:.0f}\nESS: Early = {ess_early:.0f}',
                transform=ax.transAxes, fontsize=8,
                va='top', ha='right', bbox=STATS_BOX)

save_fig(fig_d, 'R_eFig_diagnostics')

# ══════════════════════════════════════════════════════════════════════════════
# eFigure — Sensitivity analyses forest plot
# ══════════════════════════════════════════════════════════════════════════════
print("--- eFigure: Sensitivity forest plot ---")

sens_valid = sens_df.dropna(subset=['HR']).copy().reset_index(drop=True)
n_s = len(sens_valid)

fig_s, ax_s = plt.subplots(figsize=(6.0, 0.65 * n_s + 2.0))
y_s = np.arange(n_s)[::-1].astype(float)

for i, row in sens_valid.iterrows():
    yi  = y_s[i]
    est = row['HR']; cli = row['CI_lo']; chi = row['CI_hi']
    p   = row['P']
    # Color: blue if HR < 1 (favors dynamic), orange if HR ≥ 1
    col = COLOR_DYN if est < 1 else COLOR_EARLY
    is_primary = (i == 0)

    lw = 2.2 if is_primary else 1.8
    ms = 90  if is_primary else 65

    ax_s.plot([cli, chi], [yi, yi], color=col, linewidth=lw,
              solid_capstyle='round', zorder=3,
              alpha=1.0 if is_primary else 0.75)
    ax_s.scatter([est], [yi], marker='D', s=ms, color=col, zorder=5,
                 edgecolors='white', linewidths=0.8,
                 alpha=1.0 if is_primary else 0.80)

    p_str = fmt_p_val(p)
    ax_s.text(chi + 0.03, yi,
              f'{est:.2f} ({cli:.2f}, {chi:.2f})   P = {p_str}',
              va='center', fontsize=8.5, color=GRAY_DARK)

# Reference line HR = 1
ax_s.axvline(x=1.0, color='#555555', linestyle='--', linewidth=0.9, zorder=1)

ax_s.set_yticks(y_s)
ax_s.set_yticklabels(sens_valid['analysis'].tolist(), fontsize=9)
ax_s.set_xlabel('Hazard Ratio (Dynamic vs Early Combination)', fontsize=10)
ax_s.set_title('Sensitivity Analyses', fontsize=11, fontweight='bold', pad=10)

# Direction arrows
xl_s = sens_valid['CI_lo'].min() - 0.15
xr_s = sens_valid['CI_hi'].max() + max(2.0, (sens_valid['CI_hi'].max()
                                              - sens_valid['CI_lo'].min()) * 0.7)
ax_s.set_xlim(xl_s, xr_s)
ax_s.set_ylim(-0.8, n_s - 0.2)

ax_s.text(0.5, -0.12,
          '← Favors dynamic strategy          Favors early combination →',
          transform=ax_s.transAxes, ha='center', fontsize=8.5, color=GRAY_MED)

ax_s.spines['left'].set_visible(False)
ax_s.tick_params(left=False)

# Highlight primary row
ax_s.axhspan(y_s[0] - 0.45, y_s[0] + 0.45,
             color='#F5F5F5', zorder=0, alpha=0.7)
ax_s.text(xl_s + 0.02, y_s[0] + 0.5, 'Primary',
          fontsize=7.5, color=GRAY_MED, style='italic', va='bottom')

save_fig(fig_s, 'R_eFig_sensitivity_forest')

# ─── Summary ──────────────────────────────────────────────────────────────────
print(f"\n{'='*70}")
print("All figures generated successfully.")
print(f"Output directory: {OUT_DIR}")
print(f"  R_Fig1_KM_full        .pdf / .png")
print(f"  R_Fig2_forest         .pdf / .png")
print(f"  R_Fig3_flowchart      .pdf / .png")
print(f"  R_eFig_diagnostics    .pdf / .png")
print(f"  R_eFig_sensitivity_forest .pdf / .png")
print(f"{'='*70}")
