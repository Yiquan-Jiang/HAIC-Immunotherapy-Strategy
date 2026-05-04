"""
Publication-quality figures for the IT-Rules TTE applied to TWO cohorts.
Reads CSV outputs from tte_IT_R_two_cohorts.R and generates figures per cohort.
Uses the display label "Adaptive On Demand" (instead of "Dynamic").

Usage:
  # Single cohort folder
  python tte_IT_R_figures_two_cohorts.py <cohort_output_folder>

  # Or no arg: process both cohorts under IT_RULES_R_two_cohorts/
  python tte_IT_R_figures_two_cohorts.py
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

SCRIPT_DIR   = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.normpath(os.path.join(SCRIPT_DIR, '..', '..'))
DEFAULT_BASE = os.path.join(PROJECT_ROOT, 'output', 'step3_tte',
                            'IT_RULES_R_two_cohorts')

COLOR_DYN   = '#E64B35'
COLOR_EARLY = '#3C5488'
COLOR_GRAY  = '#8491B4'
COLOR_GREEN = '#00A087'
GRAY_DARK   = '#333333'
GRAY_MED    = '#666666'
GRAY_LIGHT  = '#999999'
GRAY_VLIGHT = '#CCCCCC'

STATS_BOX = dict(boxstyle='round,pad=0.45', facecolor='white',
                 edgecolor=GRAY_VLIGHT, alpha=0.96, linewidth=0.6)

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

LABEL_DYN         = 'Adaptive On Demand'
LABEL_EARLY       = 'Early Combination'
LABEL_DYN_SHORT   = 'AoD'
LABEL_EARLY_SHORT = 'Early combo'

COHORT_TITLES = {
    'cohort_3matched'     : 'HAIC then I+T on-demand (IT_RULES_v2: +SII +PIVKA +LN)',
    'cohort_7group_psm02' : 'HAIC then I on-demand (IT_RULES_v2: +NLR +PIVKA +LN)',
}


def fmt_p(p):
    if pd.isna(p):
        return ''
    if p < 0.001:
        return 'P < 0.001'
    return f'P = {p:.3f}'


def fmt_p_val(p):
    if pd.isna(p):
        return 'NA'
    if p < 0.001:
        return '<0.001'
    return f'{p:.3f}'


def save_fig(fig, out_dir, name):
    for ext in ('.pdf', '.png'):
        path = os.path.join(out_dir, name + ext)
        fig.savefig(path, dpi=300, bbox_inches='tight', pad_inches=0.05)
        print(f"  Saved: {path}")
    plt.close(fig)


def render_cohort(r_dir):
    """Generate all figures for a single cohort output folder."""
    if not os.path.isdir(r_dir):
        print(f"  [SKIP] directory does not exist: {r_dir}")
        return

    cohort_name = os.path.basename(r_dir.rstrip(os.sep))
    cohort_title = COHORT_TITLES.get(cohort_name, cohort_name)
    print("=" * 70)
    print(f"Cohort: {cohort_name}  ({cohort_title})")
    print(f"Reading from: {r_dir}")
    print("=" * 70)

    main_res = pd.read_csv(os.path.join(r_dir, 'R_main_results.csv'))
    rmst_df  = pd.read_csv(os.path.join(r_dir, 'R_Table3_RMST.csv'))
    sens_df  = pd.read_csv(os.path.join(r_dir, 'R_Table4_sensitivity.csv'))
    km_dyn   = pd.read_csv(os.path.join(r_dir, 'km_dynamic.csv'))
    km_early = pd.read_csv(os.path.join(r_dir, 'km_early_combo.csv'))
    risk_tbl = pd.read_csv(os.path.join(r_dir, 'R_risk_table.csv'))
    clone_df = pd.read_csv(os.path.join(r_dir, 'R_clone_dataset.csv'))

    hr         = main_res['HR'].iloc[0]
    hr_lo      = main_res['HR_lo'].iloc[0]
    hr_hi      = main_res['HR_hi'].iloc[0]
    hr_p       = main_res['HR_p'].iloc[0]
    n_dyn      = int(main_res['n_dyn'].iloc[0])
    n_early    = int(main_res['n_early'].iloc[0])
    ess_dyn    = main_res['ess_dyn'].iloc[0]
    ess_early  = main_res['ess_early'].iloc[0]
    ph_p = main_res['ph_test_p'].iloc[0] if 'ph_test_p' in main_res.columns else np.nan

    rmst24 = rmst_df[rmst_df['tau'] == 24].iloc[0]

    # ── Figure 1: KM curves ──────────────────────────────────────────────────
    print("--- Figure 1: KM survival curves ---")
    fig1 = plt.figure(figsize=(6.8, 5.9))
    gs1  = fig1.add_gridspec(2, 2, height_ratios=[4.2, 1.1],
                             width_ratios=[2.4, 10.0], hspace=0.12, wspace=0.05)
    ax_km = fig1.add_subplot(gs1[0, 1])
    ax_lb = fig1.add_subplot(gs1[1, 0])
    ax_rt = fig1.add_subplot(gs1[1, 1])

    mask_d = km_dyn['time']   <= KM_MAX + 0.5
    mask_e = km_early['time'] <= KM_MAX + 0.5
    t_d = km_dyn['time'][mask_d].values;  s_d = km_dyn['surv'][mask_d].values
    t_e = km_early['time'][mask_e].values; s_e = km_early['surv'][mask_e].values

    ax_km.step(t_d, s_d, where='post', color=COLOR_DYN, linewidth=2.0,
               linestyle='-',  label=f'{LABEL_DYN} (n = {n_dyn:,})')
    ax_km.step(t_e, s_e, where='post', color=COLOR_EARLY, linewidth=2.0,
               linestyle='--', label=f'{LABEL_EARLY} (n = {n_early:,})')

    for t_arr, s_arr, col in [(t_d, s_d, COLOR_DYN), (t_e, s_e, COLOR_EARLY)]:
        ax_km.fill_between(t_arr,
                           np.clip(s_arr - 0.04, 0, 1),
                           np.clip(s_arr + 0.04, 0, 1),
                           step='post', alpha=0.10, color=col)

    ax_km.set_xlim(0, KM_MAX + 0.5)
    ax_km.set_ylim(-0.02, 1.05)
    ax_km.set_ylabel('Overall Survival Probability', fontsize=11)
    ax_km.set_title(f'IPCW-Weighted Kaplan-Meier Curves\n'
                    f'{cohort_title}',
                    fontweight='bold', pad=10, fontsize=11)
    ax_km.set_xticks(RISK_TIMES)
    ax_km.yaxis.set_major_locator(mticker.MultipleLocator(0.2))
    ax_km.yaxis.set_minor_locator(mticker.MultipleLocator(0.1))
    ax_km.tick_params(axis='both', which='both', direction='out')
    ax_km.legend(loc='upper right', handlelength=2.2, fontsize=9,
                 frameon=False, borderpad=0.3)

    ph_str = f'PH test: {fmt_p(ph_p)}' if not pd.isna(ph_p) else ''
    ann_lines = [
        f'HR = {hr:.2f} (95% CI {hr_lo:.2f}-{hr_hi:.2f})',
        fmt_p(hr_p),
        f'delta RMST (24 mo) = {rmst24["delta"]:+.2f} mo '
        f'({rmst24["ci_lo"]:+.2f}, {rmst24["ci_hi"]:+.2f})',
    ]
    if ph_str:
        ann_lines.append(ph_str)
    ax_km.text(0.97, 0.05, '\n'.join(ann_lines),
               transform=ax_km.transAxes,
               fontsize=8.5, va='bottom', ha='right', bbox=STATS_BOX,
               linespacing=1.55)

    risk_d = risk_tbl[risk_tbl['arm'] == 'dynamic'].copy()
    risk_e = risk_tbl[risk_tbl['arm'] == 'early_combo'].copy()

    ax_lb.set_xlim(0, 1); ax_lb.set_ylim(0, 3.0); ax_lb.axis('off')
    y_dyn, y_early = 1.7, 0.85
    ax_lb.text(0.98, 2.55, 'No. at risk',
               ha='right', va='center', fontsize=9,
               color=GRAY_DARK, fontweight='bold')
    ax_lb.text(0.98, y_dyn, LABEL_DYN_SHORT,
               ha='right', va='center', fontsize=8.5,
               color=COLOR_DYN, fontweight='bold')
    ax_lb.text(0.98, y_early, LABEL_EARLY_SHORT,
               ha='right', va='center', fontsize=8.5,
               color=COLOR_EARLY, fontweight='bold')

    ax_rt.set_xlim(0, KM_MAX + 0.5); ax_rt.set_ylim(0, 3.0); ax_rt.axis('off')
    for t_pt in RISK_TIMES:
        row_d = risk_d[risk_d['time'] == t_pt]
        row_e = risk_e[risk_e['time'] == t_pt]
        n_d = int(row_d['n_at_risk'].iloc[0]) if len(row_d) else 0
        n_e = int(row_e['n_at_risk'].iloc[0]) if len(row_e) else 0
        ax_rt.text(t_pt, y_dyn,   f'{n_d:,}',
                   ha='center', va='center', fontsize=9,
                   color=COLOR_DYN, fontweight='bold')
        ax_rt.text(t_pt, y_early, f'{n_e:,}',
                   ha='center', va='center', fontsize=9,
                   color=COLOR_EARLY, fontweight='bold')

    ax_rt.text(KM_MAX / 2, 0.05, 'Time from First HAIC (months)',
               ha='center', va='center', fontsize=11, color=GRAY_DARK)
    save_fig(fig1, r_dir, 'R_Fig1_KM_full')

    # ── Figure 2: delta RMST forest plot ─────────────────────────────────────
    print("--- Figure 2: delta RMST forest plot ---")
    labels = [f'delta RMST  (tau = {int(r["tau"])} mo)' for _, r in rmst_df.iterrows()]
    ests  = rmst_df['delta'].tolist()
    lo_ci = rmst_df['ci_lo'].tolist()
    hi_ci = rmst_df['ci_hi'].tolist()
    pvals = rmst_df['p'].tolist()
    n_items = len(labels)

    fig2, ax_f = plt.subplots(figsize=(5.5, 0.70 * n_items + 2.2))
    y_pos = np.arange(n_items)[::-1].astype(float)

    for i in range(n_items):
        yi  = y_pos[i]
        est = ests[i]; cli = lo_ci[i]; chi = hi_ci[i]; p = pvals[i]
        col = COLOR_DYN if est >= 0 else COLOR_EARLY
        ax_f.plot([cli, chi], [yi, yi], color=col, linewidth=2.0,
                  solid_capstyle='round', zorder=3)
        ax_f.scatter([est], [yi], marker='D', s=80, color=col, zorder=5,
                     edgecolors='white', linewidths=0.8)
        p_str = fmt_p_val(p)
        ax_f.text(chi + 0.15, yi,
                  f'{est:+.2f} ({cli:+.2f}, {chi:+.2f})   P = {p_str}',
                  va='center', fontsize=8.5, color=GRAY_DARK)

    ax_f.axvline(x=0, color='#555555', linestyle='--', linewidth=0.9, zorder=1)
    ax_f.set_yticks(y_pos)
    ax_f.set_yticklabels(labels, fontsize=9.5)
    ax_f.set_xlabel(f'delta RMST: {LABEL_DYN} - {LABEL_EARLY} (months)', fontsize=10)
    ax_f.set_title(
        f'HR = {hr:.2f} (95% CI {hr_lo:.2f}-{hr_hi:.2f}),  {fmt_p(hr_p)}',
        fontsize=9.5, fontweight='bold', pad=10)

    xl = min(lo_ci) - 0.5
    xr = max(hi_ci) + max(3.5, (max(hi_ci) - min(lo_ci)) * 0.85)
    ax_f.set_xlim(xl, xr)
    ax_f.set_ylim(-1.2, n_items - 0.2)

    ax_f.text(0.02, -0.18, f'\u2190 Favors {LABEL_EARLY.lower()}',
              transform=ax_f.transAxes, ha='left', va='center',
              fontsize=8.5, color=COLOR_EARLY)
    ax_f.text(0.98, -0.18, f'Favors {LABEL_DYN.lower()} \u2192',
              transform=ax_f.transAxes, ha='right', va='center',
              fontsize=8.5, color=COLOR_DYN)
    ax_f.spines['left'].set_visible(False)
    ax_f.tick_params(left=False)
    save_fig(fig2, r_dir, 'R_Fig2_forest')

    # ── eFigure: Weight diagnostics ──────────────────────────────────────────
    print("--- eFigure: Weight diagnostics ---")
    fig_d, axes = plt.subplots(1, 2, figsize=(7.0, 3.0))
    plt.subplots_adjust(wspace=0.38)

    for ax_idx, (col_name, title, panel) in enumerate([
        ('sw_ipcw', 'IPCW Weight Distribution\n(Truncated at 99th Percentile)', 'A'),
        ('sw',      'Primary Weight Distribution\n(IPCW-Only)', 'B'),
    ]):
        ax = axes[ax_idx]
        for arm_name, color, lbl in [
            ('dynamic',     COLOR_DYN,   LABEL_DYN),
            ('early_combo', COLOR_EARLY, LABEL_EARLY),
        ]:
            ws = clone_df.loc[clone_df['arm'] == arm_name, col_name].values
            ws = ws[np.isfinite(ws)]
            if len(ws) == 0:
                continue
            ax.hist(ws, bins=60, alpha=0.55, color=color, label=lbl,
                    density=True, edgecolor='white', linewidth=0.3)
            try:
                from scipy.stats import gaussian_kde
                kde = gaussian_kde(ws, bw_method=0.25)
                x_grid = np.linspace(ws.min(), ws.max(), 300)
                ax.plot(x_grid, kde(x_grid), color=color, linewidth=1.5, alpha=0.9)
            except Exception:
                pass

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
                    f'ESS: {LABEL_DYN_SHORT} = {ess_dyn:.0f}\n'
                    f'ESS: {LABEL_EARLY_SHORT} = {ess_early:.0f}',
                    transform=ax.transAxes, fontsize=8,
                    va='top', ha='right', bbox=STATS_BOX)
    save_fig(fig_d, r_dir, 'R_eFig_diagnostics')

    # ── eFigure: Sensitivity forest ──────────────────────────────────────────
    print("--- eFigure: Sensitivity forest plot ---")
    sens_valid = sens_df.dropna(subset=['HR']).copy().reset_index(drop=True)
    n_s = len(sens_valid)

    fig_s, ax_s = plt.subplots(figsize=(6.5, 0.65 * n_s + 2.0))
    y_s = np.arange(n_s)[::-1].astype(float)

    for i, row in sens_valid.iterrows():
        yi  = y_s[i]
        est = row['HR']; cli = row['CI_lo']; chi = row['CI_hi']; p = row['P']
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

    ax_s.axvline(x=1.0, color='#555555', linestyle='--', linewidth=0.9, zorder=1)
    ax_s.set_yticks(y_s)
    ax_s.set_yticklabels(sens_valid['analysis'].tolist(), fontsize=9)
    ax_s.set_xlabel(f'Hazard Ratio ({LABEL_DYN} vs {LABEL_EARLY})', fontsize=10)
    ax_s.set_title('Sensitivity Analyses', fontsize=11, fontweight='bold', pad=10)

    xl_s = sens_valid['CI_lo'].min() - 0.15
    xr_s = sens_valid['CI_hi'].max() + max(2.0, (sens_valid['CI_hi'].max()
                                                  - sens_valid['CI_lo'].min()) * 0.7)
    ax_s.set_xlim(xl_s, xr_s)
    ax_s.set_ylim(-0.8, n_s - 0.2)
    ax_s.text(0.5, -0.12,
              f'<- Favors {LABEL_DYN.lower()}          Favors {LABEL_EARLY.lower()} ->',
              transform=ax_s.transAxes, ha='center', fontsize=8.5, color=GRAY_MED)
    ax_s.spines['left'].set_visible(False)
    ax_s.tick_params(left=False)
    ax_s.axhspan(y_s[0] - 0.45, y_s[0] + 0.45,
                 color='#F5F5F5', zorder=0, alpha=0.7)
    ax_s.text(xl_s + 0.02, y_s[0] + 0.5, 'Primary',
              fontsize=7.5, color=GRAY_MED, style='italic', va='bottom')
    save_fig(fig_s, r_dir, 'R_eFig_sensitivity_forest')

    print(f"[{cohort_name}] figures complete -> {r_dir}")


def main():
    if len(sys.argv) > 1 and str(sys.argv[1]).strip():
        target = os.path.abspath(str(sys.argv[1]).strip())
        render_cohort(target)
        return
    # No arg: process both default cohorts
    for sub in ('cohort_3matched', 'cohort_7group_psm02'):
        render_cohort(os.path.join(DEFAULT_BASE, sub))


if __name__ == '__main__':
    main()
