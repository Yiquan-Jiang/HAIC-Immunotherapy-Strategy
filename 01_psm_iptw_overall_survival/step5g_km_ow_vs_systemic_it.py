#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""step5g — ATO overlap-weighted KM: each HAIC strategy vs Systemic I+T.

Same publication style as the pairwise KM in step4_km_curves.py (two panels with risk
tables + stats boxes), but the contrast is the no-HAIC arm, so the two panels are
"Before weighting" (unweighted) | "After OW (ATO)" (balanced overlap weights). This is the
KM analogue of the step5f forest and is the balanced replacement for the step4 PSM-matched /
step4b joint-IPTW KM, which do NOT balance the Systemic arm (post-match |SMD| 0.10-0.39;
joint IPTW |SMD| 0.335). One figure per contrast.

Inputs (from step5e): results/ow_vs_systemic_it_8group/ow_weights.csv, ow_forest_data.csv
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from lifelines import KaplanMeierFitter, CoxPHFitter

BASE_DIR = "/Users/yqj/Nutstore Files/我的坚果云/Liver_tumor_big_data/FIRST_LINE_HAIC_2025-12-30/HAIC_NO_TACE_4_TIDY/update_group_7"
RES_DIR = os.path.join(BASE_DIR, "results", "ow_vs_systemic_it_8group")
FIG_DIR = os.path.join(BASE_DIR, "figures", "ow_vs_systemic_it_8group")
os.makedirs(FIG_DIR, exist_ok=True)

# ── style copied from step4_km_curves.py for visual consistency ──────────────
GROUP_COLORS = {
    'HAIC_alone': '#0072B2', 'HAIC+I_concurrent': '#E69F00', 'HAIC_then_I': '#009E73',
    'HAIC+T_concurrent': '#F0E442', 'HAIC_then_T': '#CC79A7',
    'HAIC+I+T_concurrent': '#D55E00', 'HAIC_then_I+T': '#56B4E9', 'Systemic_I+T': '#000000',
}
GROUP_LABELS = {
    'HAIC_alone': 'HAIC alone', 'HAIC+I_concurrent': 'HAIC + Immuno (concurrent)',
    'HAIC_then_I': 'HAIC → Immuno', 'HAIC+T_concurrent': 'HAIC + Antiangiogenic (concurrent)',
    'HAIC_then_T': 'HAIC → Antiangiogenic', 'HAIC+I+T_concurrent': 'HAIC + Immuno + Antiangiogenic',
    'HAIC_then_I+T': 'HAIC → Immuno + Antiangiogenic', 'Systemic_I+T': 'Systemic I+T',
}
plt.rcParams.update({
    'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 10, 'axes.labelsize': 11, 'axes.titlesize': 12,
    'xtick.labelsize': 9, 'ytick.labelsize': 9, 'legend.fontsize': 8,
    'axes.linewidth': 0.8, 'axes.spines.top': False, 'axes.spines.right': False,
    'xtick.major.width': 0.8, 'ytick.major.width': 0.8, 'xtick.major.size': 3.5,
    'ytick.major.size': 3.5, 'lines.linewidth': 2.0, 'legend.frameon': False,
    'legend.borderpad': 0.3, 'legend.handlelength': 2.0, 'savefig.dpi': 300,
    'savefig.bbox': 'tight', 'savefig.pad_inches': 0.05, 'pdf.fonttype': 42,
    'ps.fonttype': 42, 'axes.grid': False,
})
STATS_BOX = dict(boxstyle='round,pad=0.4', facecolor='white', edgecolor='#CCCCCC',
                 alpha=0.95, linewidth=0.5)
RISK_TIMES = [0, 12, 24, 36, 48, 60]


def draw_km(ax_km, ax_risk, groups_data, title, xlim=60):
    """KM curves + No.-at-risk table. 'risk_kmf' (unweighted) drives the risk counts when
    the curve kmf is weighted; defaults to the curve kmf."""
    n_groups = len(groups_data)
    for gd in groups_data:
        kmf = gd['kmf']
        rk = gd.get('risk_kmf', kmf)
        n = int(rk.event_table['at_risk'].iloc[0])
        t = kmf.survival_function_.index.values
        s = kmf.survival_function_.iloc[:, 0].values
        ci_lo = kmf.confidence_interval_.iloc[:, 0].values
        ci_hi = kmf.confidence_interval_.iloc[:, 1].values
        ax_km.step(t, s, where='post', color=gd['color'], linewidth=1.8,
                   linestyle=gd['linestyle'], label=f"{gd['label']} (n={n:,})")
        ax_km.fill_between(t, ci_lo, ci_hi, step='post', alpha=0.08, color=gd['color'])
    ax_km.set_xlim(0, xlim)
    ax_km.set_ylim(-0.02, 1.05)
    ax_km.set_ylabel('Overall Survival Probability')
    ax_km.set_title(title, fontweight='bold', pad=8)
    ax_km.set_xticks(RISK_TIMES)
    ax_km.yaxis.set_major_locator(mticker.MultipleLocator(0.2))
    ax_km.legend(loc='upper right', handlelength=2.0, fontsize=6.5)

    ax_risk.set_xlim(0, xlim + 0.5)
    ax_risk.set_ylim(-0.8, n_groups + 0.2)
    ax_risk.axis('off')
    ax_risk.text(xlim / 2, n_groups - 0.1, 'No. at risk', ha='center', va='top',
                 fontsize=9, color='#666666', style='italic')
    for i, gd in enumerate(groups_data):
        y = n_groups - i - 0.9
        rk = gd.get('risk_kmf', gd['kmf'])
        for t_pt in RISK_TIMES:
            idx = rk.event_table.index[rk.event_table.index <= t_pt]
            n = int(rk.event_table.loc[idx[-1], 'at_risk']) if len(idx) else 0
            ax_risk.text(t_pt, y, str(n), ha='center', va='center', fontsize=8,
                         color=gd['color'], fontweight='bold')
    ax_risk.text(xlim / 2, -0.65, 'Time (months)', ha='center', va='center',
                 fontsize=11, color='#333333')


def compute_cox_hr(t1, e1, t2, e2, w1=None, w2=None):
    """HR for group1 vs group2 (group2 reference). Optional weights -> robust SE."""
    df = pd.DataFrame({
        'T': np.concatenate([t1, t2]), 'E': np.concatenate([e1, e2]),
        'group': [1] * len(t1) + [0] * len(t2)})
    kw = {}
    if w1 is not None:
        df['w'] = np.concatenate([w1, w2])
        kw = dict(weights_col='w', robust=True)
    cph = CoxPHFitter()
    cph.fit(df, duration_col='T', event_col='E', **kw)
    hr = float(np.exp(cph.params_['group']))
    ci = np.exp(cph.confidence_intervals_.loc['group'].values)
    return hr, float(ci[0]), float(ci[1])


def add_stats_box(ax, p_val, hr, ci_lo, ci_hi, med1, med2, label1, label2, p_label='Log-rank'):
    p_fmt = 'P < 0.001' if p_val < 0.001 else f'P = {p_val:.3f}'
    lines = [f'{p_label} {p_fmt}',
             f'HR = {hr:.2f} (95% CI {ci_lo:.2f}–{ci_hi:.2f})',
             f'  (ref: {label2})',
             f'Median OS: {label1}: {med1:.1f} mo',
             f'  vs {label2}: {med2:.1f} mo']
    ax.text(0.97, 0.05, '\n'.join(lines), transform=ax.transAxes, fontsize=7.5,
            va='bottom', ha='right', bbox=STATS_BOX)


def save_fig(fig, name):
    base = os.path.join(FIG_DIR, name)
    fig.savefig(f'{base}.pdf', bbox_inches='tight', pad_inches=0.05)
    fig.savefig(f'{base}.png', dpi=300, bbox_inches='tight', pad_inches=0.05)
    plt.close(fig)
    print(f'  Saved: {name}.pdf/.png')


# ── data ─────────────────────────────────────────────────────────────────────
W = pd.read_csv(os.path.join(RES_DIR, 'ow_weights.csv'))
fd = pd.read_csv(os.path.join(RES_DIR, 'ow_forest_data.csv')).set_index('group')
REF = 'Systemic_I+T'
HAIC_ORDER = ['HAIC_alone', 'HAIC+I_concurrent', 'HAIC_then_I', 'HAIC+T_concurrent',
              'HAIC_then_T', 'HAIC+I+T_concurrent', 'HAIC_then_I+T']

for cid, g in enumerate(HAIC_ORDER, 1):
    sub = W[W['group'] == g]
    a = sub[sub['main_group'] == g]            # HAIC arm (treat=1)
    b = sub[sub['main_group'] == REF]          # Systemic arm (treat=0)
    n1, n2 = len(a), len(b)

    # unweighted (Before)
    k1u, k2u = KaplanMeierFitter(), KaplanMeierFitter()
    k1u.fit(a['os_months'], a['death'], label=GROUP_LABELS[g])
    k2u.fit(b['os_months'], b['death'], label=GROUP_LABELS[REF])
    from lifelines.statistics import logrank_test
    p_un = logrank_test(a['os_months'], b['os_months'], a['death'], b['death']).p_value
    hr_un, lo_un, hi_un = compute_cox_hr(a['os_months'].values, a['death'].values,
                                         b['os_months'].values, b['death'].values)

    # OW-weighted (After)
    k1w, k2w = KaplanMeierFitter(), KaplanMeierFitter()
    k1w.fit(a['os_months'], a['death'], weights=a['ow'].values, label=GROUP_LABELS[g])
    k2w.fit(b['os_months'], b['death'], weights=b['ow'].values, label=GROUP_LABELS[REF])
    r = fd.loc[g]

    fig = plt.figure(figsize=(10.5, 5.5))
    outer = fig.add_gridspec(1, 2, wspace=0.30)
    gs_l = outer[0].subgridspec(2, 1, height_ratios=[4.2, 0.8], hspace=0.04)
    ax_l, ax_lr = fig.add_subplot(gs_l[0]), fig.add_subplot(gs_l[1])
    gs_r = outer[1].subgridspec(2, 1, height_ratios=[4.2, 0.8], hspace=0.04)
    ax_r, ax_rr = fig.add_subplot(gs_r[0]), fig.add_subplot(gs_r[1])

    gd_un = [
        {'label': GROUP_LABELS[g], 'color': GROUP_COLORS[g], 'linestyle': '-', 'kmf': k1u},
        {'label': GROUP_LABELS[REF], 'color': GROUP_COLORS[REF], 'linestyle': '--', 'kmf': k2u},
    ]
    gd_ow = [
        {'label': GROUP_LABELS[g], 'color': GROUP_COLORS[g], 'linestyle': '-',
         'kmf': k1w, 'risk_kmf': k1u},
        {'label': GROUP_LABELS[REF], 'color': GROUP_COLORS[REF], 'linestyle': '--',
         'kmf': k2w, 'risk_kmf': k2u},
    ]
    draw_km(ax_l, ax_lr, gd_un, title=f'Before weighting  (n={n1}+{n2})')
    add_stats_box(ax_l, p_un, hr_un, lo_un, hi_un,
                  k1u.median_survival_time_, k2u.median_survival_time_,
                  GROUP_LABELS[g], GROUP_LABELS[REF], p_label='Log-rank')

    draw_km(ax_r, ax_rr, gd_ow,
            title=f'After OW · ATO  (ESS {int(r["ess_haic"])}+{int(r["ess_sys"])}, '
                  f'max|SMD| {r["max_smd_adj"]:.2f})')
    add_stats_box(ax_r, r['p'], r['HR'], r['CI_lower'], r['CI_upper'],
                  k1w.median_survival_time_, k2w.median_survival_time_,
                  GROUP_LABELS[g], GROUP_LABELS[REF], p_label='Weighted Cox')

    for ax, lbl in [(ax_l, 'A'), (ax_r, 'B')]:
        ax.text(-0.12, 1.08, lbl, transform=ax.transAxes, fontsize=10,
                fontweight='bold', va='top')
    fig.suptitle(f'{GROUP_LABELS[g]}  vs  Systemic I+T   (overlap-weighted, ATO)',
                 fontsize=11, fontweight='bold', y=1.01)
    save_fig(fig, f'km_ow_{cid:02d}_{g}_vs_Systemic_I+T'.replace('+', '_').replace('→', '_'))

print('Done.')
