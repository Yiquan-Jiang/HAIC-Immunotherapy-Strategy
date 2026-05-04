#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PSM 前后对比森林图 — update_group_7（7组，21组两两对比）
每个对比组展示两行：上行=PSM前（灰），下行=PSM后（蓝/朱红）
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
import pandas as pd
import os
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test

BASE_DIR   = "/Users/yqj/Nutstore Files/我的坚果云/Liver_tumor_big_data/FIRST_LINE_HAIC_2025-12-30/HAIC_NO_TACE_4_TIDY/update_group_7"
DATA_DIR   = os.path.join(BASE_DIR, 'data')
RES_DIR    = os.path.join(BASE_DIR, 'results', 'psm_balance_tables_complete')
OUTPUT_DIR = os.path.join(BASE_DIR, 'figures', 'psm_pub_quality')
os.makedirs(OUTPUT_DIR, exist_ok=True)

plt.rcParams.update({
    'font.family':        'sans-serif',
    'font.sans-serif':    ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size':           8,
    'axes.labelsize':      9,
    'axes.titlesize':     10,
    'xtick.labelsize':     7.5,
    'ytick.labelsize':     7.5,
    'axes.linewidth':      0.6,
    'axes.spines.top':     False,
    'axes.spines.right':   False,
    'xtick.major.width':   0.6,
    'ytick.major.width':   0.6,
    'xtick.major.size':    3.0,
    'ytick.major.size':    3.0,
    'savefig.dpi':        300,
    'savefig.bbox':       'tight',
    'savefig.pad_inches':  0.05,
    'pdf.fonttype':        42,
    'ps.fonttype':         42,
    'axes.grid':           False,
})

COLOR_BEFORE = '#999999'
COLOR_AFTER  = '#0072B2'
COLOR_SIG    = '#D55E00'
NULL_COLOR   = '#444444'

GROUP_ORDER = [
    'HAIC_alone', 'HAIC+I_concurrent', 'HAIC_then_I',
    'HAIC+T_concurrent', 'HAIC_then_T',
    'HAIC+I+T_concurrent', 'HAIC_then_I+T',
]
GROUP_LABELS = {
    'HAIC_alone':            'HAIC alone',
    'HAIC+I_concurrent':     'HAIC + Immuno (conc.)',
    'HAIC_then_I':           'HAIC → Immuno',
    'HAIC+T_concurrent':     'HAIC + Target (conc.)',
    'HAIC_then_T':           'HAIC → Target',
    'HAIC+I+T_concurrent':   'HAIC + I + T (conc.)',
    'HAIC_then_I+T':         'HAIC → I + T',
}

# ════════════════════════════════════════════════════════════════════
# 1. 计算 PSM 前 HR
# ════════════════════════════════════════════════════════════════════
print('1. 计算 PSM 前 HR...')

df = pd.read_csv(os.path.join(DATA_DIR, 'analysis_ready.csv'))
df = df[df['os_months'] >= 0].copy()
df['group'] = df['main_group']
df['death_status'] = df['death_status'].map(
    {'Yes': 1, 'No': 0, '1': 1, '0': 0, 1: 1, 0: 0}
).fillna(0).astype(int)

COMPARISONS = []
comp_idx = 1
for i, g1 in enumerate(GROUP_ORDER):
    for j in range(i+1, len(GROUP_ORDER)):
        g2 = GROUP_ORDER[j]
        COMPARISONS.append((comp_idx, g1, g2))
        comp_idx += 1


def calc_hr(sub1, sub2):
    lr  = logrank_test(sub1['os_months'], sub2['os_months'],
                       sub1['death_status'], sub2['death_status'])
    tmp = pd.concat([sub1, sub2])[['os_months', 'death_status']].copy()
    tmp['treat'] = ([0]*len(sub1) + [1]*len(sub2))
    cph = CoxPHFitter()
    cph.fit(tmp, duration_col='os_months', event_col='death_status')
    hr    = float(np.exp(cph.params_['treat']))
    ci_lo = float(np.exp(cph.confidence_intervals_['95% lower-bound']['treat']))
    ci_hi = float(np.exp(cph.confidence_intervals_['95% upper-bound']['treat']))
    med1  = KaplanMeierFitter().fit(sub1['os_months'], sub1['death_status']).median_survival_time_
    med2  = KaplanMeierFitter().fit(sub2['os_months'], sub2['death_status']).median_survival_time_
    return hr, ci_lo, ci_hi, lr.p_value, med1, med2


before_rows = []
for cid, g1, g2 in COMPARISONS:
    s1 = df[df['group'] == g1]
    s2 = df[df['group'] == g2]
    hr, ci_lo, ci_hi, p, med1, med2 = calc_hr(s1, s2)
    before_rows.append({
        'Group1': g1, 'Group2': g2,
        'N1': len(s1), 'N2': len(s2),
        'HR': hr, 'CI_lower': ci_lo, 'CI_upper': ci_hi,
        'P_value': p, 'Median_OS_1': med1, 'Median_OS_2': med2
    })
    print(f'  {cid:02d} Before: {g1} vs {g2} | HR={hr:.2f} ({ci_lo:.2f}-{ci_hi:.2f}) P={p:.4f}')

df_before = pd.DataFrame(before_rows)

# ════════════════════════════════════════════════════════════════════
# 2. 读取 PSM 后数据
# ════════════════════════════════════════════════════════════════════
print('\n2. 读取 PSM 后 HR...')
df_after = pd.read_csv(os.path.join(RES_DIR, 'survival_analysis_final.csv'))
df_after = df_after.rename(columns={'N1_after': 'N1', 'N2_after': 'N2'})
print(f'   共 {len(df_after)} 组对比')

# 保存汇总 CSV
summary_df = pd.merge(
    df_before.rename(columns=lambda c: f'{c}_before' if c not in ('Group1','Group2') else c),
    df_after[['Group1','Group2','N1','N2','HR','CI_lower','CI_upper','P_value',
              'Median_OS_1','Median_OS_2']].rename(
        columns=lambda c: f'{c}_after' if c not in ('Group1','Group2') else c),
    on=['Group1','Group2'], how='outer'
)
summary_df.to_csv(os.path.join(OUTPUT_DIR, 'survival_summary_before_after.csv'), index=False)

# ════════════════════════════════════════════════════════════════════
# 3. 绘制双排森林图
# ════════════════════════════════════════════════════════════════════
print('\n3. 绘制森林图...')

n_comp  = len(COMPARISONS)
GAP_IN  = 0.32
GAP_OUT = 0.55

y_before = []
y_after  = []
y_sep    = []
y = 0.0
for i in range(n_comp):
    y_before.append(y + GAP_IN / 2)
    y_after.append(y - GAP_IN / 2)
    y_sep.append(y - GAP_IN / 2 - GAP_OUT / 2)
    y -= (GAP_IN + GAP_OUT)

y_before = y_before[::-1]
y_after  = y_after[::-1]
y_sep    = y_sep[::-1]

fig_h = n_comp * (GAP_IN + GAP_OUT) + 2.5
fig, ax = plt.subplots(figsize=(9.0, fig_h))

X_MIN, X_MAX = 0.12, 6.0
ax.set_xlim(X_MIN, X_MAX)
ax.set_xscale('log')

y_all = y_before + y_after
ax.set_ylim(min(y_all) - 0.8, max(y_all) + 0.9)
ax.set_yticks([])
ax.spines['left'].set_visible(False)

ax.axvline(x=1.0, color=NULL_COLOR, linestyle='--',
           linewidth=0.8, alpha=0.7, zorder=1)

y_top = max(y_all) + 0.65
ax.text(0.01, y_top, 'Comparison', ha='left', va='center',
        fontsize=8, fontweight='bold', color='#222222',
        transform=ax.get_yaxis_transform())
ax.text(1.0, y_top, 'HR (95% CI)', ha='center', va='center',
        fontsize=8, fontweight='bold', color='#222222')
ax.text(X_MAX * 1.02, y_top, 'HR (95% CI)  P-value  Median OS',
        ha='left', va='center', fontsize=7, fontweight='bold',
        color='#222222')

ax.axhline(y=max(y_all) + 0.45, color='#333333', linewidth=0.8,
           xmin=0, xmax=1, clip_on=False)

for i in range(n_comp):
    if i % 2 == 0:
        yb = y_after[i] - GAP_OUT / 2
        yt = y_before[i] + GAP_IN / 2 + 0.05
        ax.axhspan(yb, yt, color='#F7F7F7', zorder=0, linewidth=0)

for i, (cid, g1, g2) in enumerate(COMPARISONS):
    row_b = df_before.iloc[i]

    # 匹配 df_after 中对应的行
    mask_after = (df_after['Group1'] == g1) & (df_after['Group2'] == g2)
    if mask_after.sum() == 0:
        continue
    row_a = df_after[mask_after].iloc[0]

    comp_label = (f'{GROUP_LABELS[g1]}\n'
                  f'vs {GROUP_LABELS[g2]}')

    ax.text(0.01, (y_before[i] + y_after[i]) / 2,
            comp_label, ha='left', va='center',
            fontsize=6.5, color='#222222',
            transform=ax.get_yaxis_transform())

    for row, ypos, color, label in [
        (row_b, y_before[i], COLOR_BEFORE, 'Before PSM'),
        (row_a, y_after[i],  COLOR_AFTER,  'After PSM'),
    ]:
        hr    = float(row['HR'])
        ci_lo = float(row['CI_lower'])
        ci_hi = float(row['CI_upper'])
        p_val = float(row['P_value'])
        med1  = float(row['Median_OS_1'])
        med2  = float(row['Median_OS_2'])

        dot_color = COLOR_SIG if (label == 'After PSM' and p_val < 0.05) else color

        ci_lo_plot = max(ci_lo, X_MIN * 1.01)
        ci_hi_plot = min(ci_hi, X_MAX * 0.99)
        ax.plot([ci_lo_plot, ci_hi_plot], [ypos, ypos],
                color=dot_color, linewidth=1.4,
                solid_capstyle='round', zorder=3)

        for xc in [ci_lo_plot, ci_hi_plot]:
            ax.plot([xc, xc], [ypos - 0.07, ypos + 0.07],
                    color=dot_color, linewidth=1.0, zorder=3)

        ax.scatter([hr], [ypos], marker='D', s=30,
                   color=dot_color, zorder=4, linewidths=0)

        p_str = '<0.001' if p_val < 0.001 else f'{p_val:.3f}'
        med_str = f'{med1:.1f} vs {med2:.1f}'
        right_txt = f'{hr:.2f} ({ci_lo:.2f}\u2013{ci_hi:.2f})   {p_str}   {med_str} mo'
        ax.text(X_MAX * 1.02, ypos, right_txt,
                ha='left', va='center', fontsize=6,
                color='#333333' if p_val >= 0.05 else '#B03000')

    if i < n_comp - 1:
        ax.axhline(y=y_sep[i], color='#DDDDDD', linewidth=0.5,
                   xmin=0, xmax=1, clip_on=False)

ax.axhline(y=min(y_all) - 0.45, color='#333333', linewidth=0.8,
           xmin=0, xmax=1, clip_on=False)

ax.set_xlabel('Hazard Ratio (95% CI)  [log scale]', fontsize=9, labelpad=6)
xticks = [0.15, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0]
ax.set_xticks(xticks)
ax.set_xticklabels([str(x) for x in xticks], fontsize=7)

ax.text(0.40, -0.02, '\u2190 Favors Group 2',
        ha='center', va='top', fontsize=7, color='#666666',
        transform=ax.transAxes)
ax.text(0.72, -0.02, 'Favors Group 1 \u2192',
        ha='center', va='top', fontsize=7, color='#666666',
        transform=ax.transAxes)

legend_elements = [
    mlines.Line2D([0], [0], color=COLOR_BEFORE, linewidth=1.5,
                  marker='D', markersize=4, markerfacecolor=COLOR_BEFORE,
                  label='Before PSM (unadjusted)'),
    mlines.Line2D([0], [0], color=COLOR_AFTER, linewidth=1.5,
                  marker='D', markersize=4, markerfacecolor=COLOR_AFTER,
                  label='After PSM (matched)'),
    mlines.Line2D([0], [0], color=COLOR_SIG, linewidth=1.5,
                  marker='D', markersize=4, markerfacecolor=COLOR_SIG,
                  label='After PSM, P < 0.05'),
]
ax.legend(handles=legend_elements, loc='lower left',
          bbox_to_anchor=(0.0, -0.06), ncol=3,
          fontsize=7, frameon=False, handlelength=1.8)

ax.set_title('Overall Survival \u2014 Hazard Ratio Summary (21 Pairwise Comparisons)\nBefore vs After PSM',
             fontsize=10, fontweight='bold', pad=10, loc='left')

stats_box = dict(boxstyle='round,pad=0.4', facecolor='white',
                 edgecolor='#CCCCCC', alpha=0.95, linewidth=0.5)
ax.text(0.99, 0.01,
        'PSM: 1:1 nearest-neighbor matching\nReference: HR = 1.0 (dashed line)',
        transform=ax.transAxes, fontsize=6, ha='right', va='bottom',
        bbox=stats_box, color='#555555')

plt.tight_layout()

out_base = os.path.join(OUTPUT_DIR, 'HR_forest_plot_summary')
fig.savefig(f'{out_base}.pdf', bbox_inches='tight', pad_inches=0.05)
fig.savefig(f'{out_base}.png', dpi=300, bbox_inches='tight', pad_inches=0.05)
plt.close(fig)
print(f'\n✅ 已保存: {out_base}.pdf / .png')
