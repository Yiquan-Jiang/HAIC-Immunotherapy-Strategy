#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
发表级 KM 生存曲线 — update_group_7（7组，21组两两对比）
数据源: update_group_7/data/analysis_ready.csv
PSM 匹配 ID: results/psm_balance_tables_complete/matched_ids_*.csv
PSM 生存分析: results/psm_balance_tables_complete/survival_analysis_final.csv

输出:
  figures/km/01_km_all_before_psm.pdf/png          7组整体
  figures/km/02_km_comp01_*.pdf/png ... 22_km_comp21_*.pdf/png   21组两两
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
import os
import glob
from itertools import combinations
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test, multivariate_logrank_test

BASE_DIR = "/Users/yqj/Nutstore Files/我的坚果云/Liver_tumor_big_data/FIRST_LINE_HAIC_2025-12-30/HAIC_NO_TACE_4_TIDY/update_group_7"
EIGHT_GROUP = os.environ.get("EIGHT_GROUP", "0") == "1"
SFX = "_8group" if EIGHT_GROUP else ""
DATA_CSV = "analysis_ready_8group.csv" if EIGHT_GROUP else "analysis_ready.csv"
DATA_DIR = os.path.join(BASE_DIR, 'data')
FIG_DIR  = os.path.join(BASE_DIR, 'figures', 'km' + SFX)
RES_DIR  = os.path.join(BASE_DIR, 'results', 'psm_balance_tables_complete' + SFX)
os.makedirs(FIG_DIR, exist_ok=True)

GROUP_COLORS = {
    'HAIC_alone':            '#0072B2',
    'HAIC+I_concurrent':     '#E69F00',
    'HAIC_then_I':           '#009E73',
    'HAIC+T_concurrent':     '#F0E442',
    'HAIC_then_T':           '#CC79A7',
    'HAIC+I+T_concurrent':   '#D55E00',
    'HAIC_then_I+T':         '#56B4E9',
    'Systemic_I+T':          '#117733',
}
GROUP_LABELS = {
    'HAIC_alone':            'HAIC alone',
    'HAIC+I_concurrent':     'HAIC + Immuno (concurrent)',
    'HAIC_then_I':           'HAIC → Immuno',
    'HAIC+T_concurrent':     'HAIC + Antiangiogenic (concurrent)',
    'HAIC_then_T':           'HAIC → Antiangiogenic',
    'HAIC+I+T_concurrent':   'HAIC + Immuno + Antiangiogenic',
    'HAIC_then_I+T':         'HAIC → Immuno + Antiangiogenic',
    'Systemic_I+T':          'Systemic Immuno + Antiangiogenic',
}
GROUP_ORDER = [
    'HAIC_alone', 'HAIC+I_concurrent', 'HAIC_then_I', 'HAIC+T_concurrent',
    'HAIC_then_T', 'HAIC+I+T_concurrent', 'HAIC_then_I+T',
] + (["Systemic_I+T"] if EIGHT_GROUP else [])
LINESTYLES  = ['-', '--', '-.', ':', (0, (3, 1, 1, 1)), (0, (5, 2)), (0, (1, 1)), (0, (3, 1, 1, 1, 1, 1))]

plt.rcParams.update({
    'font.family':        'sans-serif',
    'font.sans-serif':    ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size':          10,
    'axes.labelsize':     11,
    'axes.titlesize':     12,
    'xtick.labelsize':     9,
    'ytick.labelsize':     9,
    'legend.fontsize':     8,
    'axes.linewidth':      0.8,
    'axes.spines.top':     False,
    'axes.spines.right':   False,
    'xtick.major.width':   0.8,
    'ytick.major.width':   0.8,
    'xtick.major.size':    3.5,
    'ytick.major.size':    3.5,
    'lines.linewidth':     2.0,
    'legend.frameon':      False,
    'legend.borderpad':    0.3,
    'legend.handlelength': 2.0,
    'savefig.dpi':        300,
    'savefig.bbox':       'tight',
    'savefig.pad_inches':  0.05,
    'pdf.fonttype':        42,
    'ps.fonttype':         42,
    'axes.grid':           False,
})

STATS_BOX  = dict(boxstyle='round,pad=0.4', facecolor='white',
                  edgecolor='#CCCCCC', alpha=0.95, linewidth=0.5)
RISK_TIMES = [0, 12, 24, 36, 48, 60]


def draw_km(ax_km, ax_risk, groups_data, title, xlim=60):
    n_groups = len(groups_data)
    for gd in groups_data:
        kmf   = gd['kmf']
        color = gd['color']
        ls    = gd['linestyle']
        label = gd['label']
        n     = int(kmf.event_table['at_risk'].iloc[0])
        t     = kmf.survival_function_.index.values
        s     = kmf.survival_function_.iloc[:, 0].values
        ci_lo = kmf.confidence_interval_.iloc[:, 0].values
        ci_hi = kmf.confidence_interval_.iloc[:, 1].values
        ax_km.step(t, s, where='post', color=color, linewidth=1.8,
                   linestyle=ls, label=f'{label} (n={n:,})')
        ax_km.fill_between(t, ci_lo, ci_hi, step='post',
                           alpha=0.08, color=color)
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
    ax_risk.text(xlim / 2, n_groups - 0.1, 'No. at risk',
                 ha='center', va='top', fontsize=9,
                 color='#666666', style='italic')
    for i, gd in enumerate(groups_data):
        y   = n_groups - i - 0.9
        kmf = gd['kmf']
        for t_pt in RISK_TIMES:
            idx = kmf.event_table.index[kmf.event_table.index <= t_pt]
            n   = int(kmf.event_table.loc[idx[-1], 'at_risk']) if len(idx) else 0
            ax_risk.text(t_pt, y, str(n), ha='center', va='center',
                         fontsize=8, color=gd['color'], fontweight='bold')
    ax_risk.text(xlim / 2, -0.65, 'Time (months)',
                 ha='center', va='center', fontsize=11, color='#333333')


def compute_cox_hr(t1, e1, t2, e2):
    """Compute HR for group1 vs group2 (group2 as reference).

    "A_vs_B" → HR represents A relative to B (B is reference).
    """
    t = pd.concat([t1.reset_index(drop=True), t2.reset_index(drop=True)],
                  ignore_index=True)
    e = pd.concat([e1.reset_index(drop=True), e2.reset_index(drop=True)],
                  ignore_index=True)
    g = pd.Series([1] * len(t1) + [0] * len(t2))
    cph_df = pd.DataFrame({'T': t, 'E': e, 'group': g})
    cph = CoxPHFitter()
    cph.fit(cph_df, duration_col='T', event_col='E')
    hr = float(np.exp(cph.params_['group']))
    ci = np.exp(cph.confidence_intervals_.loc['group'].values)
    return hr, float(ci[0]), float(ci[1])


def add_stats_box(ax, p_val, hr=None, ci_lo=None, ci_hi=None,
                  med1=None, med2=None, label1=None, label2=None):
    """统计信息标注框，HR 标注参考组，Median 标注各组归属。"""
    p_fmt = 'P < 0.001' if p_val < 0.001 else f'P = {p_val:.3f}'
    lines = [f'Log-rank {p_fmt}']
    if hr is not None:
        lines.append(f'HR = {hr:.2f} (95% CI {ci_lo:.2f}\u2013{ci_hi:.2f})')
        if label2:
            lines.append(f'  (ref: {label2})')
    if med1 is not None and med2 is not None:
        if label1 and label2:
            lines.append(f'Median OS: {label1}: {med1:.1f} mo')
            lines.append(f'  vs {label2}: {med2:.1f} mo')
        else:
            lines.append(f'Median OS: {med1:.1f} vs {med2:.1f} mo')
    ax.text(0.97, 0.05, '\n'.join(lines), transform=ax.transAxes,
            fontsize=7.5, va='bottom', ha='right', bbox=STATS_BOX)


def save_fig(fig, name):
    base = os.path.join(FIG_DIR, name)
    fig.savefig(f'{base}.pdf', bbox_inches='tight', pad_inches=0.05)
    fig.savefig(f'{base}.png', dpi=300, bbox_inches='tight', pad_inches=0.05)
    plt.close(fig)
    print(f'  Saved: {name}.pdf/.png')


# ════════════════════════════════════════════════════════════════════
# 1. 读取数据
# ════════════════════════════════════════════════════════════════════
print('1. 读取数据...')

df = pd.read_csv(os.path.join(DATA_DIR, DATA_CSV))
df = df[df['os_months'] >= 0].copy()
df['group'] = pd.Categorical(df['main_group'], categories=GROUP_ORDER, ordered=True)
df['death_status'] = df['death_status'].map(
    {'Yes': 1, 'No': 0, '1': 1, '0': 0, 1: 1, 0: 0}
).fillna(0).astype(int)

print(f'   总样本量: {len(df)} 患者')
print(df['group'].value_counts().sort_index().to_string())


# ════════════════════════════════════════════════════════════════════
# 2. PSM 前整体 KM 曲线（7组）
# ════════════════════════════════════════════════════════════════════
print('\n2. PSM 前整体 KM 曲线...')

fig = plt.figure(figsize=(6.0, 6.5))
gs  = fig.add_gridspec(2, 1, height_ratios=[4.0, 1.5], hspace=0.04)
ax_km   = fig.add_subplot(gs[0])
ax_risk = fig.add_subplot(gs[1])

groups_data_all = []
for i, grp in enumerate(GROUP_ORDER):
    sub = df[df['group'] == grp]
    if len(sub) == 0:
        continue
    kmf = KaplanMeierFitter()
    kmf.fit(sub['os_months'], sub['death_status'], label=GROUP_LABELS[grp])
    groups_data_all.append({
        'label':     GROUP_LABELS[grp],
        'color':     GROUP_COLORS[grp],
        'linestyle': LINESTYLES[i],
        'kmf':       kmf,
    })

draw_km(ax_km, ax_risk, groups_data_all,
        title=f'Overall Survival \u2014 All Groups (Before PSM, N={len(df):,})')

results_all = multivariate_logrank_test(df['os_months'], df['group'], df['death_status'])
p_all = results_all.p_value
p_fmt = 'P < 0.001' if p_all < 0.001 else f'P = {p_all:.3f}'
ax_km.text(0.97, 0.05, f'Log-rank {p_fmt}', transform=ax_km.transAxes,
           fontsize=7, va='bottom', ha='right', bbox=STATS_BOX)

save_fig(fig, '01_km_all_before_psm')


# ════════════════════════════════════════════════════════════════════
# 3. 读取 R（MatchIt）匹配患者 ID 列表 + 生存分析结果
# ════════════════════════════════════════════════════════════════════
print('\n3. 读取 R PSM 匹配结果...')

r_results_path = os.path.join(RES_DIR, 'survival_analysis_final.csv')
r_results = None
if os.path.exists(r_results_path):
    r_results = pd.read_csv(r_results_path)
    print(f'   已读取 R 生存分析结果: {len(r_results)} 组')
else:
    print(f'   ⚠ 未找到 R 生存分析结果，请先运行 step3_psm_analysis.R')

matched_id_files = sorted(glob.glob(os.path.join(RES_DIR, 'matched_ids_*.csv')))
r_matched = {}
for f in matched_id_files:
    try:
        mdf = pd.read_csv(f)
        key = mdf['comparison_key'].iloc[0]
        r_matched[key] = mdf
        print(f'   已读取: {os.path.basename(f)}  ({len(mdf)} 行)')
    except Exception as e:
        print(f'   ⚠ 读取失败 {f}: {e}')

if not r_matched:
    print('   ⚠ 未找到任何 matched_ids_*.csv，请先运行 step3_psm_analysis.R')


# ════════════════════════════════════════════════════════════════════
# 4. 21组两两对比 KM 曲线（PSM前后并排）
# ════════════════════════════════════════════════════════════════════
print('\n4. 21组两两对比 KM 曲线...')

COMPARISONS = []
comp_idx = 1
for i, g1 in enumerate(GROUP_ORDER):
    for j in range(i+1, len(GROUP_ORDER)):
        g2 = GROUP_ORDER[j]
        key = f'{g1}_vs_{g2}'
        COMPARISONS.append((comp_idx, g1, g2, key))
        comp_idx += 1

for cid, g1, g2, comp_key in COMPARISONS:
    print(f'  Comp {cid:02d}: {g1} vs {g2}')

    raw  = df[df['group'].isin([g1, g2])].copy()
    raw1 = raw[raw['group'] == g1]
    raw2 = raw[raw['group'] == g2]
    n1_raw, n2_raw = len(raw1), len(raw2)

    kmf1_raw = KaplanMeierFitter()
    kmf2_raw = KaplanMeierFitter()
    kmf1_raw.fit(raw1['os_months'], raw1['death_status'], label=GROUP_LABELS[g1])
    kmf2_raw.fit(raw2['os_months'], raw2['death_status'], label=GROUP_LABELS[g2])

    lr_raw = logrank_test(raw1['os_months'], raw2['os_months'],
                          raw1['death_status'], raw2['death_status'])
    p_raw  = lr_raw.p_value

    if comp_key not in r_matched:
        print(f'    ⚠ 未找到 R 匹配结果 ({comp_key})，跳过')
        continue

    matched_ids = r_matched[comp_key]['patient_id'].tolist()
    matched_groups = r_matched[comp_key].set_index('patient_id')['group_label']

    matched = df[df['patient_id'].isin(matched_ids)].copy()
    matched['group_r'] = matched['patient_id'].map(matched_groups)
    matched['group_r'] = matched['group_r'].fillna(matched['group'].astype(str))

    sub1 = matched[matched['group_r'] == g1]
    sub2 = matched[matched['group_r'] == g2]
    n1, n2 = len(sub1), len(sub2)

    if n1 == 0 or n2 == 0:
        print(f'    ⚠ 匹配后某组为空，跳过')
        continue

    kmf1 = KaplanMeierFitter()
    kmf2 = KaplanMeierFitter()
    kmf1.fit(sub1['os_months'], sub1['death_status'], label=GROUP_LABELS[g1])
    kmf2.fit(sub2['os_months'], sub2['death_status'], label=GROUP_LABELS[g2])

    lr_res = logrank_test(sub1['os_months'], sub2['os_months'],
                          sub1['death_status'], sub2['death_status'])
    p_val  = lr_res.p_value

    hr = ci_lo = ci_hi = med1 = med2 = None
    if r_results is not None:
        row = r_results[r_results['Comparison'].str.startswith(f'{cid:02d}_')]
        if len(row):
            hr    = float(row['HR'].iloc[0])
            ci_lo = float(row['CI_lower'].iloc[0])
            ci_hi = float(row['CI_upper'].iloc[0])
            med1  = float(row['Median_OS_1'].iloc[0])
            med2  = float(row['Median_OS_2'].iloc[0])

    if med1 is None:
        med1 = kmf1.median_survival_time_
        med2 = kmf2.median_survival_time_

    print(f'    Before PSM: n={n1_raw}+{n2_raw} | P={p_raw:.4f}')
    print(f'    After  PSM: n={n1}+{n2} | P={p_val:.4f}' +
          (f' | HR={hr:.2f}' if hr else ''))

    fig = plt.figure(figsize=(10.5, 5.5))
    outer = fig.add_gridspec(1, 2, wspace=0.30)

    gs_l = outer[0].subgridspec(2, 1, height_ratios=[4.2, 0.8], hspace=0.04)
    ax_km_l   = fig.add_subplot(gs_l[0])
    ax_risk_l = fig.add_subplot(gs_l[1])

    gs_r = outer[1].subgridspec(2, 1, height_ratios=[4.2, 0.8], hspace=0.04)
    ax_km_r   = fig.add_subplot(gs_r[0])
    ax_risk_r = fig.add_subplot(gs_r[1])

    gd_raw = [
        {'label': GROUP_LABELS[g1], 'color': GROUP_COLORS[g1],
         'linestyle': '-',  'kmf': kmf1_raw},
        {'label': GROUP_LABELS[g2], 'color': GROUP_COLORS[g2],
         'linestyle': '--', 'kmf': kmf2_raw},
    ]
    gd_psm = [
        {'label': GROUP_LABELS[g1], 'color': GROUP_COLORS[g1],
         'linestyle': '-',  'kmf': kmf1},
        {'label': GROUP_LABELS[g2], 'color': GROUP_COLORS[g2],
         'linestyle': '--', 'kmf': kmf2},
    ]

    # Before PSM: 计算 HR 和 Median OS
    try:
        hr_raw, ci_lo_raw, ci_hi_raw = compute_cox_hr(
            raw1['os_months'], raw1['death_status'],
            raw2['os_months'], raw2['death_status'])
    except Exception:
        hr_raw = ci_lo_raw = ci_hi_raw = None
    med1_raw = kmf1_raw.median_survival_time_
    med2_raw = kmf2_raw.median_survival_time_

    draw_km(ax_km_l, ax_risk_l, gd_raw,
            title=f'Before PSM  (n={n1_raw}+{n2_raw})')
    add_stats_box(ax_km_l, p_raw, hr_raw, ci_lo_raw, ci_hi_raw,
                  med1_raw, med2_raw,
                  label1=GROUP_LABELS[g1], label2=GROUP_LABELS[g2])

    draw_km(ax_km_r, ax_risk_r, gd_psm,
            title=f'After PSM  (n={n1}+{n2})')
    add_stats_box(ax_km_r, p_val, hr, ci_lo, ci_hi, med1, med2,
                  label1=GROUP_LABELS[g1], label2=GROUP_LABELS[g2])

    for ax, lbl in [(ax_km_l, 'A'), (ax_km_r, 'B')]:
        ax.text(-0.12, 1.08, lbl, transform=ax.transAxes,
                fontsize=10, fontweight='bold', va='top')

    fig.suptitle(
        f'Comparison {cid:02d}: {GROUP_LABELS[g1]}  vs  {GROUP_LABELS[g2]}',
        fontsize=11, fontweight='bold', y=1.01
    )

    fname = (f'{cid+1:02d}_km_comp{cid:02d}_'
             f'{g1.replace("+", "_").replace(" ", "_")}_vs_'
             f'{g2.replace("+", "_").replace(" ", "_")}')
    save_fig(fig, fname)

print(f'\n全部完成！输出目录: {FIG_DIR}')
