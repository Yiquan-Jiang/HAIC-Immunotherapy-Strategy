#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CBPS-IPTW 加权后 7 组叠加 KM 生存曲线
KM 估计及 CI 来自 R survey::svykm（正确的加权方差估计）
全局 P 值来自加权 Cox 全模型 Wald 检验
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
import os

BASE_DIR = "/Users/yqj/Nutstore Files/我的坚果云/Liver_tumor_big_data/FIRST_LINE_HAIC_2025-12-30/HAIC_NO_TACE_4_TIDY/update_group_7"
EIGHT_GROUP = os.environ.get("EIGHT_GROUP", "0") == "1"
SFX = "_8group" if EIGHT_GROUP else ""
DATA_CSV = "analysis_ready_8group.csv" if EIGHT_GROUP else "analysis_ready.csv"
FIG_DIR  = os.path.join(BASE_DIR, 'figures', 'km' + SFX)
RES_DIR  = os.path.join(BASE_DIR, 'results', 'psm_vs_template' + SFX)
os.makedirs(FIG_DIR, exist_ok=True)

GROUP_ORDER = [
    'HAIC_alone', 'HAIC+I_concurrent', 'HAIC_then_I',
    'HAIC+T_concurrent', 'HAIC_then_T',
    'HAIC+I+T_concurrent', 'HAIC_then_I+T',
] + (["Systemic_I+T"] if EIGHT_GROUP else [])
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
# 与 step4_km_curves.py（01_km_all_before_psm）图例文案一致
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
LINESTYLES = ['-', '--', '-.', ':', (0, (3, 1, 1, 1)), (0, (5, 2)), (0, (1, 1)), (0, (3, 1, 1, 1, 1, 1))]

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
XLIM = 60


def save_fig(fig, name):
    base = os.path.join(FIG_DIR, name)
    fig.savefig(f'{base}.pdf', bbox_inches='tight', pad_inches=0.05)
    fig.savefig(f'{base}.png', dpi=300, bbox_inches='tight', pad_inches=0.05)
    plt.close(fig)
    print(f'  Saved: {name}.pdf/.png')


# ════════════════════════════════════════════════════════════════════
# 1. 读取 R 导出的加权 KM 数据 + 权重数据
# ════════════════════════════════════════════════════════════════════
print('1. 读取数据...')

km_file = os.path.join(RES_DIR, 'iptw_km_data.csv')
wt_file = os.path.join(RES_DIR, 'iptw_weights.csv')
gt_file = os.path.join(RES_DIR, 'iptw_global_test.csv')

for f, desc in [(km_file, 'KM 数据'), (wt_file, '权重'), (gt_file, '全局检验')]:
    if not os.path.exists(f):
        raise FileNotFoundError(f'未找到 {desc} 文件: {f}\n请先运行 step3b_psm_vs_template.R')

km_df = pd.read_csv(km_file)
wt_df = pd.read_csv(wt_file)
gt_df = pd.read_csv(gt_file)

wt_df['death_status'] = wt_df['death_status'].map(
    {'Yes': 1, 'No': 0, '1': 1, '0': 0, 1: 1, 0: 0}
).fillna(0).astype(int)

print(f'   KM 数据: {len(km_df)} 行')
print(f'   总样本量: {len(wt_df)} 患者')

# ════════════════════════════════════════════════════════════════════
# 2. 绘制 7 组叠加 KM 曲线
# ════════════════════════════════════════════════════════════════════
print('\n2. 绘制 7 组 IPTW 加权 KM 曲线 (survey::svykm CI)...')
# 画布与栅格与 01_km_all_before_psm（step4_km_curves.py）一致
fig = plt.figure(figsize=(6.0, 6.5))
gs  = fig.add_gridspec(2, 1, height_ratios=[4.0, 1.5], hspace=0.04)
ax_km   = fig.add_subplot(gs[0])
ax_risk = fig.add_subplot(gs[1])

n_groups = 0
group_info = []

for i, grp in enumerate(GROUP_ORDER):
    sub_km = km_df[km_df['group'] == grp].sort_values('time')
    sub_wt = wt_df[wt_df['main_group'] == grp]

    if len(sub_km) == 0:
        print(f'   ⚠ {grp}: 无 KM 数据，跳过')
        continue

    n_raw = len(sub_wt)

    t = sub_km['time'].values
    s = sub_km['surv'].values
    ci_lo = sub_km['surv_lower'].values
    ci_hi = sub_km['surv_upper'].values

    color = GROUP_COLORS[grp]
    ls = LINESTYLES[i]
    # 图例格式与 01_km_all_before_psm 一致：仅 (n=...)
    label = f'{GROUP_LABELS[grp]} (n={n_raw:,})'

    ax_km.step(t, s, where='post', color=color, linewidth=1.8,
               linestyle=ls, label=label)
    ax_km.fill_between(t, ci_lo, ci_hi, step='post', alpha=0.08, color=color)

    group_info.append({
        'grp': grp, 'color': color, 'n_raw': n_raw,
        'time': t, 'surv': s, 'wt_df': sub_wt,
    })
    n_groups += 1
    w = sub_wt['iptw_weight'].values
    ess = w.sum()**2 / (w**2).sum() if len(w) > 0 else 0
    print(f'   {grp}: n={n_raw}, ESS={ess:.0f}')

ax_km.set_xlim(0, XLIM)
ax_km.set_ylim(-0.02, 1.05)
ax_km.set_ylabel('Overall Survival Probability')
ax_km.set_title(
    f'Overall Survival — All Groups (After CBPS-IPTW, N={len(wt_df):,})',
    fontweight='bold', pad=8)
ax_km.set_xticks(RISK_TIMES)
ax_km.yaxis.set_major_locator(mticker.MultipleLocator(0.2))
ax_km.legend(loc='upper right', handlelength=2.0, fontsize=6.5)

# [TIER1] 全局检验注记
# 由于 Schoenfeld 全局 P 严重违背 PH 假设（χ²=870 on 6 df, P<<0.001），
# Cox Wald 的 HR 解释失效。改以 RMST(τ=36mo) 全局 χ² 作为主要疗效假设检验，
# 并显式标注 PH 违背状态供读者参考。
def _fmt_p(p):
    if pd.isna(p):
        return 'NA'
    return 'P < 0.001' if p < 0.001 else f'P = {p:.3f}'

ph_row    = gt_df[gt_df['Test'] == 'PH_global']
wald_row  = gt_df[gt_df['Test'] == 'Wald_global']
rmst_row  = gt_df[gt_df['Test'] == 'RMST_global_tau36']

lines = []
if len(rmst_row) > 0:
    lines.append(f"RMST (τ=36 mo) global {_fmt_p(float(rmst_row['P_value'].iloc[0]))}")
if len(ph_row) > 0:
    ph_p = float(ph_row['P_value'].iloc[0])
    lines.append(f"PH violated (Schoenfeld {_fmt_p(ph_p)})")
if len(wald_row) > 0:
    lines.append(f"Cox Wald {_fmt_p(float(wald_row['P_value'].iloc[0]))}  (HR not interpretable)")

if lines:
    ax_km.text(0.97, 0.05, '\n'.join(lines),
               transform=ax_km.transAxes,
               fontsize=6.8, va='bottom', ha='right',
               linespacing=1.35, bbox=STATS_BOX)

# 风险表布局与 01_km_all_before_psm 的 draw_km 一致
ax_risk.set_xlim(0, XLIM + 0.5)
ax_risk.set_ylim(-0.8, n_groups + 0.2)
ax_risk.axis('off')
ax_risk.text(XLIM / 2, n_groups - 0.1, 'No. at risk',
             ha='center', va='top', fontsize=9,
             color='#666666', style='italic')

for idx, gi in enumerate(group_info):
    y = n_groups - idx - 0.9
    sub_wt = gi['wt_df']
    for t_pt in RISK_TIMES:
        n_at_risk = int((sub_wt['os_months'] >= t_pt).sum())
        ax_risk.text(t_pt, y, str(n_at_risk), ha='center', va='center',
                     fontsize=8, color=gi['color'], fontweight='bold')

ax_risk.text(XLIM / 2, -0.65, 'Time (months)',
             ha='center', va='center', fontsize=11, color='#333333')

save_fig(fig, 'km_7groups_template_matched')

print(f'\n完成！输出: {FIG_DIR}/km_7groups_template_matched.pdf/.png')
