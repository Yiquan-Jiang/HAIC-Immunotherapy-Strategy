#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
[TIER1] PH 假设违背时的主要疗效图

输出:
  figures/km/rmst_forest.{pdf,png}        — RMSTD forest (τ=24/36/60 vs HAIC_alone)
  figures/km/landmark_survival.{pdf,png}  — Landmark survival probabilities (12/24/36/48/60 mo)

数据来源 (来自 step3b_psm_vs_template.R):
  results/psm_vs_template/iptw_rmst.csv         per-group RMST + SE
  results/psm_vs_template/iptw_rmst_diff.csv    pairwise RMSTD vs HAIC_alone
  results/psm_vs_template/iptw_landmark.csv     per-group landmark surv + log-log CI
  results/psm_vs_template/iptw_global_test.csv  global χ²
"""

import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd

BASE_DIR = "/Users/yqj/Nutstore Files/我的坚果云/Liver_tumor_big_data/FIRST_LINE_HAIC_2025-12-30/HAIC_NO_TACE_4_TIDY/update_group_7"
EIGHT_GROUP = os.environ.get("EIGHT_GROUP", "0") == "1"
SFX = "_8group" if EIGHT_GROUP else ""
DATA_CSV = "analysis_ready_8group.csv" if EIGHT_GROUP else "analysis_ready.csv"
FIG_DIR  = os.path.join(BASE_DIR, 'figures', 'km' + SFX)
RES_DIR  = os.path.join(BASE_DIR, 'results', 'psm_vs_template' + SFX)
os.makedirs(FIG_DIR, exist_ok=True)

REF_GROUP = 'HAIC_alone'

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
    'lines.linewidth':     1.5,
    'legend.frameon':      False,
    'savefig.dpi':        300,
    'savefig.bbox':       'tight',
    'savefig.pad_inches':  0.05,
    'pdf.fonttype':        42,
    'ps.fonttype':         42,
})


def save_fig(fig, name):
    base = os.path.join(FIG_DIR, name)
    fig.savefig(f'{base}.pdf', bbox_inches='tight', pad_inches=0.05)
    fig.savefig(f'{base}.png', dpi=300, bbox_inches='tight', pad_inches=0.05)
    plt.close(fig)
    print(f'  Saved: {name}.pdf/.png')


def _fmt_p(p):
    if pd.isna(p):
        return 'NA'
    return '<0.001' if p < 0.001 else f'{p:.3f}'


# ════════════════════════════════════════════════════════════════════
# 读取数据
# ════════════════════════════════════════════════════════════════════
print('1. 读取 RMST / Landmark 数据...')
rmst_df     = pd.read_csv(os.path.join(RES_DIR, 'iptw_rmst.csv'))
rmstd_df    = pd.read_csv(os.path.join(RES_DIR, 'iptw_rmst_diff.csv'))
landmark_df = pd.read_csv(os.path.join(RES_DIR, 'iptw_landmark.csv'))
global_df   = pd.read_csv(os.path.join(RES_DIR, 'iptw_global_test.csv'))

available_taus = sorted(rmst_df['tau'].unique().tolist())
print(f'   RMST τ values: {available_taus}')
print(f'   Landmark times: {sorted(landmark_df["time"].unique().tolist())}')

# ════════════════════════════════════════════════════════════════════
# Figure 1: RMSTD forest (one panel per τ)
# ════════════════════════════════════════════════════════════════════
print('\n2. 绘制 RMSTD forest plot...')

non_ref_groups = [g for g in GROUP_ORDER if g != REF_GROUP]
y_pos = np.arange(len(non_ref_groups))[::-1]   # 顶部 = HAIC+I_concurrent

n_panels = len(available_taus)

# 布局: 紧凑型 — 左侧 1 列 group label + 每个 τ 给 (forest, table) 两个子列
LABEL_W   = 2.5
FOREST_W  = 1.3
TABLE_W   = 2.0
fig_w     = LABEL_W + n_panels * (FOREST_W + TABLE_W)
fig = plt.figure(figsize=(fig_w, 3.4))
width_ratios = [LABEL_W] + sum([[FOREST_W, TABLE_W]] * n_panels, [])
gs  = fig.add_gridspec(1, len(width_ratios),
                       width_ratios=width_ratios,
                       wspace=0.08,
                       left=0.01, right=0.995,
                       top=0.82, bottom=0.16)

# 左：组名列
ax_lab = fig.add_subplot(gs[0])
ax_lab.set_xlim(0, 1); ax_lab.set_ylim(-0.55, len(non_ref_groups) - 0.4)
ax_lab.axis('off')
ax_lab.text(0.02, len(non_ref_groups) - 0.15, 'Comparison vs HAIC alone',
            ha='left', va='bottom', fontsize=8.6, fontweight='bold')
for y, g in zip(y_pos, non_ref_groups):
    ax_lab.text(0.02, y, GROUP_LABELS[g], ha='left', va='center',
                fontsize=7.8, color=GROUP_COLORS[g], fontweight='bold')

# 计算所有 τ 共用的 forest x 轴范围（保证可视比较）
all_lo = rmstd_df['ci_lo'].min()
all_hi = rmstd_df['ci_hi'].max()
margin = 0.08 * (all_hi - all_lo)
x_lo   = min(0, all_lo - margin)
x_hi   = all_hi + margin

for i, tau in enumerate(available_taus):
    sub = rmstd_df[rmstd_df['tau'] == tau]
    sub = sub.set_index('group').reindex(non_ref_groups).reset_index()

    # ─── Forest 子列 ───
    ax = fig.add_subplot(gs[1 + 2 * i])
    ax.axvline(0, color='#888888', lw=0.7, ls='--', zorder=1)
    for y, (_, row) in zip(y_pos, sub.iterrows()):
        c = GROUP_COLORS[row['group']]
        ax.errorbar(row['rmstd'], y,
                    xerr=[[row['rmstd'] - row['ci_lo']],
                          [row['ci_hi'] - row['rmstd']]],
                    fmt='o', color=c, ecolor=c, elinewidth=1.2,
                    markersize=4.5, capsize=2.5, capthick=0.9, zorder=3)
    ax.set_xlim(x_lo, x_hi)
    ax.set_ylim(-0.55, len(non_ref_groups) - 0.4)
    ax.set_yticks([])
    ax.xaxis.set_major_locator(mticker.MaxNLocator(integer=False, nbins=4))
    ax.tick_params(axis='x', labelsize=7.5, pad=2)
    ax.set_xlabel('RMSTD (mo)', fontsize=8.0, labelpad=2)

    # 标题不放在 ax 上 (会超出 sub-axis)，改用 fig.text 跨越 forest+table
    pass

    # 移除 favours 文本（紧凑布局下空间不足；方向由零线 + 主标题已表达）

    # ─── Table 子列 (single combined row text, mono) ───
    ax_t = fig.add_subplot(gs[2 + 2 * i])
    ax_t.set_xlim(0, 1); ax_t.set_ylim(-0.55, len(non_ref_groups) - 0.4)
    ax_t.axis('off')
    for y, (_, row) in zip(y_pos, sub.iterrows()):
        p_str = _fmt_p(row['P_holm'])
        if not p_str.startswith('<'):
            p_str = '=' + p_str
        line = (f"{row['rmstd']:+5.2f} ({row['ci_lo']:+5.2f}, {row['ci_hi']:+5.2f})"
                f"  {p_str:>7s}")
        ax_t.text(0.0, y, line, ha='left', va='center', fontsize=7.0,
                  color='#222222', family='DejaVu Sans Mono')

# 在 forest+table 列对的中部写 τ 横幅（跨 axis 边界，避免标题挤压相邻 panel）
fig.canvas.draw()
for i, tau in enumerate(available_taus):
    ax_f = fig.axes[1 + 2 * i + 1]   # +1 因为第 0 个是 ax_lab；这里其实 forest 是 1+2i, table 是 2+2i
    # 重新拿 forest 与 table 实际 axes
    ax_forest = None; ax_table = None
    cnt = 0
    for k, a in enumerate(fig.axes):
        if k == 0: continue
        # axes 顺序: lab, forest_0, table_0, forest_1, table_1, forest_2, table_2
        pass
    ax_forest = fig.axes[1 + 2 * i]
    ax_table  = fig.axes[2 + 2 * i]
    bb_f = ax_forest.get_position()
    bb_t = ax_table.get_position()
    cx = 0.5 * (bb_f.x0 + bb_t.x1)
    y_top = max(bb_f.y1, bb_t.y1)

    grow = global_df[global_df['Test'] == f'RMST_global_tau{tau}']
    if len(grow):
        banner = (f'τ = {tau} mo   |   global χ² = {float(grow["Statistic"].iloc[0]):.0f},'
                  f' P {_fmt_p(float(grow["P_value"].iloc[0]))}')
    else:
        banner = f'τ = {tau} mo'
    fig.text(cx, y_top + 0.075, banner,
             ha='center', va='bottom', fontsize=8.8, fontweight='bold')
    cx_table = 0.5 * (bb_t.x0 + bb_t.x1)
    fig.text(cx_table, y_top + 0.020, 'RMSTD (95% CI)        P (Holm)',
             ha='center', va='bottom', fontsize=7.0, fontweight='bold',
             color='#444444', family='DejaVu Sans Mono')

fig.suptitle('Restricted Mean Survival Time Difference vs HAIC alone (CBPS-IPTW)',
             fontsize=10.5, fontweight='bold', y=0.99)
save_fig(fig, 'rmst_forest')

# ════════════════════════════════════════════════════════════════════
# Figure 2: Landmark survival probabilities
# ════════════════════════════════════════════════════════════════════
print('\n3. 绘制 landmark survival 图...')

landmark_times = sorted(landmark_df['time'].unique().tolist())
n_groups = len(GROUP_ORDER)
n_t      = len(landmark_times)

fig = plt.figure(figsize=(7.5, 5.0))
ax  = fig.add_subplot(111)

bar_w   = 0.8 / n_groups
x_base  = np.arange(n_t)

for gi, g in enumerate(GROUP_ORDER):
    sub = landmark_df[landmark_df['group'] == g].set_index('time').reindex(landmark_times)
    s   = sub['surv'].values
    ci_lo = sub['ci_lo'].values
    ci_hi = sub['ci_hi'].values
    yerr = np.vstack([s - ci_lo, ci_hi - s])
    yerr = np.where(yerr < 0, 0, yerr)
    xs = x_base + (gi - n_groups / 2 + 0.5) * bar_w
    ax.bar(xs, s, width=bar_w * 0.92,
           color=GROUP_COLORS[g], edgecolor='white', linewidth=0.4,
           label=GROUP_LABELS[g], zorder=2)
    ax.errorbar(xs, s, yerr=yerr, fmt='none',
                ecolor='#333333', elinewidth=0.6, capsize=1.6, zorder=3)

ax.set_xticks(x_base)
ax.set_xticklabels([f'{t} mo' for t in landmark_times])
ax.set_ylim(0, 1.02)
ax.yaxis.set_major_locator(mticker.MultipleLocator(0.2))
ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f'{v:.0%}'))
ax.set_ylabel('Survival probability (95% CI)')
ax.set_title('Landmark Overall Survival by Treatment Strategy (CBPS-IPTW, N=3,885)',
             fontsize=11, fontweight='bold', pad=6)
ax.legend(loc='upper right', ncol=2, fontsize=7.5,
          handlelength=1.2, columnspacing=0.8, labelspacing=0.25)
ax.grid(axis='y', linestyle=':', linewidth=0.5, color='#CCCCCC', zorder=0)

save_fig(fig, 'landmark_survival')

print(f'\n完成！输出: {FIG_DIR}/rmst_forest.pdf/png 与 landmark_survival.pdf/png')
