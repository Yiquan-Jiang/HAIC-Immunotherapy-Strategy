#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
OW 加权后基线平衡表 — 发表级图片渲染
每个亚组生成一张表格图（未加权 vs OW 加权，含 SMD 热图色条）
支持两组对比: I (HAIC+I) 和 IT (HAIC+I+T)
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import os

BASE_DIR = "/Users/yqj/Nutstore Files/我的坚果云/Liver_tumor_big_data/FIRST_LINE_HAIC_2025-12-30/HAIC_NO_TACE_4_TIDY/update_group_7"
RES_DIR  = os.path.join(BASE_DIR, 'results', 'subgroup_analysis')
FIG_DIR  = os.path.join(BASE_DIR, 'figures', 'subgroup')
os.makedirs(FIG_DIR, exist_ok=True)

plt.rcParams.update({
    'font.family':     'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'pdf.fonttype':    42,
    'ps.fonttype':     42,
})

# NPG 配色
C_BLUE  = '#3C5488'
C_CYAN  = '#4DBBD5'
C_RED   = '#E64B35'
C_TEAL  = '#00A087'
C_GRAY  = '#8491B4'

# SMD 颜色映射：0→绿，0.4→黄，1.0→红（输入已归一化到[0,1]）
SMD_CMAP = mcolors.LinearSegmentedColormap.from_list(
    'smd', [(0.0, '#27AE60'), (0.4, '#F39C12'), (1.0, '#E74C3C')], N=256)

# --- 对比组配置 ---
COMPARISONS = [
    {
        'tag': 'I',
        'treat_label': 'HAIC+I (concurrent)',
        'control_label': 'HAIC\u2192I (delayed)',
        'title_short': 'HAIC+Immuno Concurrent vs Delayed',
    },
    {
        'tag': 'IT',
        'treat_label': 'HAIC+I+T (concurrent)',
        'control_label': 'HAIC\u2192I+T (delayed)',
        'title_short': 'HAIC+Immuno+TKI Concurrent vs Delayed',
    },
]

for COMP in COMPARISONS:
    CMP_TAG     = COMP['tag']
    TREAT_LABEL = COMP['treat_label']
    CTRL_LABEL  = COMP['control_label']
    TITLE_SHORT = COMP['title_short']

    CMP_RES_DIR = os.path.join(RES_DIR, CMP_TAG)
    CMP_FIG_DIR = os.path.join(FIG_DIR, CMP_TAG)
    os.makedirs(CMP_FIG_DIR, exist_ok=True)

    print(f'\n{"#"*70}')
    print(f'#  对比组: {TITLE_SHORT}  (tag={CMP_TAG})')
    print(f'{"#"*70}')

    df_full = pd.read_csv(os.path.join(CMP_RES_DIR, f'ow_balance_table_full_{CMP_TAG}.csv'))

    SUBGROUPS = df_full['Subgroup'].unique().tolist()

    # ── 渲染单个亚组表格 ──────────────────────────────────────────
    def render_balance_table(sg_name, df_sg, out_path):
        n_vars = len(df_sg)
        n_conc = int(df_sg['N_concurrent'].iloc[0])
        n_then = int(df_sg['N_then_I'].iloc[0])

        # 列定义
        COL_WIDTHS = [3.2, 1.4, 1.4, 0.55, 1.4, 1.4, 0.55]
        COL_LABELS = [
            'Variable',
            f'Concurrent\n(n={n_conc})',
            f'Delayed\n(n={n_then})',
            'SMD',
            f'Concurrent\n(n={n_conc})',
            f'Delayed\n(n={n_then})',
            'SMD',
        ]
        total_w = sum(COL_WIDTHS)
        fig_w   = total_w + 0.4
        row_h   = 0.28
        header_h = 0.70
        fig_h   = header_h + n_vars * row_h + 0.35

        fig, ax = plt.subplots(figsize=(fig_w, fig_h))
        ax.set_xlim(0, total_w)
        ax.set_ylim(0, fig_h)
        ax.axis('off')

        # ── 表头背景 ──────────────────────────────────────────────
        ax.add_patch(mpatches.FancyBboxPatch(
            (0, fig_h - header_h), total_w, header_h,
            boxstyle='square,pad=0', facecolor='#2C3E50', edgecolor='none'))

        # 分组标题
        x_unw_mid = sum(COL_WIDTHS[:1]) + sum(COL_WIDTHS[1:4]) / 2
        x_ow_mid  = sum(COL_WIDTHS[:4]) + sum(COL_WIDTHS[4:]) / 2
        ax.text(x_unw_mid, fig_h - 0.18, 'Before OW',
                ha='center', va='center', fontsize=9, fontweight='bold',
                color='white')
        ax.text(x_ow_mid, fig_h - 0.18, 'After OW',
                ha='center', va='center', fontsize=9, fontweight='bold',
                color='white')

        # 分隔线
        ax.plot([sum(COL_WIDTHS[:1]), sum(COL_WIDTHS[:1])],
                [fig_h - header_h, fig_h], color='#566573', linewidth=0.5)
        ax.plot([sum(COL_WIDTHS[:4]), sum(COL_WIDTHS[:4])],
                [fig_h - header_h, fig_h], color='#566573', linewidth=0.5)
        ax.plot([sum(COL_WIDTHS[:1]), total_w],
                [fig_h - 0.32, fig_h - 0.32], color='#566573', linewidth=0.3)

        # 列标题
        x_pos = 0
        for ci, (cw, cl) in enumerate(zip(COL_WIDTHS, COL_LABELS)):
            ha = 'left' if ci == 0 else 'center'
            pad = 0.08 if ci == 0 else 0
            ax.text(x_pos + pad + (0 if ci == 0 else cw / 2),
                    fig_h - 0.52, cl,
                    ha=ha, va='center', fontsize=7.5, fontweight='bold',
                    color='white', linespacing=1.3)
            x_pos += cw

        # ── 数据行 ────────────────────────────────────────────────
        for row_idx, (_, row) in enumerate(df_sg.iterrows()):
            y_top = fig_h - header_h - row_idx * row_h
            y_mid = y_top - row_h / 2

            # 交替背景
            bg_color = '#F8F9FA' if row_idx % 2 == 0 else 'white'
            ax.add_patch(mpatches.FancyBboxPatch(
                (0, y_top - row_h), total_w, row_h,
                boxstyle='square,pad=0', facecolor=bg_color, edgecolor='none'))

            # 变量名
            ax.text(0.08, y_mid, row['Variable'],
                    ha='left', va='center', fontsize=7.5, color='#1A1A2E')

            # 数据列
            vals = [
                row['Concurrent_unw'], row['Then_I_unw'], row['SMD_unweighted'],
                row['Concurrent_OW'],  row['Then_I_OW'],  row['SMD_OW'],
            ]
            x_pos = COL_WIDTHS[0]
            for ci, (cw, val) in enumerate(zip(COL_WIDTHS[1:], vals)):
                x_mid = x_pos + cw / 2
                is_smd = (ci == 2 or ci == 5)

                if is_smd:
                    try:
                        smd_val = float(val)
                    except:
                        smd_val = np.nan

                    if not np.isnan(smd_val):
                        # SMD 彩色方块
                        norm_v = min(smd_val / 0.25, 1.0)
                        smd_color = SMD_CMAP(norm_v)
                        bar_w = cw * 0.55
                        bar_h = row_h * 0.45
                        ax.add_patch(mpatches.FancyBboxPatch(
                            (x_mid - bar_w / 2, y_mid - bar_h / 2),
                            bar_w, bar_h,
                            boxstyle='round,pad=0.01',
                            facecolor=smd_color, edgecolor='none', alpha=0.85))
                        txt_color = 'white' if norm_v > 0.5 else '#1A1A2E'
                        ax.text(x_mid, y_mid, f'{smd_val:.2f}',
                                ha='center', va='center', fontsize=6.5,
                                color=txt_color, fontweight='bold')
                    else:
                        ax.text(x_mid, y_mid, '—',
                                ha='center', va='center', fontsize=7, color='#999')
                else:
                    # 数值文字
                    is_ow_col = (ci >= 3)
                    txt_color = C_BLUE if is_ow_col else '#555555'
                    fw = 'bold' if is_ow_col else 'normal'
                    ax.text(x_mid, y_mid, str(val),
                            ha='center', va='center', fontsize=7,
                            color=txt_color, fontweight=fw)
                x_pos += cw

            # 行底线
            ax.plot([0, total_w], [y_top - row_h, y_top - row_h],
                    color='#E8E8E8', linewidth=0.3)

        # ── 垂直分隔线 ────────────────────────────────────────────
        x_pos = COL_WIDTHS[0]
        for ci, cw in enumerate(COL_WIDTHS[1:]):
            lc = '#CCCCCC' if ci not in (2, 5) else '#AAAAAA'
            lw = 0.4 if ci not in (2, 5) else 0.6
            ax.plot([x_pos, x_pos], [0.2, fig_h - header_h],
                    color=lc, linewidth=lw)
            x_pos += cw

        # ── 底部框线 ──────────────────────────────────────────────
        ax.plot([0, total_w], [0.2, 0.2], color='#2C3E50', linewidth=0.8)
        ax.plot([0, total_w], [fig_h - header_h, fig_h - header_h],
                color='#2C3E50', linewidth=0.8)

        # ── SMD 图例 ──────────────────────────────────────────────
        legend_x = 0.08
        legend_y = 0.06
        for val, label in [(0.0, '0'), (0.05, '0.05'), (0.10, '0.10'), (0.20, '≥0.20')]:
            norm_v = min(val / 0.25, 1.0)
            c = SMD_CMAP(norm_v)
            ax.add_patch(mpatches.FancyBboxPatch(
                (legend_x, legend_y - 0.04), 0.18, 0.08,
                boxstyle='round,pad=0.01', facecolor=c, edgecolor='none', alpha=0.85))
            ax.text(legend_x + 0.09, legend_y + 0.06, label,
                    ha='center', va='bottom', fontsize=6, color='#444')
            legend_x += 0.32

        ax.text(0.08, legend_y + 0.10, 'SMD scale:',
                ha='left', va='bottom', fontsize=6.5, color='#666', style='italic')
        ax.text(total_w - 0.08, legend_y + 0.06,
                'Bold blue = OW-adjusted values',
                ha='right', va='bottom', fontsize=6.5, color=C_BLUE, style='italic')

        # ── 标题 ──────────────────────────────────────────────────
        ax.set_title(f'Baseline Characteristics — {sg_name}\n'
                     f'{TREAT_LABEL} vs {CTRL_LABEL}',
                     fontsize=9.5, fontweight='bold', pad=8, loc='left',
                     color='#1A1A2E')

        plt.tight_layout(pad=0.3)
        fig.savefig(f'{out_path}.pdf', bbox_inches='tight', pad_inches=0.08)
        fig.savefig(f'{out_path}.png', dpi=300, bbox_inches='tight', pad_inches=0.08)
        plt.close(fig)
        print(f'   已保存: {os.path.basename(out_path)}.pdf/.png')


    # ── 逐亚组渲染 ────────────────────────────────────────────────
    print('渲染 OW 加权基线平衡表...\n')
    for sg_name in SUBGROUPS:
        df_sg = df_full[df_full['Subgroup'] == sg_name].copy()
        safe  = sg_name.replace(' ', '_').replace('>', 'gt').replace('/', '_')
        out   = os.path.join(CMP_FIG_DIR, f'ow_balance_table_{safe}')
        render_balance_table(sg_name, df_sg, out)

    # ── 合并图（所有亚组纵向拼接）────────────────────────────────
    print('\n生成合并总表...')

    n_sg   = len(SUBGROUPS)
    n_vars = len(df_full[df_full['Subgroup'] == SUBGROUPS[0]])

    COL_WIDTHS = [3.0, 1.3, 1.3, 0.5, 1.3, 1.3, 0.5]
    total_w    = sum(COL_WIDTHS)
    row_h      = 0.26
    header_h   = 0.65
    sg_gap     = 0.30
    fig_w      = total_w + 0.4
    fig_h      = (header_h + n_vars * row_h + sg_gap) * n_sg + 0.5

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.set_xlim(0, total_w)
    ax.set_ylim(0, fig_h)
    ax.axis('off')

    y_cursor = fig_h

    for sg_idx, sg_name in enumerate(SUBGROUPS):
        df_sg   = df_full[df_full['Subgroup'] == sg_name].copy()
        n_conc  = int(df_sg['N_concurrent'].iloc[0])
        n_then  = int(df_sg['N_then_I'].iloc[0])
        n_vars_sg = len(df_sg)

        block_h = header_h + n_vars_sg * row_h

        # 亚组标题栏
        ax.add_patch(mpatches.FancyBboxPatch(
            (0, y_cursor - header_h), total_w, header_h,
            boxstyle='square,pad=0', facecolor='#2C3E50', edgecolor='none'))

        # 分组标题
        x_unw_mid = COL_WIDTHS[0] + sum(COL_WIDTHS[1:4]) / 2
        x_ow_mid  = sum(COL_WIDTHS[:4]) + sum(COL_WIDTHS[4:]) / 2
        ax.text(x_unw_mid, y_cursor - 0.16, 'Before OW',
                ha='center', va='center', fontsize=8.5, fontweight='bold', color='white')
        ax.text(x_ow_mid, y_cursor - 0.16, 'After OW',
                ha='center', va='center', fontsize=8.5, fontweight='bold', color='white')

        ax.plot([COL_WIDTHS[0], COL_WIDTHS[0]],
                [y_cursor - header_h, y_cursor], color='#566573', linewidth=0.5)
        ax.plot([sum(COL_WIDTHS[:4]), sum(COL_WIDTHS[:4])],
                [y_cursor - header_h, y_cursor], color='#566573', linewidth=0.5)
        ax.plot([COL_WIDTHS[0], total_w],
                [y_cursor - 0.30, y_cursor - 0.30], color='#566573', linewidth=0.3)

        COL_LABELS = [
            sg_name,
            f'Concurrent\n(n={n_conc})',
            f'Delayed\n(n={n_then})',
            'SMD',
            f'Concurrent\n(n={n_conc})',
            f'Delayed\n(n={n_then})',
            'SMD',
        ]
        x_pos = 0
        for ci, (cw, cl) in enumerate(zip(COL_WIDTHS, COL_LABELS)):
            ha = 'left' if ci == 0 else 'center'
            pad = 0.08 if ci == 0 else 0
            ax.text(x_pos + pad + (0 if ci == 0 else cw / 2),
                    y_cursor - 0.50, cl,
                    ha=ha, va='center', fontsize=7, fontweight='bold',
                    color='white', linespacing=1.3)
            x_pos += cw

        # 数据行
        for row_idx, (_, row) in enumerate(df_sg.iterrows()):
            y_top = y_cursor - header_h - row_idx * row_h
            y_mid = y_top - row_h / 2

            bg = '#F8F9FA' if row_idx % 2 == 0 else 'white'
            ax.add_patch(mpatches.FancyBboxPatch(
                (0, y_top - row_h), total_w, row_h,
                boxstyle='square,pad=0', facecolor=bg, edgecolor='none'))

            ax.text(0.07, y_mid, row['Variable'],
                    ha='left', va='center', fontsize=7, color='#1A1A2E')

            vals = [
                row['Concurrent_unw'], row['Then_I_unw'], row['SMD_unweighted'],
                row['Concurrent_OW'],  row['Then_I_OW'],  row['SMD_OW'],
            ]
            x_pos = COL_WIDTHS[0]
            for ci, (cw, val) in enumerate(zip(COL_WIDTHS[1:], vals)):
                x_mid = x_pos + cw / 2
                is_smd = (ci == 2 or ci == 5)
                if is_smd:
                    try:
                        smd_val = float(val)
                    except:
                        smd_val = np.nan
                    if not np.isnan(smd_val):
                        norm_v = min(smd_val / 0.25, 1.0)
                        c = SMD_CMAP(norm_v)
                        bw = cw * 0.55; bh = row_h * 0.45
                        ax.add_patch(mpatches.FancyBboxPatch(
                            (x_mid - bw / 2, y_mid - bh / 2), bw, bh,
                            boxstyle='round,pad=0.01', facecolor=c,
                            edgecolor='none', alpha=0.85))
                        tc = 'white' if norm_v > 0.5 else '#1A1A2E'
                        ax.text(x_mid, y_mid, f'{smd_val:.2f}',
                                ha='center', va='center', fontsize=6,
                                color=tc, fontweight='bold')
                    else:
                        ax.text(x_mid, y_mid, '—',
                                ha='center', va='center', fontsize=6.5, color='#999')
                else:
                    is_ow = (ci >= 3)
                    ax.text(x_mid, y_mid, str(val),
                            ha='center', va='center', fontsize=6.5,
                            color=C_BLUE if is_ow else '#555',
                            fontweight='bold' if is_ow else 'normal')
                x_pos += cw

            ax.plot([0, total_w], [y_top - row_h, y_top - row_h],
                    color='#E8E8E8', linewidth=0.25)

        # 框线
        ax.plot([0, total_w], [y_cursor - header_h, y_cursor - header_h],
                color='#2C3E50', linewidth=0.7)
        ax.plot([0, total_w],
                [y_cursor - header_h - n_vars_sg * row_h,
                 y_cursor - header_h - n_vars_sg * row_h],
                color='#AAAAAA', linewidth=0.5)

        y_cursor -= (block_h + sg_gap)

    # 底线
    ax.plot([0, total_w], [y_cursor + sg_gap * 0.5, y_cursor + sg_gap * 0.5],
            color='#2C3E50', linewidth=0.8)

    # SMD 图例
    legend_x = 0.08
    legend_y = y_cursor + sg_gap * 0.15
    for val, label in [(0.0, '0'), (0.05, '0.05'), (0.10, '0.10'), (0.20, '≥0.20')]:
        norm_v = min(val / 0.25, 1.0)
        c = SMD_CMAP(norm_v)
        ax.add_patch(mpatches.FancyBboxPatch(
            (legend_x, legend_y - 0.04), 0.16, 0.07,
            boxstyle='round,pad=0.01', facecolor=c, edgecolor='none', alpha=0.85))
        ax.text(legend_x + 0.08, legend_y + 0.05, label,
                ha='center', va='bottom', fontsize=5.5, color='#444')
        legend_x += 0.30

    ax.text(0.08, legend_y + 0.12, 'SMD:',
            ha='left', va='bottom', fontsize=6, color='#666', style='italic')
    ax.text(total_w - 0.08, legend_y + 0.12,
            'Bold blue = OW-adjusted  |  SMD <0.10 indicates good balance',
            ha='right', va='bottom', fontsize=6, color='#555', style='italic')

    ax.set_title(
        'Baseline Characteristics Before and After Overlap Weighting\n'
        f'{TREAT_LABEL} vs {CTRL_LABEL}',
        fontsize=10, fontweight='bold', pad=10, loc='left', color='#1A1A2E')

    plt.tight_layout(pad=0.3)
    out_combined = os.path.join(CMP_FIG_DIR, f'ow_balance_table_combined_{CMP_TAG}')
    fig.savefig(f'{out_combined}.pdf', bbox_inches='tight', pad_inches=0.1)
    fig.savefig(f'{out_combined}.png', dpi=300, bbox_inches='tight', pad_inches=0.1)
    plt.close(fig)
    print(f'\n合并总表已保存: {out_combined}.pdf/.png')
    print('\n✅ 全部完成！')
