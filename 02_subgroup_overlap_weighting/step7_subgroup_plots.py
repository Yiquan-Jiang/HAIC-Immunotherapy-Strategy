#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
亚组分析可视化 — 发表级森林图 + KM 曲线
读取 R 端 OW 分析结果，绘制 Nature/Lancet 级图表
支持两组对比: I (HAIC+I) 和 IT (HAIC+I+T)
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
import os
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test

BASE_DIR = "/Users/yqj/Nutstore Files/我的坚果云/Liver_tumor_big_data/FIRST_LINE_HAIC_2025-12-30/HAIC_NO_TACE_4_TIDY/update_group_7"
RES_DIR  = os.path.join(BASE_DIR, 'results', 'subgroup_analysis')
FIG_DIR  = os.path.join(BASE_DIR, 'figures', 'subgroup')
DATA_DIR = os.path.join(BASE_DIR, 'data')
os.makedirs(FIG_DIR, exist_ok=True)

# ════════════════════════════════════════════════════════════════
# 全局 rcParams — 发表级规范
# ════════════════════════════════════════════════════════════════
plt.rcParams.update({
    'font.family':        'sans-serif',
    'font.sans-serif':    ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size':           8,
    'axes.labelsize':      9,
    'axes.titlesize':     10,
    'xtick.labelsize':     7,
    'ytick.labelsize':     7,
    'legend.fontsize':     7,
    'axes.linewidth':      0.6,
    'axes.spines.top':     False,
    'axes.spines.right':   False,
    'xtick.major.width':   0.6,
    'ytick.major.width':   0.6,
    'xtick.major.size':    3.0,
    'ytick.major.size':    3.0,
    'lines.linewidth':     1.5,
    'lines.markersize':    4,
    'legend.frameon':      False,
    'legend.borderpad':    0.3,
    'legend.handlelength': 1.5,
    'figure.dpi':         150,
    'savefig.dpi':        300,
    'savefig.bbox':       'tight',
    'savefig.pad_inches':  0.05,
    'pdf.fonttype':        42,
    'ps.fonttype':         42,
    'axes.grid':           False,
})

# NPG 配色
C_RED   = '#E64B35'
C_CYAN  = '#4DBBD5'
C_TEAL  = '#00A087'
C_BLUE  = '#3C5488'
C_CORAL = '#F39B7F'
C_GRAY  = '#8491B4'
C_DRED  = '#DC0000'

# --- 对比组配置 ---
COMPARISONS = [
    {
        'tag': 'I',
        'treat_group': 'HAIC+I_concurrent',
        'control_group': 'HAIC_then_I',
        'treat_label': 'HAIC+I (concurrent)',
        'control_label': 'HAIC\u2192I (delayed)',
        'title_short': 'HAIC+Immuno Concurrent vs Delayed',
    },
    {
        'tag': 'IT',
        'treat_group': 'HAIC+I+T_concurrent',
        'control_group': 'HAIC_then_I+T',
        'treat_label': 'HAIC+I+T (concurrent)',
        'control_label': 'HAIC\u2192I+T (delayed)',
        'title_short': 'HAIC+Immuno+TKI Concurrent vs Delayed',
    },
]

for COMP in COMPARISONS:
    CMP_TAG       = COMP['tag']
    TREAT_GROUP   = COMP['treat_group']
    CONTROL_GROUP = COMP['control_group']
    TREAT_LABEL   = COMP['treat_label']
    CTRL_LABEL    = COMP['control_label']
    TITLE_SHORT   = COMP['title_short']

    CMP_RES_DIR = os.path.join(RES_DIR, CMP_TAG)
    CMP_FIG_DIR = os.path.join(FIG_DIR, CMP_TAG)
    os.makedirs(CMP_FIG_DIR, exist_ok=True)

    print(f'\n{"#"*70}')
    print(f'#  对比组: {TREAT_GROUP} vs {CONTROL_GROUP}  (tag={CMP_TAG})')
    print(f'{"#"*70}')

    # ════════════════════════════════════════════════════════════════
    # 1. 读取数据
    # ════════════════════════════════════════════════════════════════
    print('1. 读取 OW 分析结果...')
    df_ow    = pd.read_csv(os.path.join(CMP_RES_DIR, f'ow_subgroup_results_{CMP_TAG}.csv'))
    df_inter = pd.read_csv(os.path.join(CMP_RES_DIR, f'ow_interaction_tests_{CMP_TAG}.csv'))

    df_all = pd.read_csv(os.path.join(DATA_DIR, 'analysis_ready.csv'))
    df_all = df_all[df_all['os_months'] >= 0].copy()
    df_pop = df_all[df_all['main_group'].isin([TREAT_GROUP, CONTROL_GROUP])].copy()
    df_pop['event'] = df_pop['death_status'].map(
        {'Yes': 1, 'No': 0, '1': 1, '0': 0, 1: 1, 0: 0}).fillna(0).astype(int)
    df_pop['treat'] = (df_pop['main_group'] == TREAT_GROUP).astype(int)

    cph_all = CoxPHFitter()
    cph_all.fit(df_pop[['os_months', 'event', 'treat']],
                duration_col='os_months', event_col='event')
    hr_all = float(np.exp(cph_all.params_['treat']))
    ci_all = np.exp(cph_all.confidence_intervals_.loc['treat']).values
    p_all  = float(cph_all.summary.loc['treat', 'p'])
    n_conc = int(df_pop['treat'].sum())
    n_then = int((df_pop['treat'] == 0).sum())

    # PSM results — try to load (may not exist for all comparisons)
    has_psm = False
    try:
        psm_csv = os.path.join(BASE_DIR, 'results', 'psm_balance_tables_complete',
                               'survival_analysis_final.csv')
        df_psm  = pd.read_csv(psm_csv)
        psm_rows = df_psm[(df_psm['Group1'] == TREAT_GROUP) &
                          (df_psm['Group2'] == CONTROL_GROUP)]
        if len(psm_rows) > 0:
            row_psm = psm_rows.iloc[0]
            hr_psm    = 1.0 / float(row_psm['HR'])
            ci_psm_lo = 1.0 / float(row_psm['CI_upper'])
            ci_psm_hi = 1.0 / float(row_psm['CI_lower'])
            p_psm     = float(row_psm['P_value'])
            n_psm     = int(row_psm['N1_after']) + int(row_psm['N2_after'])
            has_psm = True
    except Exception:
        pass

    print(f'   Overall unadjusted: HR={hr_all:.3f} ({ci_all[0]:.3f}-{ci_all[1]:.3f}) P={p_all:.4f}')
    if has_psm:
        print(f'   Overall PSM:        HR={hr_psm:.3f} ({ci_psm_lo:.3f}-{ci_psm_hi:.3f}) P={p_psm:.4f}')
    else:
        print(f'   Overall PSM:        N/A (no PSM results for this comparison)')

    # ════════════════════════════════════════════════════════════════
    # 2. 发表级森林图
    # ════════════════════════════════════════════════════════════════
    print('\n2. 绘制发表级森林图...')

    ROWS = [
        # (label, csv_key, row_type, n_conc, n_then)
        # row_type: 'overall_unadj', 'overall_psm', 'sep', 'header', 'sub', 'sub_comp'
        ('Overall (unadjusted)',    '__unadj__',    'overall_unadj'),
    ]
    if has_psm:
        ROWS.append(('Overall (PSM matched)',   '__psm__',      'overall_psm'))
    ROWS += [
        ('SEP',                     '__sep1__',     'sep'),
        ('Composite high-risk',     'Composite high-risk', 'header'),
        ('SEP',                     '__sep2__',     'sep'),
        ('Tumor count',             '__hdr_tc__',   'section_hdr'),
        ('  >3',                    'Tumor count >3',            'sub'),
        ('  \u22643',               'Non-Tumor count >3',        'sub_comp'),
        ('Tumor diameter',          '__hdr_td__',   'section_hdr'),
        ('  >10 cm',                'Tumor diameter >10 cm',     'sub'),
        ('  \u226410 cm',           'Non-Tumor diameter >10 cm', 'sub_comp'),
        ('PVTT',                    '__hdr_pvtt__', 'section_hdr'),
        ('  Vp3/4',                 'PVTT Vp3/4',                'sub'),
        ('  Absent or Vp1\u20132',  'Non-PVTT Vp3/4',           'sub_comp'),
        ('Extrahepatic metastasis', '__hdr_em__',   'section_hdr'),
        ('  Yes',                   'Extrahepatic metastasis',   'sub'),
        ('  No',                    'Non-Extrahepatic metastasis','sub_comp'),
    ]

    EQUIV_LO, EQUIV_HI = 0.60, 1.67

    COLOR_MAP = {
        'overall_unadj': '#888888',
        'overall_psm':   C_TEAL,
        'header':        C_BLUE,
        'sub':           C_BLUE,
        'sub_comp':      C_BLUE,
    }

    ROW_H       = 0.42
    SEP_H       = 0.15
    SECTION_H   = 0.30
    FIG_W       = 7.2

    y_positions = []
    y = 0
    for r in ROWS:
        rtype = r[2]
        if rtype == 'sep':
            y -= SEP_H
            y_positions.append(y)
        elif rtype == 'section_hdr':
            y -= SECTION_H
            y_positions.append(y)
        else:
            y -= ROW_H
            y_positions.append(y)

    y_min = min(y_positions) - ROW_H * 0.6
    y_max = 0 + ROW_H * 0.3

    fig_h = (y_max - y_min) / 2.54 * 2.54 * 0.38
    fig_h = max(fig_h, 5.5)
    fig_h = min(fig_h, 7.5)

    fig, ax = plt.subplots(figsize=(FIG_W, 6.2))

    X_MIN, X_MAX = 0.20, 6.5
    ax.set_xlim(X_MIN, X_MAX)
    ax.set_xscale('log')
    ax.set_ylim(y_min, y_max)
    ax.set_yticks([])
    ax.spines['left'].set_visible(False)

    ax.axvline(x=1.0, color='#555555', linestyle='--', linewidth=0.8, alpha=0.6, zorder=1)

    ax.axvspan(EQUIV_LO, EQUIV_HI, alpha=0.04, color=C_TEAL, zorder=0)
    ax.axvline(x=EQUIV_LO, color=C_TEAL, linestyle=':', linewidth=0.5, alpha=0.4, zorder=1)
    ax.axvline(x=EQUIV_HI, color=C_TEAL, linestyle=':', linewidth=0.5, alpha=0.4, zorder=1)

    HEADER_Y = y_max - ROW_H * 0.05
    COL_LABEL = 0.01
    COL_HR    = X_MAX * 1.15
    COL_P     = X_MAX * 3.2
    COL_PINT  = X_MAX * 5.5

    ax.text(COL_LABEL, HEADER_Y, 'Subgroup', ha='left', va='bottom',
            fontsize=8.5, fontweight='bold', color='#111111',
            transform=ax.get_yaxis_transform())
    ax.text(COL_HR, HEADER_Y, 'HR (95% CI)', ha='left', va='bottom',
            fontsize=7.5, fontweight='bold', color='#111111')
    ax.text(COL_P, HEADER_Y, 'P', ha='center', va='bottom',
            fontsize=7.5, fontweight='bold', color='#111111')
    ax.text(COL_PINT, HEADER_Y, 'P int', ha='center', va='bottom',
            fontsize=7.5, fontweight='bold', color='#111111')

    ax.axhline(y=HEADER_Y - 0.06, color='#333333', linewidth=0.8,
               xmin=0, xmax=1, clip_on=False)

    def fmt_p(p):
        if p < 0.001:
            return '<0.001'
        return f'{p:.3f}'

    for idx, (label, key, rtype) in enumerate(ROWS):
        ypos = y_positions[idx]

        if rtype == 'sep':
            ax.axhline(y=ypos + SEP_H * 0.5, color='#E0E0E0', linewidth=0.4,
                       xmin=0, xmax=1, clip_on=False)
            continue

        if rtype == 'section_hdr':
            ax.text(COL_LABEL, ypos, label, ha='left', va='center',
                    fontsize=8, fontweight='bold', fontstyle='italic', color='#333333',
                    transform=ax.get_yaxis_transform())
            continue

        if rtype == 'overall_unadj':
            hr, ci_lo, ci_hi, p_val = hr_all, ci_all[0], ci_all[1], p_all
            n_total = n_conc + n_then
            p_int_str = ''
        elif rtype == 'overall_psm':
            hr, ci_lo, ci_hi, p_val = hr_psm, ci_psm_lo, ci_psm_hi, p_psm
            n_total = n_psm
            p_int_str = ''
        else:
            row = df_ow[df_ow['Subgroup'] == key]
            if row.empty:
                continue
            row = row.iloc[0]
            hr      = float(row['HR_OW'])
            ci_lo   = float(row['CI95_lo_OW'])
            ci_hi   = float(row['CI95_hi_OW'])
            p_val   = float(row['P_OW'])
            n_total = int(row['N_concurrent']) + int(row['N_then_I'])

            inter_row = df_inter[df_inter['Subgroup'] == key]
            if not inter_row.empty:
                p_int_str = fmt_p(float(inter_row.iloc[0]['P_interaction']))
            else:
                base_name = key.replace('Non-', '')
                inter_row2 = df_inter[df_inter['Subgroup'] == base_name]
                p_int_str = ''

        color = COLOR_MAP.get(rtype, C_BLUE)
        is_sig = p_val < 0.05
        dot_color = C_DRED if is_sig else color

        lbl_weight = 'bold' if rtype in ('overall_unadj', 'overall_psm', 'header') else 'normal'
        lbl_size   = 8.5 if rtype in ('overall_unadj', 'overall_psm', 'header') else 7.5
        lbl_x      = COL_LABEL if rtype in ('overall_unadj', 'overall_psm', 'header') else 0.04

        display_label = f'{label}  (n={n_total})'
        ax.text(lbl_x, ypos, display_label, ha='left', va='center',
                fontsize=lbl_size, fontweight=lbl_weight, color='#111111',
                transform=ax.get_yaxis_transform())

        ci_lo_plot = max(ci_lo, X_MIN * 1.02)
        ci_hi_plot = min(ci_hi, X_MAX * 0.98)

        ax.plot([ci_lo_plot, ci_hi_plot], [ypos, ypos],
                color=dot_color, linewidth=1.6, solid_capstyle='round', zorder=3)
        for xc in [ci_lo_plot, ci_hi_plot]:
            ax.plot([xc, xc], [ypos - ROW_H * 0.10, ypos + ROW_H * 0.10],
                    color=dot_color, linewidth=0.9, zorder=3)

        d_size = 55 if rtype in ('overall_unadj', 'overall_psm', 'header') else 38
        ax.scatter([hr], [ypos], marker='D', s=d_size,
                   color=dot_color, zorder=4, linewidths=0)

        hr_txt = f'{hr:.2f} ({ci_lo:.2f}\u2013{ci_hi:.2f})'
        txt_color = '#B03000' if is_sig else '#333333'
        ax.text(COL_HR, ypos, hr_txt, ha='left', va='center',
                fontsize=7, color=txt_color, fontfamily='monospace')
        ax.text(COL_P, ypos, fmt_p(p_val), ha='center', va='center',
                fontsize=7, color=txt_color, fontweight='bold' if is_sig else 'normal')
        if p_int_str:
            ax.text(COL_PINT, ypos, p_int_str, ha='center', va='center',
                    fontsize=7, color='#555555')

    ax.axhline(y=y_min + ROW_H * 0.3, color='#333333', linewidth=0.8,
               xmin=0, xmax=1, clip_on=False)

    ax.set_xlabel('Hazard Ratio (95% CI)', fontsize=9, labelpad=8)
    xticks = [0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0]
    ax.set_xticks(xticks)
    ax.set_xticklabels([str(x) for x in xticks], fontsize=7.5)

    ax.text(0.25, -0.04, '\u2190 Favors concurrent',
            ha='center', va='top', fontsize=7.5, color='#666666',
            transform=ax.transAxes)
    ax.text(0.65, -0.04, 'Favors delayed \u2192',
            ha='center', va='top', fontsize=7.5, color='#666666',
            transform=ax.transAxes)

    legend_elements = [
        mlines.Line2D([0], [0], color='#888888', linewidth=1.5,
                      marker='D', markersize=4.5, markerfacecolor='#888888',
                      label='Unadjusted'),
        mlines.Line2D([0], [0], color=C_TEAL, linewidth=1.5,
                      marker='D', markersize=4.5, markerfacecolor=C_TEAL,
                      label='PSM matched'),
        mlines.Line2D([0], [0], color=C_BLUE, linewidth=1.5,
                      marker='D', markersize=4.5, markerfacecolor=C_BLUE,
                      label='Overlap weighted'),
        mlines.Line2D([0], [0], color=C_DRED, linewidth=1.5,
                      marker='D', markersize=4.5, markerfacecolor=C_DRED,
                      label='P < 0.05'),
        mpatches.Patch(facecolor=C_TEAL, alpha=0.12, edgecolor=C_TEAL,
                       linewidth=0.5, label='Equivalence [0.60\u20131.67]'),
    ]
    ax.legend(handles=legend_elements, loc='lower center',
              bbox_to_anchor=(0.42, -0.12), ncol=5,
              fontsize=6.5, frameon=False, handlelength=1.8,
              columnspacing=0.8, handletextpad=0.4)

    plt.tight_layout()
    out_forest = os.path.join(CMP_FIG_DIR, f'ow_subgroup_forest_plot_{CMP_TAG}')
    fig.savefig(f'{out_forest}.pdf', bbox_inches='tight', pad_inches=0.08)
    fig.savefig(f'{out_forest}.png', dpi=300, bbox_inches='tight', pad_inches=0.08)
    plt.close(fig)
    print(f'   森林图已保存: {out_forest}.pdf/.png')

    # ════════════════════════════════════════════════════════════════
    # 3. 发表级 KM 曲线 (2×2 多面板 + 风险表)
    # ════════════════════════════════════════════════════════════════
    print('\n3. 绘制发表级 KM 曲线...')

    plt.rcParams.update({
        'font.size':          10,
        'axes.labelsize':     11,
        'axes.titlesize':     12,
        'xtick.labelsize':     9,
        'ytick.labelsize':     9,
        'legend.fontsize':     9,
        'axes.linewidth':      0.8,
        'xtick.major.width':   0.8,
        'ytick.major.width':   0.8,
        'xtick.major.size':    3.5,
        'ytick.major.size':    3.5,
        'lines.linewidth':     2.0,
    })

    RISK_TIMES = [0, 6, 12, 18, 24, 30, 36]
    STATS_BOX  = dict(boxstyle='round,pad=0.4', facecolor='white',
                      edgecolor='#CCCCCC', alpha=0.95, linewidth=0.5)

    KM_SUBGROUPS = [
        ('Composite high-risk', 'Composite_high-risk'),
        ('PVTT Vp3/4',          'PVTT_Vp3_4'),
        ('Extrahepatic metastasis', 'Extrahepatic_metastasis'),
        ('Tumor diameter >10 cm',   'Tumor_diameter__10_cm'),
    ]

    fig = plt.figure(figsize=(10.5, 11.0))
    outer_gs = fig.add_gridspec(2, 2, hspace=0.32, wspace=0.30)

    COLORS_KM = [C_CYAN, C_RED]
    LABELS_KM = {0: CTRL_LABEL, 1: TREAT_LABEL}
    LS_KM     = {0: '-', 1: '--'}

    for panel_idx, (sg_name, file_suffix) in enumerate(KM_SUBGROUPS):
        csv_path = os.path.join(CMP_RES_DIR, f'ow_weighted_ids_{file_suffix}.csv')
        if not os.path.exists(csv_path):
            print(f'   Warning: {csv_path} not found')
            continue

        sg_data = pd.read_csv(csv_path)

        inner_gs = outer_gs[panel_idx // 2, panel_idx % 2].subgridspec(
            2, 1, height_ratios=[4.5, 1.0], hspace=0.05)
        ax_km   = fig.add_subplot(inner_gs[0])
        ax_risk = fig.add_subplot(inner_gs[1])

        kmf_dict = {}
        for grp_val in [0, 1]:
            sub = sg_data[sg_data['treat'] == grp_val]
            kmf = KaplanMeierFitter()
            kmf.fit(sub['os_months'], sub['event'], label=LABELS_KM[grp_val])
            kmf_dict[grp_val] = kmf

            t     = kmf.survival_function_.index.values
            s     = kmf.survival_function_.iloc[:, 0].values
            ci_lo = kmf.confidence_interval_.iloc[:, 0].values
            ci_hi = kmf.confidence_interval_.iloc[:, 1].values
            n     = len(sub)

            ax_km.step(t, s, where='post', color=COLORS_KM[grp_val],
                       linewidth=2.0, linestyle=LS_KM[grp_val],
                       label=f'{LABELS_KM[grp_val]} (n={n})')
            ax_km.fill_between(t, ci_lo, ci_hi, step='post',
                               alpha=0.08, color=COLORS_KM[grp_val])

        max_t = sg_data['os_months'].quantile(0.95)
        xlim  = min(max(max_t, 24), 48)

        ax_km.set_xlim(0, xlim)
        ax_km.set_ylim(-0.02, 1.05)
        ax_km.set_ylabel('Overall Survival Probability')
        ax_km.yaxis.set_major_locator(mticker.MultipleLocator(0.2))
        ax_km.set_xticks(RISK_TIMES)
        ax_km.tick_params(labelbottom=False)
        ax_km.legend(loc='upper right', fontsize=8, handlelength=2.0, frameon=False)

        panel_label = chr(65 + panel_idx)
        ax_km.text(-0.12, 1.08, panel_label, transform=ax_km.transAxes,
                   fontsize=12, fontweight='bold', va='top')
        ax_km.set_title(sg_name, fontweight='bold', pad=8, fontsize=11)

        ow_row = df_ow[df_ow['Subgroup'] == sg_name]
        if not ow_row.empty:
            r = ow_row.iloc[0]
            hr_ow   = float(r['HR_OW'])
            ci95_lo = float(r['CI95_lo_OW'])
            ci95_hi = float(r['CI95_hi_OW'])
            p_ow    = float(r['P_OW'])
            med_c   = float(r['Median_OS_conc'])
            med_t   = float(r['Median_OS_then'])
            p_str   = 'P < 0.001' if p_ow < 0.001 else f'P = {p_ow:.3f}'
            stats_lines = [
                f'OW HR = {hr_ow:.2f} (95% CI {ci95_lo:.2f}\u2013{ci95_hi:.2f})',
                p_str,
                f'Median OS: {med_c:.1f} vs {med_t:.1f} mo',
            ]
            ax_km.text(0.97, 0.05, '\n'.join(stats_lines),
                       transform=ax_km.transAxes,
                       fontsize=8, va='bottom', ha='right', bbox=STATS_BOX,
                       linespacing=1.4)

        # ── 风险表 ──
        ax_risk.set_xlim(0, xlim)
        n_groups = 2
        ax_risk.set_ylim(-0.8, n_groups + 0.3)
        ax_risk.axis('off')

        ax_risk.text(xlim / 2, n_groups + 0.05, 'No. at risk',
                     ha='center', va='top', fontsize=8.5,
                     color='#666666', style='italic')

        for grp_val in [0, 1]:
            kmf = kmf_dict[grp_val]
            y_row = n_groups - grp_val - 0.8
            for t_pt in RISK_TIMES:
                if t_pt > xlim:
                    break
                idx = kmf.event_table.index[kmf.event_table.index <= t_pt]
                n_at_risk = int(kmf.event_table.loc[idx[-1], 'at_risk']) if len(idx) else 0
                ax_risk.text(t_pt, y_row, str(n_at_risk),
                             ha='center', va='center', fontsize=8.5,
                             color=COLORS_KM[grp_val], fontweight='bold')

        ax_risk.text(xlim / 2, -0.65, 'Time (months)',
                     ha='center', va='center', fontsize=11, color='#333333')

    plt.savefig(os.path.join(CMP_FIG_DIR, f'ow_subgroup_km_curves_{CMP_TAG}.pdf'),
                bbox_inches='tight', pad_inches=0.1)
    plt.savefig(os.path.join(CMP_FIG_DIR, f'ow_subgroup_km_curves_{CMP_TAG}.png'),
                dpi=300, bbox_inches='tight', pad_inches=0.1)
    plt.close(fig)
    print(f'   KM 曲线已保存: {os.path.join(CMP_FIG_DIR, f"ow_subgroup_km_curves_{CMP_TAG}")}.pdf/.png')

    print('\n\u2705 全部完成！')
