#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
高危亚组分析 — 两组对比
================================================
对比1: HAIC+I_concurrent vs HAIC_then_I
对比2: HAIC+I+T_concurrent vs HAIC_then_I+T

目的: 证明在高危亚组中，延迟免疫(+靶向)治疗与同步免疫(+靶向)治疗的OS无临床有意义差异
方法: 交互作用检验 + 多变量Cox + RMST + 等价性检验(TOST) + 森林图 + KM曲线

输出:
  results/subgroup_analysis/  — 汇总表格CSV (按对比组加后缀)
  figures/subgroup/           — 森林图 + KM曲线 (PDF+PNG, 按对比组加后缀)
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import os
import warnings
warnings.filterwarnings('ignore')

from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test

BASE_DIR = "/Users/yqj/Nutstore Files/我的坚果云/Liver_tumor_big_data/FIRST_LINE_HAIC_2025-12-30/HAIC_NO_TACE_4_TIDY/update_group_7"
DATA_DIR = os.path.join(BASE_DIR, 'data')
RES_DIR  = os.path.join(BASE_DIR, 'results', 'subgroup_analysis')
FIG_DIR  = os.path.join(BASE_DIR, 'figures', 'subgroup')
os.makedirs(RES_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)

# ── rcParams ─────────────────────────────────────────────────────────
plt.rcParams.update({
    'font.family':        'sans-serif',
    'font.sans-serif':    ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size':          10,
    'axes.labelsize':     11,
    'axes.titlesize':     12,
    'xtick.labelsize':     9,
    'ytick.labelsize':     9,
    'legend.fontsize':     9,
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

COLOR_CONC  = '#4DBBD5'
COLOR_THEN  = '#E64B35'

EQUIV_LO = 0.60
EQUIV_HI = 1.67

# ════════════════════════════════════════════════════════════════════
# 0. 数据准备
# ════════════════════════════════════════════════════════════════════

# --- 对比组配置 ---
COMPARISONS = [
    {
        'tag': 'I',
        'treat_group': 'HAIC+I_concurrent',
        'control_group': 'HAIC_then_I',
        'treat_label': 'HAIC+I (concurrent)',
        'control_label': 'HAIC→I (delayed)',
        'title_short': 'HAIC+Immuno Concurrent vs Delayed',
    },
    {
        'tag': 'IT',
        'treat_group': 'HAIC+I+T_concurrent',
        'control_group': 'HAIC_then_I+T',
        'treat_label': 'HAIC+I+T (concurrent)',
        'control_label': 'HAIC→I+T (delayed)',
        'title_short': 'HAIC+Immuno+TKI Concurrent vs Delayed',
    },
]

df_all = pd.read_csv(os.path.join(DATA_DIR, 'analysis_ready.csv'))
df_all = df_all[df_all['os_months'] >= 0].copy()
df_all['event'] = df_all['death_status'].map(
    {'Yes': 1, 'No': 0, '1': 1, '0': 0, 1: 1, 0: 0}
).fillna(0).astype(int)

def prepare_two_group(df_all, treat_group, control_group):
    """筛选并派生两组对比所需的变量"""
    two = df_all[df_all['main_group'].isin([treat_group, control_group])].copy()
    two['treat'] = (two['main_group'] == treat_group).astype(int)

    two['sex_male']         = (two['sex'] == 'Male').astype(int)
    two['pvtt_vp34']        = (two['pvtt_classification'] == 'Vp3/4').astype(int)
    two['pvtt_vp12']        = (two['pvtt_classification'] == 'Vp1/2').astype(int)
    two['pvtt_present']     = (two['pvtt_classification'] != 'Absent').astype(int)
    two['hvtt_yes']         = (two['hvtt'] == 'Yes').astype(int)
    two['ivc_ra_yes']       = (two['ivc_or_ra_thrombus'] == 'Yes').astype(int)
    two['dist_meta_yes']    = (two['distant_metastasis'] == 'Yes').astype(int)
    two['lymph_meta_yes']   = (two['lymph_node_metastasis'] == 'Yes').astype(int)
    two['ascites_yes']      = (two['ascites'] != 'Absent').astype(int)
    two['varices_yes']      = (two['varices'] == 'Yes').astype(int)
    two['tumor_gt10']       = (two['tumor_max_diameter_cm'] > 10).astype(int)
    two['tumor_multi']      = (two['tumor_count_category'] == '>3').astype(int)
    two['afp_high_bin']     = (two['afp_high'] == 'Yes').astype(int)
    two['pivka_high_bin']   = (two['pivka_high'] == 'Yes').astype(int)
    two['albi_grade_num']   = two['albi_grade'].astype(int)

    for col in ['age', 'albi_score', 'nlr', 'tumor_max_diameter_cm', 'alb', 'tbil', 'plt']:
        two[f'{col}_std'] = (two[col] - two[col].mean()) / two[col].std()

    two['high_risk_composite'] = (
        two['ivc_ra_yes'] | two['tumor_multi'] | two['pvtt_vp34'] |
        two['dist_meta_yes'] | two['tumor_gt10']
    ).astype(int)
    return two


# ════════════════════════════════════════════════════════════════════
# 1. 定义亚组 & 协变量（两组对比共用）
# ════════════════════════════════════════════════════════════════════
SUBGROUPS = [
    {
        'name': 'Tumor count >3',
        'col': 'tumor_multi',
        'val': 1,
        'exclude_from_covars': ['tumor_multi'],
    },
    {
        'name': 'Tumor diameter >10 cm',
        'col': 'tumor_gt10',
        'val': 1,
        'exclude_from_covars': ['tumor_gt10', 'tumor_max_diameter_cm_std'],
    },
    {
        'name': 'PVTT Vp3/4',
        'col': 'pvtt_vp34',
        'val': 1,
        'exclude_from_covars': ['pvtt_vp34', 'pvtt_vp12', 'pvtt_present'],
    },
    {
        'name': 'Extrahepatic metastasis',
        'col': 'dist_meta_yes',
        'val': 1,
        'exclude_from_covars': ['dist_meta_yes'],
    },
]

BASE_COVARS = [
    'age_std', 'sex_male', 'albi_score_std', 'afp_high_bin', 'pivka_high_bin',
    'tumor_max_diameter_cm_std', 'tumor_multi',
    'pvtt_vp34', 'pvtt_vp12', 'dist_meta_yes',
    'ascites_yes', 'nlr_std',
]


def select_covars(base, exclude, n_events):
    """按 EPV ≥ 10 规则选择调整变量"""
    covars = [c for c in base if c not in exclude]
    max_vars = max(1, n_events // 10 - 1)
    return covars[:max_vars]


# ════════════════════════════════════════════════════════════════════
# 2. RMST 计算函数
# ════════════════════════════════════════════════════════════════════
def compute_rmst(time, event, tau):
    """计算 RMST 及其标准误（Greenwood 方法）"""
    kmf = KaplanMeierFitter()
    kmf.fit(time, event)
    sf = kmf.survival_function_.copy()
    sf = sf[sf.index <= tau]
    if sf.index.max() < tau:
        last_s = sf.iloc[-1, 0]
        sf.loc[tau] = last_s
    sf = sf.sort_index()

    t_vals = sf.index.values
    s_vals = sf.iloc[:, 0].values
    dt = np.diff(t_vals)
    s_mid = s_vals[:-1]
    rmst = np.sum(s_mid * dt)

    et = kmf.event_table
    et = et[et.index <= tau]
    d = et['observed'].values
    n = et['at_risk'].values
    var_terms = np.zeros(len(et))
    for i in range(len(et)):
        if n[i] > 0 and n[i] > d[i]:
            var_terms[i] = d[i] / (n[i] * (n[i] - d[i]))

    t_et = et.index.values
    area_after = np.zeros(len(t_et))
    for i in range(len(t_et)):
        mask = sf.index >= t_et[i]
        t_sub = sf.index[mask].values
        s_sub = sf.iloc[:, 0][mask].values
        if len(t_sub) > 1:
            dt_sub = np.diff(t_sub)
            area_after[i] = np.sum(s_sub[:-1] * dt_sub)

    var_rmst = np.sum(area_after**2 * var_terms)
    se_rmst = np.sqrt(max(var_rmst, 0))
    return rmst, se_rmst


def rmst_diff(t1, e1, t2, e2, tau):
    """两组 RMST 差值及 95% CI"""
    r1, se1 = compute_rmst(t1, e1, tau)
    r2, se2 = compute_rmst(t2, e2, tau)
    diff = r1 - r2
    se_diff = np.sqrt(se1**2 + se2**2)
    ci_lo = diff - 1.96 * se_diff
    ci_hi = diff + 1.96 * se_diff
    z = diff / se_diff if se_diff > 0 else 0
    p = 2 * (1 - __import__('scipy').stats.norm.cdf(abs(z)))
    return {'rmst1': r1, 'rmst2': r2, 'diff': diff,
            'se': se_diff, 'ci_lo': ci_lo, 'ci_hi': ci_hi, 'p': p}


# ════════════════════════════════════════════════════════════════════
# 遍历每组对比
# ════════════════════════════════════════════════════════════════════
for COMP in COMPARISONS:
    CMP_TAG       = COMP['tag']
    TREAT_GROUP   = COMP['treat_group']
    CONTROL_GROUP = COMP['control_group']
    TREAT_LABEL   = COMP['treat_label']
    CTRL_LABEL    = COMP['control_label']
    TITLE_SHORT   = COMP['title_short']

    # 输出目录加后缀
    CMP_RES_DIR = os.path.join(RES_DIR, CMP_TAG)
    CMP_FIG_DIR = os.path.join(FIG_DIR, CMP_TAG)
    os.makedirs(CMP_RES_DIR, exist_ok=True)
    os.makedirs(CMP_FIG_DIR, exist_ok=True)

    print('\n' + '#' * 70)
    print(f'#  对比组: {TREAT_GROUP} vs {CONTROL_GROUP}  (tag={CMP_TAG})')
    print('#' * 70)

    two = prepare_two_group(df_all, TREAT_GROUP, CONTROL_GROUP)

    print(f'\n总样本: {len(two)} (Concurrent={two["treat"].sum()}, Delayed={(1-two["treat"]).sum()})')
    print(f'事件数: {two["event"].sum()} ({two["event"].mean()*100:.1f}%)')

    # ════════════════════════════════════════════════════════════════
    # 3. 交互作用检验（全人群）
    # ════════════════════════════════════════════════════════════════
    print('\n' + '=' * 70)
    print('  Step 1: 交互作用检验（全人群 n=%d）' % len(two))
    print('=' * 70)

    interaction_results = []
    for sg in SUBGROUPS:
        sg_col = sg['col']
        sg_name = sg['name']

        cph_data = two[['os_months', 'event', 'treat', sg_col,
                         'age_std', 'albi_score_std', 'afp_high_bin',
                         'nlr_std']].copy()
        cph_data['treat_x_sg'] = cph_data['treat'] * cph_data[sg_col]

        cph = CoxPHFitter()
        cph.fit(cph_data, duration_col='os_months', event_col='event')

        inter_hr = np.exp(cph.params_['treat_x_sg'])
        inter_p  = cph.summary.loc['treat_x_sg', 'p']
        inter_ci_lo = np.exp(cph.confidence_intervals_.loc['treat_x_sg'].iloc[0])
        inter_ci_hi = np.exp(cph.confidence_intervals_.loc['treat_x_sg'].iloc[1])

        interaction_results.append({
            'Subgroup': sg_name,
            'Interaction_HR': round(inter_hr, 3),
            'Interaction_CI_lo': round(inter_ci_lo, 3),
            'Interaction_CI_hi': round(inter_ci_hi, 3),
            'P_interaction': round(inter_p, 4),
        })
        sig = '***' if inter_p < 0.01 else ('**' if inter_p < 0.05 else ('*' if inter_p < 0.10 else ''))
        print(f'  {sg_name:30s} | P_interaction = {inter_p:.4f} {sig}')

    # 复合高危
    cph_data = two[['os_months', 'event', 'treat', 'high_risk_composite',
                     'age_std', 'albi_score_std', 'afp_high_bin', 'nlr_std']].copy()
    cph_data['treat_x_hr'] = cph_data['treat'] * cph_data['high_risk_composite']
    cph = CoxPHFitter()
    cph.fit(cph_data, duration_col='os_months', event_col='event')
    inter_p_comp = cph.summary.loc['treat_x_hr', 'p']
    print(f'  {"Composite high-risk":30s} | P_interaction = {inter_p_comp:.4f}')

    interaction_df = pd.DataFrame(interaction_results)
    interaction_df.to_csv(os.path.join(CMP_RES_DIR, f'interaction_tests_{CMP_TAG}.csv'), index=False)

    # ════════════════════════════════════════════════════════════════════
    # 4. 亚组内生存分析
    # ════════════════════════════════════════════════════════════════════
    print('\n' + '=' * 70)
    print('  Step 2: 亚组内生存分析')
    print('=' * 70)

    all_results = []

    def analyze_subgroup(data, sg_name, sg_exclude_covars, is_complement=False):
        """对单个亚组执行完整分析"""
        label = sg_name if not is_complement else f'Non-{sg_name}'
        grp1 = data[data['treat'] == 1]
        grp0 = data[data['treat'] == 0]
        n1, n0 = len(grp1), len(grp0)
        ev1, ev0 = grp1['event'].sum(), grp0['event'].sum()
        n_events = ev1 + ev0

        if n1 < 5 or n0 < 5 or n_events < 10:
            print(f'    ⚠ {label}: 样本量不足 (n={n1}+{n0}, events={n_events})，跳过')
            return None

        print(f'\n  --- {label} (n={n1}+{n0}, events={n_events}) ---')

        # KM + Log-rank
        kmf1 = KaplanMeierFitter().fit(grp1['os_months'], grp1['event'], label='Concurrent')
        kmf0 = KaplanMeierFitter().fit(grp0['os_months'], grp0['event'], label='Then_I')
        lr = logrank_test(grp1['os_months'], grp0['os_months'], grp1['event'], grp0['event'])
        med1 = kmf1.median_survival_time_
        med0 = kmf0.median_survival_time_

        # 未调整 Cox
        cox_data = data[['os_months', 'event', 'treat']].copy()
        cph_unadj = CoxPHFitter()
        cph_unadj.fit(cox_data, duration_col='os_months', event_col='event')
        hr_unadj = np.exp(cph_unadj.params_['treat'])
        ci_lo_unadj = np.exp(cph_unadj.confidence_intervals_.loc['treat'].iloc[0])
        ci_hi_unadj = np.exp(cph_unadj.confidence_intervals_.loc['treat'].iloc[1])
        p_unadj = cph_unadj.summary.loc['treat', 'p']

        # 90% CI for equivalence
        cph_90 = CoxPHFitter(alpha=0.10)
        cph_90.fit(cox_data, duration_col='os_months', event_col='event')
        ci90_lo = np.exp(cph_90.confidence_intervals_.loc['treat'].iloc[0])
        ci90_hi = np.exp(cph_90.confidence_intervals_.loc['treat'].iloc[1])

        # 调整 Cox
        covars = select_covars(BASE_COVARS, sg_exclude_covars, n_events)
        hr_adj = ci_lo_adj = ci_hi_adj = p_adj = np.nan
        ci90_lo_adj = ci90_hi_adj = np.nan
        if n_events >= 20 and len(covars) >= 1:
            adj_cols = ['os_months', 'event', 'treat'] + covars
            adj_data = data[adj_cols].dropna()
            if len(adj_data) > 20:
                try:
                    cph_adj = CoxPHFitter()
                    cph_adj.fit(adj_data, duration_col='os_months', event_col='event')
                    hr_adj = np.exp(cph_adj.params_['treat'])
                    ci_lo_adj = np.exp(cph_adj.confidence_intervals_.loc['treat'].iloc[0])
                    ci_hi_adj = np.exp(cph_adj.confidence_intervals_.loc['treat'].iloc[1])
                    p_adj = cph_adj.summary.loc['treat', 'p']

                    cph_adj90 = CoxPHFitter(alpha=0.10)
                    cph_adj90.fit(adj_data, duration_col='os_months', event_col='event')
                    ci90_lo_adj = np.exp(cph_adj90.confidence_intervals_.loc['treat'].iloc[0])
                    ci90_hi_adj = np.exp(cph_adj90.confidence_intervals_.loc['treat'].iloc[1])
                    print(f'    调整Cox ({len(covars)} vars): HR={hr_adj:.2f} ({ci_lo_adj:.2f}-{ci_hi_adj:.2f}) P={p_adj:.4f}')
                except Exception as e:
                    print(f'    ⚠ 调整Cox失败: {e}')

        # RMST (tau=24)
        rmst_24 = rmst_diff(grp1['os_months'].values, grp1['event'].values,
                             grp0['os_months'].values, grp0['event'].values, tau=24)
        # RMST (tau=12, 敏感性)
        rmst_12 = rmst_diff(grp1['os_months'].values, grp1['event'].values,
                             grp0['os_months'].values, grp0['event'].values, tau=12)

        # 等价性判断 (TOST)
        equiv_unadj = (ci90_lo >= EQUIV_LO and ci90_hi <= EQUIV_HI)
        equiv_adj = False
        if not np.isnan(ci90_lo_adj):
            equiv_adj = (ci90_lo_adj >= EQUIV_LO and ci90_hi_adj <= EQUIV_HI)

        # E-value
        e_val = hr_unadj + np.sqrt(hr_unadj * (hr_unadj - 1)) if hr_unadj >= 1 else \
                1/hr_unadj + np.sqrt(1/hr_unadj * (1/hr_unadj - 1))

        print(f'    未调整: HR={hr_unadj:.2f} ({ci_lo_unadj:.2f}-{ci_hi_unadj:.2f}) P={p_unadj:.4f}')
        print(f'    90% CI: ({ci90_lo:.2f}-{ci90_hi:.2f}) → 等价性: {"✓" if equiv_unadj else "✗"}')
        print(f'    RMST@24m: Δ={rmst_24["diff"]:.2f} ({rmst_24["ci_lo"]:.2f} to {rmst_24["ci_hi"]:.2f}) P={rmst_24["p"]:.4f}')
        print(f'    RMST@12m: Δ={rmst_12["diff"]:.2f} ({rmst_12["ci_lo"]:.2f} to {rmst_12["ci_hi"]:.2f}) P={rmst_12["p"]:.4f}')
        print(f'    Median OS: {med1:.1f} vs {med0:.1f} mo')

        result = {
            'Subgroup': label,
            'N_concurrent': n1, 'N_then_I': n0,
            'Events_concurrent': ev1, 'Events_then_I': ev0, 'Events_total': n_events,
            'Median_OS_concurrent': round(med1, 1), 'Median_OS_then_I': round(med0, 1),
            'HR_unadjusted': round(hr_unadj, 3),
            'CI95_lo_unadj': round(ci_lo_unadj, 3), 'CI95_hi_unadj': round(ci_hi_unadj, 3),
            'P_unadjusted': round(p_unadj, 4),
            'CI90_lo_unadj': round(ci90_lo, 3), 'CI90_hi_unadj': round(ci90_hi, 3),
            'HR_adjusted': round(hr_adj, 3) if not np.isnan(hr_adj) else np.nan,
            'CI95_lo_adj': round(ci_lo_adj, 3) if not np.isnan(ci_lo_adj) else np.nan,
            'CI95_hi_adj': round(ci_hi_adj, 3) if not np.isnan(ci_hi_adj) else np.nan,
            'P_adjusted': round(p_adj, 4) if not np.isnan(p_adj) else np.nan,
            'CI90_lo_adj': round(ci90_lo_adj, 3) if not np.isnan(ci90_lo_adj) else np.nan,
            'CI90_hi_adj': round(ci90_hi_adj, 3) if not np.isnan(ci90_hi_adj) else np.nan,
            'N_adjust_vars': len(covars) if n_events >= 20 else 0,
            'Equivalence_unadj': equiv_unadj,
            'Equivalence_adj': equiv_adj,
            'RMST24_diff': round(rmst_24['diff'], 2),
            'RMST24_CI_lo': round(rmst_24['ci_lo'], 2),
            'RMST24_CI_hi': round(rmst_24['ci_hi'], 2),
            'RMST24_P': round(rmst_24['p'], 4),
            'RMST12_diff': round(rmst_12['diff'], 2),
            'RMST12_CI_lo': round(rmst_12['ci_lo'], 2),
            'RMST12_CI_hi': round(rmst_12['ci_hi'], 2),
            'RMST12_P': round(rmst_12['p'], 4),
            'E_value': round(e_val, 2),
            'Logrank_P': round(lr.p_value, 4),
            'kmf_conc': kmf1, 'kmf_then': kmf0,
        }
        return result

    # 整体（未匹配）
    print('\n  === Overall (unmatched) ===')
    res_overall = analyze_subgroup(two, 'Overall (unmatched)', [])
    if res_overall:
        all_results.append(res_overall)

    # 复合高危
    print('\n  === Composite high-risk ===')
    res_comp = analyze_subgroup(two[two['high_risk_composite'] == 1],
                                'Composite high-risk (≥1 factor)', [])
    if res_comp:
        all_results.append(res_comp)

    # 各亚组 + 互补亚组
    for sg in SUBGROUPS:
        sg_name = sg['name']
        sg_col  = sg['col']
        sg_val  = sg['val']
        excl    = sg['exclude_from_covars']

        sub_hi = two[two[sg_col] == sg_val]
        sub_lo = two[two[sg_col] != sg_val]

        res_hi = analyze_subgroup(sub_hi, sg_name, excl)
        res_lo = analyze_subgroup(sub_lo, sg_name, excl, is_complement=True)

        if res_hi:
            all_results.append(res_hi)
        if res_lo:
            all_results.append(res_lo)

    # 保存结果
    cols_save = [c for c in all_results[0].keys() if not c.startswith('kmf_')]
    results_df = pd.DataFrame([{k: v for k, v in r.items() if k in cols_save} for r in all_results])
    results_df.to_csv(os.path.join(CMP_RES_DIR, f'subgroup_survival_results_{CMP_TAG}.csv'), index=False)
    _save_path = os.path.join(CMP_RES_DIR, f'subgroup_survival_results_{CMP_TAG}.csv')
    print(f'\n结果已保存: {_save_path}')

    # ════════════════════════════════════════════════════════════════════
    # 5. 森林图
    # ════════════════════════════════════════════════════════════════════
    print('\n' + '=' * 70)
    print('  Step 3: 绘制森林图')
    print('=' * 70)

    forest_rows = []

    # Overall
    r = all_results[0]
    forest_rows.append({
        'label': 'Overall (unmatched)',
        'n_str': f'{r["N_concurrent"]}+{r["N_then_I"]}',
        'events': r['Events_total'],
        'hr': r['HR_unadjusted'], 'ci_lo': r['CI95_lo_unadj'], 'ci_hi': r['CI95_hi_unadj'],
        'p': r['P_unadjusted'],
        'ci90_lo': r['CI90_lo_unadj'], 'ci90_hi': r['CI90_hi_unadj'],
        'equiv': r['Equivalence_unadj'],
        'rmst_diff': r['RMST24_diff'], 'rmst_ci': f"({r['RMST24_CI_lo']:.1f}, {r['RMST24_CI_hi']:.1f})",
        'p_inter': np.nan, 'is_header': False, 'is_overall': True, 'indent': 0,
    })

    # PSM matched (from step3 results) — only for I comparison
    psm_path = os.path.join(BASE_DIR, 'results', 'psm_balance_tables_complete', 'survival_analysis_final.csv')
    if CMP_TAG == 'I' and os.path.exists(psm_path):
        psm_df = pd.read_csv(psm_path)
        psm_row = psm_df[(psm_df['Group1'] == TREAT_GROUP) & (psm_df['Group2'] == CONTROL_GROUP)]
        if len(psm_row):
            pr = psm_row.iloc[0]
            forest_rows.append({
                'label': 'Overall (PSM matched)',
                'n_str': f'{int(pr["N1_after"])}+{int(pr["N2_after"])}',
                'events': '',
                'hr': pr['HR'], 'ci_lo': pr['CI_lower'], 'ci_hi': pr['CI_upper'],
                'p': pr['P_value'],
                'ci90_lo': np.nan, 'ci90_hi': np.nan,
                'equiv': '',
                'rmst_diff': '', 'rmst_ci': '',
                'p_inter': np.nan, 'is_header': False, 'is_overall': True, 'indent': 0,
            })

    # Composite high-risk
    r = all_results[1]
    forest_rows.append({
        'label': 'Composite high-risk',
        'n_str': f'{r["N_concurrent"]}+{r["N_then_I"]}',
        'events': r['Events_total'],
        'hr': r['HR_unadjusted'], 'ci_lo': r['CI95_lo_unadj'], 'ci_hi': r['CI95_hi_unadj'],
        'p': r['P_unadjusted'],
        'ci90_lo': r['CI90_lo_unadj'], 'ci90_hi': r['CI90_hi_unadj'],
        'equiv': r['Equivalence_unadj'],
        'rmst_diff': r['RMST24_diff'], 'rmst_ci': f"({r['RMST24_CI_lo']:.1f}, {r['RMST24_CI_hi']:.1f})",
        'p_inter': round(inter_p_comp, 4), 'is_header': False, 'is_overall': False, 'indent': 0,
    })

    # 各亚组
    idx = 2
    for sg_i, sg in enumerate(SUBGROUPS):
        p_inter = interaction_results[sg_i]['P_interaction']

        for j in range(2):
            if idx >= len(all_results):
                break
            r = all_results[idx]
            is_hi = not r['Subgroup'].startswith('Non-')
            forest_rows.append({
                'label': r['Subgroup'],
                'n_str': f'{r["N_concurrent"]}+{r["N_then_I"]}',
                'events': r['Events_total'],
                'hr': r['HR_unadjusted'], 'ci_lo': r['CI95_lo_unadj'], 'ci_hi': r['CI95_hi_unadj'],
                'p': r['P_unadjusted'],
                'ci90_lo': r['CI90_lo_unadj'], 'ci90_hi': r['CI90_hi_unadj'],
                'equiv': r['Equivalence_unadj'],
                'rmst_diff': r['RMST24_diff'],
                'rmst_ci': f"({r['RMST24_CI_lo']:.1f}, {r['RMST24_CI_hi']:.1f})",
                'p_inter': p_inter if is_hi else np.nan,
                'is_header': False, 'is_overall': False,
                'indent': 1,
            })
            idx += 1

    n_rows = len(forest_rows)
    fig_h = max(5.5, n_rows * 0.48 + 2.5)
    fig, ax = plt.subplots(figsize=(12.0, fig_h))

    y_positions = list(range(n_rows - 1, -1, -1))
    X_MIN, X_MAX = 0.20, 5.0
    ax.set_xlim(X_MIN, X_MAX)
    ax.set_xscale('log')
    ax.set_ylim(-1.0, n_rows + 0.5)
    ax.set_yticks([])
    ax.spines['left'].set_visible(False)

    ax.axvspan(EQUIV_LO, EQUIV_HI, color='#E8F4E8', alpha=0.5, zorder=0,
               label=f'Equivalence zone [{EQUIV_LO}–{EQUIV_HI}]')
    ax.axvline(x=1.0, color='#444444', linestyle='--', linewidth=0.8, zorder=1)

    y_top = n_rows + 0.2
    ax.text(0.005, y_top, 'Subgroup', ha='left', va='center',
            fontsize=8.5, fontweight='bold', color='#222222',
            transform=ax.get_yaxis_transform())
    ax.text(1.0, y_top, 'HR (95% CI)', ha='center', va='center',
            fontsize=8.5, fontweight='bold', color='#222222')
    col_x = X_MAX * 1.03
    ax.text(col_x, y_top, 'HR (95% CI)         P        P_int',
            ha='left', va='center', fontsize=7.5, fontweight='bold', color='#222222')

    ax.axhline(y=n_rows - 0.3, color='#333333', linewidth=0.8, xmin=0, xmax=1, clip_on=False)

    prev_was_overall = False
    for i, (row, ypos) in enumerate(zip(forest_rows, y_positions)):
        hr    = float(row['hr'])
        ci_lo = float(row['ci_lo'])
        ci_hi = float(row['ci_hi'])
        p_val = float(row['p']) if row['p'] != '' else np.nan

        if row['is_overall'] and not prev_was_overall:
            pass
        if not row['is_overall'] and prev_was_overall:
            ax.axhline(y=ypos + 0.55, color='#BBBBBB', linewidth=0.6, xmin=0, xmax=1, clip_on=False)
        prev_was_overall = row['is_overall']

        if i % 2 == 0 and not row['is_overall']:
            ax.axhspan(ypos - 0.4, ypos + 0.4, color='#F7F7F7', zorder=0, linewidth=0)

        indent = '  ' if row.get('indent', 0) else ''
        lbl_weight = 'bold' if row['is_overall'] or row.get('indent', 0) == 0 else 'normal'
        lbl_size = 8 if row.get('indent', 0) else 8.5
        ax.text(0.005, ypos,
                f'{indent}{row["label"]}  (n={row["n_str"]})',
                ha='left', va='center', fontsize=lbl_size,
                fontweight=lbl_weight, color='#222222',
                transform=ax.get_yaxis_transform())

        dot_color = '#3C5488' if row['is_overall'] else '#E64B35'
        ci_lo_plot = max(ci_lo, X_MIN * 1.02)
        ci_hi_plot = min(ci_hi, X_MAX * 0.98)

        ax.plot([ci_lo_plot, ci_hi_plot], [ypos, ypos],
                color=dot_color, linewidth=1.6, solid_capstyle='round', zorder=3)
        for xc in [ci_lo_plot, ci_hi_plot]:
            ax.plot([xc, xc], [ypos - 0.12, ypos + 0.12],
                    color=dot_color, linewidth=1.0, zorder=3)
        marker_size = 45 if row['is_overall'] else 30
        ax.scatter([hr], [ypos], marker='D', s=marker_size,
                   color=dot_color, zorder=4, linewidths=0)

        p_str = '<0.001' if (not np.isnan(p_val) and p_val < 0.001) else \
                f'{p_val:.3f}' if not np.isnan(p_val) else ''
        p_inter = row['p_inter']
        pi_str = '' if (isinstance(p_inter, float) and np.isnan(p_inter)) else \
                 f'{p_inter:.3f}' if p_inter >= 0.001 else '<0.001'

        right_txt = f'{hr:.2f} ({ci_lo:.2f}\u2013{ci_hi:.2f})   {p_str:>7s}   {pi_str:>7s}'
        ax.text(col_x, ypos, right_txt, ha='left', va='center',
                fontsize=7, color='#333333', family='monospace')

    ax.axhline(y=-0.6, color='#333333', linewidth=0.8, xmin=0, xmax=1, clip_on=False)

    xticks = [0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0]
    ax.set_xticks(xticks)
    ax.set_xticklabels([str(x) for x in xticks], fontsize=7.5)
    ax.set_xlabel('Hazard Ratio (95% CI)  [log scale]', fontsize=9, labelpad=6)

    ax.text(0.32, -0.03, f'\u2190 Favors {TREAT_LABEL}', ha='center', va='top',
            fontsize=7.5, color='#666666', transform=ax.transAxes)
    ax.text(0.68, -0.03, f'Favors {CTRL_LABEL} \u2192', ha='center', va='top',
            fontsize=7.5, color='#666666', transform=ax.transAxes)

    legend_elements = [
        mpatches.Patch(facecolor='#E8F4E8', edgecolor='#AACCAA', alpha=0.7,
                       label=f'Equivalence zone (HR {EQUIV_LO}\u2013{EQUIV_HI})'),
        mlines.Line2D([0], [0], color='#444444', linestyle='--', linewidth=0.8,
                      label='HR = 1.0 (null)'),
    ]
    ax.legend(handles=legend_elements, loc='lower right',
              bbox_to_anchor=(0.98, -0.07), fontsize=7, frameon=False)

    ax.set_title(f'Subgroup Analysis: {TITLE_SHORT}\n'
                 'Overall Survival — High-Risk Subgroups',
                 fontsize=11, fontweight='bold', pad=10, loc='left')

    plt.tight_layout()
    out_base = os.path.join(CMP_FIG_DIR, f'subgroup_forest_plot_{CMP_TAG}')
    fig.savefig(f'{out_base}.pdf', bbox_inches='tight', pad_inches=0.05)
    fig.savefig(f'{out_base}.png', dpi=300, bbox_inches='tight', pad_inches=0.05)
    plt.close(fig)
    print(f'  森林图已保存: {out_base}.pdf/.png')

    # ════════════════════════════════════════════════════════════════════
    # 6. 亚组 KM 曲线（多面板）
    # ════════════════════════════════════════════════════════════════════
    print('\n' + '=' * 70)
    print('  Step 4: 绘制亚组 KM 曲线')
    print('=' * 70)

    km_results = [r for r in all_results if 'kmf_conc' in r and not r['Subgroup'].startswith('Non-')]

    n_panels = len(km_results)
    ncols = 3
    nrows = (n_panels + ncols - 1) // ncols

    fig = plt.figure(figsize=(5.0 * ncols, 5.5 * nrows))
    outer_gs = fig.add_gridspec(nrows, ncols, hspace=0.45, wspace=0.30)

    for panel_i, r in enumerate(km_results):
        row_i = panel_i // ncols
        col_i = panel_i % ncols

        inner_gs = outer_gs[row_i, col_i].subgridspec(2, 1, height_ratios=[4.2, 0.8], hspace=0.04)
        ax_km   = fig.add_subplot(inner_gs[0])
        ax_risk = fig.add_subplot(inner_gs[1])

        kmf_c = r['kmf_conc']
        kmf_t = r['kmf_then']

        for kmf, color, ls, lbl_key in [
            (kmf_c, COLOR_CONC, '-',  'Concurrent'),
            (kmf_t, COLOR_THEN, '--', 'Then_I'),
        ]:
            n_at_risk = int(kmf.event_table['at_risk'].iloc[0])
            t = kmf.survival_function_.index.values
            s = kmf.survival_function_.iloc[:, 0].values
            ci_lo_km = kmf.confidence_interval_.iloc[:, 0].values
            ci_hi_km = kmf.confidence_interval_.iloc[:, 1].values
            ax_km.step(t, s, where='post', color=color, linewidth=1.8,
                       linestyle=ls, label=f'{lbl_key} (n={n_at_risk})')
            ax_km.fill_between(t, ci_lo_km, ci_hi_km, step='post', alpha=0.08, color=color)

        ax_km.set_xlim(0, 60)
        ax_km.set_ylim(-0.02, 1.05)
        ax_km.set_ylabel('OS Probability', fontsize=9)
        ax_km.set_xticks(RISK_TIMES)
        ax_km.yaxis.set_major_locator(mticker.MultipleLocator(0.2))
        ax_km.legend(loc='upper right', fontsize=7.5, handlelength=2.0)

        panel_letter = chr(65 + panel_i)
        title_str = f'{r["Subgroup"]}'
        ax_km.set_title(title_str, fontsize=10, fontweight='bold', pad=6)
        ax_km.text(-0.12, 1.08, panel_letter, transform=ax_km.transAxes,
                   fontsize=11, fontweight='bold', va='top')

        p_fmt = 'P < 0.001' if r['Logrank_P'] < 0.001 else f'P = {r["Logrank_P"]:.3f}'
        hr_val = r['HR_unadjusted']
        ci_l = r['CI95_lo_unadj']
        ci_h = r['CI95_hi_unadj']
        stats_lines = [
            p_fmt,
            f'HR = {hr_val:.2f} ({ci_l:.2f}\u2013{ci_h:.2f})',
            f'ΔRMST@24m = {r["RMST24_diff"]:.1f} mo',
        ]
        equiv_str = '✓ Equivalent' if r['Equivalence_unadj'] else '? Inconclusive'
        stats_lines.append(equiv_str)
        ax_km.text(0.97, 0.05, '\n'.join(stats_lines), transform=ax_km.transAxes,
                   fontsize=7, va='bottom', ha='right', bbox=STATS_BOX)

        # Risk table
        ax_risk.set_xlim(0, 60.5)
        ax_risk.set_ylim(-0.8, 2.2)
        ax_risk.axis('off')
        ax_risk.text(30, 2.0, 'No. at risk', ha='center', va='top',
                     fontsize=7.5, color='#666666', style='italic')

        for ki, (kmf, color) in enumerate([(kmf_c, COLOR_CONC), (kmf_t, COLOR_THEN)]):
            y_rt = 1.0 - ki * 0.9
            for t_pt in RISK_TIMES:
                idx = kmf.event_table.index[kmf.event_table.index <= t_pt]
                n_val = int(kmf.event_table.loc[idx[-1], 'at_risk']) if len(idx) else 0
                ax_risk.text(t_pt, y_rt, str(n_val), ha='center', va='center',
                             fontsize=7.5, color=color, fontweight='bold')

        ax_risk.text(30, -0.65, 'Time (months)', ha='center', va='center',
                     fontsize=9, color='#333333')

    # 隐藏多余面板
    for extra_i in range(n_panels, nrows * ncols):
        row_i = extra_i // ncols
        col_i = extra_i % ncols
        inner_gs = outer_gs[row_i, col_i].subgridspec(2, 1, height_ratios=[4.2, 0.8], hspace=0.04)
        ax1 = fig.add_subplot(inner_gs[0])
        ax2 = fig.add_subplot(inner_gs[1])
        ax1.axis('off')
        ax2.axis('off')

    fig.suptitle(f'Subgroup KM Curves: {TITLE_SHORT}',
                 fontsize=13, fontweight='bold', y=1.01)

    out_km = os.path.join(CMP_FIG_DIR, f'subgroup_km_curves_{CMP_TAG}')
    fig.savefig(f'{out_km}.pdf', bbox_inches='tight', pad_inches=0.05)
    fig.savefig(f'{out_km}.png', dpi=300, bbox_inches='tight', pad_inches=0.05)
    plt.close(fig)
    print(f'  KM 曲线已保存: {out_km}.pdf/.png')

    # ════════════════════════════════════════════════════════════════════
    # 7. 汇总打印
    # ════════════════════════════════════════════════════════════════════
    print('\n' + '=' * 70)
    print('  分析完成！结果汇总')
    print('=' * 70)

    print(f'\n等价性界值: HR [{EQUIV_LO}, {EQUIV_HI}]')
    print(f'RMST 截断时间: τ = 24 月 (敏感性: τ = 12 月)')
    print(f'\n{"Subgroup":<35s} {"HR (95%CI)":>22s} {"P":>8s} {"90%CI":>18s} {"Equiv":>6s} {"ΔRMST@24":>10s}')
    print('-' * 105)
    for r in all_results:
        sg = r['Subgroup']
        hr = r['HR_unadjusted']
        cl = r['CI95_lo_unadj']
        ch = r['CI95_hi_unadj']
        p  = r['P_unadjusted']
        c90l = r['CI90_lo_unadj']
        c90h = r['CI90_hi_unadj']
        eq = '✓' if r['Equivalence_unadj'] else '✗'
        rd = r['RMST24_diff']
        p_str = '<0.001' if p < 0.001 else f'{p:.3f}'
        print(f'{sg:<35s} {hr:.2f} ({cl:.2f}-{ch:.2f})  {p_str:>7s}  ({c90l:.2f}-{c90h:.2f})  {eq:>5s}  {rd:>+8.2f} mo')

    print(f'\n交互作用检验:')
    for ir in interaction_results:
        print(f'  {ir["Subgroup"]:<30s} P_interaction = {ir["P_interaction"]:.4f}')

    print('\n输出文件:')
    print(f'  {CMP_RES_DIR}/subgroup_survival_results_{CMP_TAG}.csv')
    print(f'  {CMP_RES_DIR}/interaction_tests_{CMP_TAG}.csv')
    print(f'  {CMP_FIG_DIR}/subgroup_forest_plot_{CMP_TAG}.pdf/.png')
    print(f'  {CMP_FIG_DIR}/subgroup_km_curves_{CMP_TAG}.pdf/.png')
