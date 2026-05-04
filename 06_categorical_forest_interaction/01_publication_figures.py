#!/usr/bin/env python3
"""
============================================================================
PSM matched_ids_02：HAIC alone vs HAIC then I — 发表级森林图（Fig4 + Landmark 补充）
============================================================================
队列：仅 matched_ids_02 中患者；数据源自 HAIC_NO_TACE_4_TIDY_baseline（=pre-HAIC-1）与
  longitudinal 的 pre_haic_3；分组 = CSV 的 treatment（0/1）。非加权 Cox / KM / LRT。

本目录数据根 ../../；RCS 交叉点 CSV：../nonlinear_analysis/method1_rcs/

产出：Fig4_forest_v4；OUT_PREFIX 的 SuppFig_lm_forest；Table_lm_forest.csv。
============================================================================
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
from lifelines import KaplanMeierFitter, CoxPHFitter
from scipy import stats
import warnings
import os
from collections import defaultdict
warnings.filterwarnings('ignore')

# ============================================================================
# 全局配置 — NPG 配色（Nature Publishing Group 标准）
# ============================================================================
# NPG 主色
CB_BLUE   = '#4DBBD5'   # NPG 青色（对照/基线）
CB_ORANGE = '#E64B35'   # NPG 红橙（治疗/实验）
CB_GREEN  = '#00A087'   # NPG 蓝绿（改善/阳性）
CB_PINK   = '#F39B7F'   # NPG 浅珊瑚
CB_LBLUE  = '#8491B4'   # NPG 灰蓝
CB_GRAY   = '#999999'   # 中性灰
CB_RED    = '#DC0000'   # NPG 深红（警示/高风险）
CB_YELLOW = '#3C5488'   # NPG 深蓝（第4组）

# Fig4：按 treatment×分层因子 LRT 的 P_interaction 给 HR/ΔRMST 图形与 P 列上色
P_INTERACT_BRIGHT_RED = '#FF1744'   # 鲜红：P_interaction < 0.05
P_INTERACT_DARK_RED   = '#7F1D1D'   # 暗红：0.05 ≤ P_interaction < 0.10

# NPG 完整调色板
NPG_PALETTE = ['#E64B35', '#4DBBD5', '#00A087', '#3C5488',
               '#F39B7F', '#8491B4', '#91D1C2', '#DC0000', '#7E6148', '#B09C85']

# 2组对比标准配色
COLOR_IMMUNO = '#E64B35'   # 免疫治疗组（红橙）
COLOR_HAIC   = '#4DBBD5'   # HAIC单药组（青色）

# 统计框样式（KM图右下角）
STATS_BOX = dict(boxstyle='round,pad=0.4', facecolor='white',
                 edgecolor='#CCCCCC', alpha=0.95, linewidth=0.5)
# 风险表时间点
RISK_TIMES = [0, 12, 24, 36, 48, 60]

import matplotlib.ticker as mticker

plt.rcParams.update({
    # ── 字体 ──
    'font.family':        'sans-serif',
    'font.sans-serif':    ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size':           8,
    'axes.labelsize':      9,
    'axes.titlesize':     10,
    'xtick.labelsize':     7,
    'ytick.labelsize':     7,
    'legend.fontsize':     7,
    # ── 线条与边框 ──
    'axes.linewidth':      0.6,
    'axes.spines.top':     False,
    'axes.spines.right':   False,
    'xtick.major.width':   0.6,
    'ytick.major.width':   0.6,
    'xtick.major.size':    3.0,
    'ytick.major.size':    3.0,
    'xtick.minor.size':    1.5,
    'ytick.minor.size':    1.5,
    'lines.linewidth':     1.5,
    'lines.markersize':    4,
    # ── 图例 ──
    'legend.frameon':      False,
    'legend.borderpad':    0.3,
    'legend.handlelength': 1.5,
    # ── 导出 ──
    'figure.dpi':         150,
    'savefig.dpi':        300,
    'savefig.bbox':       'tight',
    'savefig.pad_inches':  0.05,
    'pdf.fonttype':        42,
    'ps.fonttype':         42,
    # ── 网格 ──
    'axes.grid':           False,
    'grid.alpha':          0.15,
    'grid.linewidth':      0.4,
})

# 路径（脚本位于 true_pre_haic_3/publication_relaxed/；原始 xlsx/csv 在数据根目录）
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR   = os.path.abspath(os.path.join(SCRIPT_DIR, '..', '..'))
PROJECT_ROOT = os.path.normpath(os.path.join(SCRIPT_DIR, '..', '..'))
OUT_DIR    = os.path.join(PROJECT_ROOT, 'output', 'step2_interaction_forest', 'psm02_HAIC_then_I')
os.makedirs(OUT_DIR, exist_ok=True)

MATCHED_IDS_02_CSV = os.path.normpath(os.path.join(
    DATA_DIR, '..', '..', 'HAIC_NO_TACE_4_TIDY', 'update_group_7', 'results',
    'psm_balance_tables_complete', 'matched_ids_02_HAIC_alone_vs_HAIC_then_I.csv'))
MATCHED_IDS_02_LOCAL = os.path.join(SCRIPT_DIR, 'matched_ids_02_HAIC_alone_vs_HAIC_then_I.csv')

# PSM：自洽 TIDY 队列（不读 trace_planB）。baseline = pre-HAIC-1；longitudinal.pre_haic_3 → pre3 与变化率。
_HAIC_DELAYED_DATA_DIR = os.path.normpath(os.path.join(DATA_DIR, '..', 'HAIC_ALONE_AND _DELAYED_IMMUNE'))
TIDY_BASELINE_FOR_PSM_SUPP = [
    os.path.join(SCRIPT_DIR, 'HAIC_NO_TACE_4_TIDY_baseline.csv'),
    os.path.join(_HAIC_DELAYED_DATA_DIR, 'HAIC_NO_TACE_4_TIDY_baseline.csv'),
]
TIDY_LONGITUDINAL_FOR_PSM_SUPP = [
    os.path.join(SCRIPT_DIR, 'HAIC_NO_TACE_4_TIDY_longitudinal.csv'),
    os.path.join(_HAIC_DELAYED_DATA_DIR, 'HAIC_NO_TACE_4_TIDY_longitudinal.csv'),
]

# 与历史 trace_planB_relaxed 同构的分析表列顺序（无参考 CSV 时用）
TRACE_PLANB_COLUMNS_FALLBACK = [
    'patient_id', 'immune_added', 'haic_index_at_immune', 'os_months', 'death_status',
    'afp', 'pivka', 'nlr_bl', 'neut_bl', 'albi_bl', 'alb_bl', 'tbil_bl', 'plr_bl', 'mono_bl',
    'sii_bl', 'piv_bl', 'pvtt_grade', 'hvtt_binary', 'metastasis_binary', 'tumor_count_enc',
    'afp_pre3', 'pivka_pre3', 'nlr_pre3', 'neut_pre3', 'albi_pre3', 'plr_pre3', 'mono_pre3',
    'plt_pre3', 'tbil_pre3', 'alb_pre3', 'sii_pre3', 'piv_pre3',
    'afp_change_pre3', 'pivka_change_pre3', 'nlr_change_pre3', 'neut_change_pre3',
    'albi_change_pre3', 'plr_change_pre3', 'lymph_change_pre3', 'mono_change_pre3',
    'plt_change_pre3', 'sii_change_pre3', 'piv_change_pre3', 'alb_change_pre3', 'tbil_change_pre3',
    'tumor_max_diameter_cm', 'inr', 'plt', 'creatinine', 'ivc_ra_binary', 'ascites_score_enc',
    'log_afp_bl', 'log_pivka_bl', 'lymph_node_binary', 'lymph_bl', 'egv_binary', 'lymph_pre3',
]

OUT_PREFIX = 'psm02_tidymatched_unweighted'

HR_AXIS_LABEL = 'Hazard ratio (unweighted)'
RMST_FIG4_NOTE = '[τ=18m, unweighted KM]'
HR_COL_CSV = 'HR (95% CI)'

# Landmark 敏感性分析参数
LANDMARK_DAYS   = 42             # 2 × 21天（第3次HAIC前）
LANDMARK_MONTHS = LANDMARK_DAYS / 30.44   # ≈ 1.38 个月
TAU_LM          = 16             # landmark后RMST tau（18 - 1.38 ≈ 16.6，取整16）
RMST_LM_SUBCAP = f'[τ={TAU_LM}m from landmark, unweighted KM]'

kmf = KaplanMeierFitter()

def format_p(p):
    if p < 0.001: return '<0.001'
    return f'{p:.3f}'


def _first_existing_path(candidates):
    for p in candidates:
        if p and os.path.isfile(p):
            return os.path.abspath(p)
    return None


def _trace_planb_relaxed_column_order():
    """分析表列名顺序（优先读参考 trace CSV 表头，否则 FALLBACK）。"""
    _cands = [
        os.path.join(SCRIPT_DIR, 'trace_planB_relaxed.csv'),
        os.path.normpath(os.path.join(
            SCRIPT_DIR, '..', 'nonlinear_analysis', 'rms_rcs_relaxed_dual', 'trace_planB_relaxed.csv')),
    ]
    for p in _cands:
        if os.path.isfile(p):
            return list(pd.read_csv(p, nrows=0).columns)
    return list(TRACE_PLANB_COLUMNS_FALLBACK)


def _pct_change_pre3(new_v, old_v):
    """与 00_build_cohort_relaxed_pre3.pct_change 一致（相对变化 %）。"""
    if pd.notna(old_v) and old_v != 0 and pd.notna(new_v):
        return (new_v - old_v) / abs(old_v) * 100.0
    return np.nan


def _tid_baseline_series_to_trace_row(row, pid, psm_treatment_int, pre3_row=None):
    """将 HAIC_NO_TACE_4_TIDY_baseline 映射为 trace_planB_relaxed 列。

    baseline 行视为 pre-HAIC-1；若提供 longitudinal 的 pre_haic_3 行 ``pre3_row``，则按主队列脚本
    计算 pre3 字段及 *_change_pre3（否则 pre3 与变化率为 NaN）。
    haic_index_at_immune 恒为空（本分析不使用 HAIC 周期序号）。
    """
    def num(x):
        v = pd.to_numeric(x, errors='coerce')
        return float(v) if pd.notna(v) else np.nan

    pvtt_raw = str(row.get('pvtt_classification', '') or '').strip().lower()
    pvtt_map = {
        'absent': 0, 'vp1/2': 1, 'vp3/4': 2,
        '无': 0, 'vp1或vp2': 1, 'vp3或vp4': 2,
    }
    pvtt_grade = pvtt_map.get(pvtt_raw, 0)

    asc_raw = str(row.get('ascites', '') or '').strip().lower()
    asc_map = {
        'absent': 0, 'mild': 1, 'moderate-severe': 2,
        '无': 0, '少量': 1, '中-大量': 2, '中或大量': 2,
    }
    ascites_score_enc = asc_map.get(asc_raw, 0)

    def yn_bin(v):
        s = str(v).strip().lower()
        if s in ('yes', '是', '1', '1.0', 'true'):
            return 1.0
        if s in ('no', '否', '0', '0.0', 'false'):
            return 0.0
        return np.nan

    tum_raw = str(row.get('tumor_count_category', '') or '').strip().lower()
    tum_map = {
        'solitary': 0, '>3': 2, '2-3': 1,
        '单发 (1个)': 0, '多发 (>3个)': 2,
    }
    tumor_count_enc = tum_map.get(tum_raw, 0)

    death_raw = str(row.get('death_status', '') or '').strip().lower()
    if death_raw in ('yes', '是', '1', 'dead'):
        death_status = 1
    elif death_raw in ('no', '否', '0', 'alive'):
        death_status = 0
    else:
        death_status = int(num(row.get('death_status'))) if pd.notna(num(row.get('death_status'))) else np.nan

    os_m = num(row.get('os_months'))
    if not np.isfinite(os_m) or os_m <= 0:
        od = num(row.get('os_days'))
        os_m = od / 30.44 if np.isfinite(od) and od > 0 else np.nan

    afp = num(row.get('afp'))
    pivka = num(row.get('pivka'))
    nlr_bl = num(row.get('nlr'))
    neut_bl = num(row.get('neut'))
    alb_bl = num(row.get('alb'))
    tbil_bl = num(row.get('tbil'))
    # Johnson 标准 ALBI：0.66×log₁₀(总胆红素 μmol/L) − 0.085×白蛋白(g/L)（np.log10 = log₁₀）
    # TBil 下限 0.1 仅为避免 log(0)；与 pre3、albi_change_pre3 分母一致。
    albi_bl_lon = (
        0.66 * np.log10(np.maximum(tbil_bl, 0.1)) - 0.085 * alb_bl
        if np.isfinite(tbil_bl) and np.isfinite(alb_bl) else np.nan
    )
    albi_bl = albi_bl_lon if np.isfinite(albi_bl_lon) else num(row.get('albi_score'))
    mono_bl = num(row.get('mono'))
    plt_v = num(row.get('plt'))
    lymph_v = num(row.get('lymph'))
    inr = num(row.get('inr'))
    cr = num(row.get('creatinine'))
    tmax = num(row.get('tumor_max_diameter_cm'))

    plr_bl = plt_v / lymph_v if np.isfinite(lymph_v) and lymph_v > 0 and np.isfinite(plt_v) else np.nan
    sii_bl = plt_v * nlr_bl if np.isfinite(plt_v) and np.isfinite(nlr_bl) else np.nan
    piv_bl = (mono_bl * plt_v * nlr_bl
              if all(np.isfinite(x) for x in (mono_bl, plt_v, nlr_bl)) else np.nan)

    imm = int(psm_treatment_int)
    hi = np.nan  # PSM+TIDY：不按 HAIC 周期定义免疫时机

    log_afp = np.log1p(max(afp, 0)) if np.isfinite(afp) else np.nan
    log_piv = np.log1p(max(pivka, 0)) if np.isfinite(pivka) else np.nan

    out = {
        'patient_id': pid,
        'immune_added': imm,
        'haic_index_at_immune': hi,
        'os_months': os_m,
        'death_status': death_status,
        'afp': afp, 'pivka': pivka,
        'nlr_bl': nlr_bl, 'neut_bl': neut_bl,
        'albi_bl': albi_bl, 'alb_bl': alb_bl, 'tbil_bl': tbil_bl,
        'plr_bl': plr_bl, 'mono_bl': mono_bl, 'sii_bl': sii_bl, 'piv_bl': piv_bl,
        'pvtt_grade': pvtt_grade,
        'hvtt_binary': yn_bin(row.get('hvtt')),
        'metastasis_binary': yn_bin(row.get('distant_metastasis')),
        'tumor_count_enc': tumor_count_enc,
        'tumor_max_diameter_cm': tmax,
        'inr': inr, 'plt': plt_v, 'creatinine': cr,
        'ivc_ra_binary': yn_bin(row.get('ivc_or_ra_thrombus')),
        'ascites_score_enc': ascites_score_enc,
        'log_afp_bl': log_afp, 'log_pivka_bl': log_piv,
        'lymph_node_binary': yn_bin(row.get('lymph_node_metastasis')),
        'lymph_bl': lymph_v,
        'egv_binary': yn_bin(row.get('varices')),
    }
    # pre-HAIC-1 侧衍生（与 00_build 中 pre1 纵向切片一致，用于变化率分母）
    plr_bl_lon = plr_bl
    sii_bl_lon = plt_v * nlr_bl if np.isfinite(plt_v) and np.isfinite(nlr_bl) else np.nan
    piv_bl_lon = (
        mono_bl * plt_v * nlr_bl
        if all(np.isfinite(x) for x in (mono_bl, plt_v, nlr_bl)) else np.nan
    )

    pre3_keys = [
        'afp_pre3', 'pivka_pre3', 'nlr_pre3', 'neut_pre3', 'albi_pre3', 'plr_pre3',
        'mono_pre3', 'plt_pre3', 'tbil_pre3', 'alb_pre3', 'sii_pre3', 'piv_pre3',
        'afp_change_pre3', 'pivka_change_pre3', 'nlr_change_pre3', 'neut_change_pre3',
        'albi_change_pre3', 'plr_change_pre3', 'lymph_change_pre3', 'mono_change_pre3',
        'plt_change_pre3', 'sii_change_pre3', 'piv_change_pre3', 'alb_change_pre3',
        'tbil_change_pre3', 'lymph_pre3',
    ]
    for k in pre3_keys:
        out[k] = np.nan

    if pre3_row is not None:
        pr = pre3_row.iloc[0] if isinstance(pre3_row, pd.DataFrame) else pre3_row
        afp_p3 = num(pr.get('afp'))
        pivka_p3 = num(pr.get('pivka'))
        nlr_p3 = num(pr.get('nlr'))
        neut_p3 = num(pr.get('neut'))
        lymph_p3 = num(pr.get('lymph'))
        plt_p3 = num(pr.get('plt'))
        mono_p3 = num(pr.get('mono'))
        tbil_p3 = num(pr.get('tbil'))
        alb_p3 = num(pr.get('alb'))

        plr_p3 = plt_p3 / lymph_p3 if np.isfinite(lymph_p3) and lymph_p3 > 0 and np.isfinite(plt_p3) else np.nan
        albi_p3 = (
            0.66 * np.log10(np.maximum(tbil_p3, 0.1)) - 0.085 * alb_p3
            if np.isfinite(tbil_p3) and np.isfinite(alb_p3) else np.nan
        )
        sii_p3 = plt_p3 * nlr_p3 if np.isfinite(plt_p3) and np.isfinite(nlr_p3) else np.nan
        piv_p3 = (
            mono_p3 * plt_p3 * nlr_p3
            if all(np.isfinite(x) for x in (mono_p3, plt_p3, nlr_p3)) else np.nan
        )

        out['afp_pre3'] = afp_p3
        out['pivka_pre3'] = pivka_p3
        out['nlr_pre3'] = nlr_p3
        out['neut_pre3'] = neut_p3
        out['lymph_pre3'] = lymph_p3
        out['albi_pre3'] = albi_p3
        out['plr_pre3'] = plr_p3
        out['mono_pre3'] = mono_p3
        out['plt_pre3'] = plt_p3
        out['tbil_pre3'] = tbil_p3
        out['alb_pre3'] = alb_p3
        out['sii_pre3'] = sii_p3
        out['piv_pre3'] = piv_p3

        out['afp_change_pre3'] = _pct_change_pre3(afp_p3, afp)
        out['pivka_change_pre3'] = _pct_change_pre3(pivka_p3, pivka)
        out['nlr_change_pre3'] = _pct_change_pre3(nlr_p3, nlr_bl)
        out['neut_change_pre3'] = _pct_change_pre3(neut_p3, neut_bl)
        out['albi_change_pre3'] = _pct_change_pre3(albi_p3, albi_bl_lon)
        out['plr_change_pre3'] = _pct_change_pre3(plr_p3, plr_bl_lon)
        out['lymph_change_pre3'] = _pct_change_pre3(lymph_p3, lymph_v)
        out['mono_change_pre3'] = _pct_change_pre3(mono_p3, mono_bl)
        out['plt_change_pre3'] = _pct_change_pre3(plt_p3, plt_v)
        out['alb_change_pre3'] = _pct_change_pre3(alb_p3, alb_bl)
        out['tbil_change_pre3'] = _pct_change_pre3(tbil_p3, tbil_bl)
        out['sii_change_pre3'] = _pct_change_pre3(sii_p3, sii_bl_lon)
        out['piv_change_pre3'] = _pct_change_pre3(piv_p3, piv_bl_lon)

    return out


# ============================================================================
# Phase 0: 数据加载（matched_ids_02 + TIDY baseline≡pre1 + longitudinal pre_haic_3）
# ============================================================================
def _normalize_pid(series):
    return series.astype(str).str.strip()


print("=" * 70)
print("加载数据 [PSM matched_ids_02 + TIDY baseline≡pre-HAIC-1 + longitudinal]")
print("=" * 70)

mid_path = MATCHED_IDS_02_CSV if os.path.isfile(MATCHED_IDS_02_CSV) else MATCHED_IDS_02_LOCAL
if not os.path.isfile(mid_path):
    raise FileNotFoundError(
        f"未找到 PSM 匹配 ID 文件。已尝试:\n  {MATCHED_IDS_02_CSV}\n  {MATCHED_IDS_02_LOCAL}\n"
        "请将 matched_ids_02_HAIC_alone_vs_HAIC_then_I.csv 放到 publication_relaxed/ 或 HAIC_NO_TACE_4_TIDY 路径。"
    )
mid = pd.read_csv(mid_path)
if 'patient_id' not in mid.columns:
    raise ValueError(f"{mid_path} 需含列 patient_id")
if 'treatment' not in mid.columns:
    raise ValueError(
        f"{mid_path} 需含列 treatment（0=HAIC alone, 1=HAIC then I）"
    )
mid = mid.copy()
mid['patient_id'] = _normalize_pid(mid['patient_id'])
mid['treatment'] = mid['treatment'].astype(int)

trace_col_order = _trace_planb_relaxed_column_order()
tid_path = _first_existing_path(TIDY_BASELINE_FOR_PSM_SUPP)
if not tid_path:
    raise FileNotFoundError(
        "需要 HAIC_NO_TACE_4_TIDY_baseline.csv，请将文件放在 "
        f"publication_relaxed/ 或 {_HAIC_DELAYED_DATA_DIR}"
    )
long_path = _first_existing_path(TIDY_LONGITUDINAL_FOR_PSM_SUPP)

tid_bl = pd.read_csv(tid_path)
tid_bl['patient_id'] = _normalize_pid(tid_bl['patient_id'])
tid_bl = tid_bl.drop_duplicates(subset=['patient_id'], keep='first').set_index('patient_id')

p3_indexed = None
if long_path:
    tlong = pd.read_csv(long_path)
    tlong['patient_id'] = _normalize_pid(tlong['patient_id'])
    p3_indexed = (
        tlong[tlong['timepoint_type'] == 'pre_haic_3']
        .drop_duplicates(subset=['patient_id'], keep='first')
        .set_index('patient_id')
    )
else:
    print(
        "  注意: 未找到 HAIC_NO_TACE_4_TIDY_longitudinal.csv，"
        "pre3 与变化率列为 NaN。"
    )

rows = []
missing_pid = []
for _, mr in mid.iterrows():
    pid = mr['patient_id']
    trt = int(mr['treatment'])
    if pid not in tid_bl.index:
        missing_pid.append(pid)
        continue
    br = tid_bl.loc[pid]
    if isinstance(br, pd.DataFrame):
        br = br.iloc[0]
    pre3_row = None
    if p3_indexed is not None and pid in p3_indexed.index:
        pre3_row = p3_indexed.loc[pid]
        if isinstance(pre3_row, pd.DataFrame):
            pre3_row = pre3_row.iloc[0]
    sub = _tid_baseline_series_to_trace_row(br, pid, trt, pre3_row=pre3_row)
    rows.append(sub)

if missing_pid:
    raise KeyError(
        f"{len(missing_pid)} 个 PSM patient_id 不在 TIDY baseline（{tid_path}），"
        f"示例: {missing_pid[:5]}"
    )
if len(rows) != len(mid):
    raise ValueError("PSM 行数与 matched_ids 不一致（逻辑错误）")

cohort = pd.DataFrame(rows)
for c in trace_col_order:
    if c not in cohort.columns:
        cohort[c] = np.nan
trace_planB = cohort[trace_col_order].copy().reset_index(drop=True)

n_pre3 = int(trace_planB['nlr_change_pre3'].notna().sum())
_long_msg = f"\n  longitudinal: {long_path}" if long_path else ""
print(
    f"PSM + TIDY: {len(trace_planB)} 人（matched_ids 全数）\n"
    f"  baseline (=pre-HAIC-1): {tid_path}{_long_msg}\n"
    f"  分组 = treatment：HAIC_then_I(1)={int(trace_planB['immune_added'].sum())}，"
    f"HAIC_alone(0)={int((trace_planB['immune_added'] == 0).sum())}\n"
    f"  nlr_change_pre3 非空: {n_pre3}"
)

baseline = pd.read_csv(os.path.join(DATA_DIR, 'HAIC_NO_TACE_4_TIDY_baseline_imputed.csv'))

print(f"队列表（TIDY+matched_ids 骨架）: {len(trace_planB)} 患者")
print(f"  免疫组: {trace_planB['immune_added'].sum()}")
print(f"  HAIC单药: {(trace_planB['immune_added']==0).sum()}")

# ============================================================================
# Phase 1-2: TIDY 分析表 + baseline_imputed 补充列 → full
# ============================================================================
print("\n--- Phase 1-2: 构建分析数据框 ---")

bl_extra = baseline[baseline['patient_id'].isin(trace_planB['patient_id'])].copy()
bl_extra_keep = [c for c in bl_extra.columns
                 if c == 'patient_id' or c not in trace_planB.columns]
bl_extra = bl_extra[bl_extra_keep]

full = trace_planB.merge(bl_extra, on='patient_id', how='left')

full['sex_binary'] = (full['sex'] == '男').astype(float) if 'sex' in full.columns else 0.0
full['afp_high']     = (full['afp'] > 1000).astype(float)
full['pivka_high']   = (full['pivka'] > 8000).astype(float)
full['pvtt_advanced'] = (full['pvtt_grade'] >= 2).astype(int)
full['ascites_binary'] = (full['ascites_score_enc'] > 0).astype(float)
full['tumor_large']    = (full['tumor_max_diameter_cm'] > 10).astype(float)
full['tumor_multiple'] = (full['tumor_count_enc'] >= 1).astype(float)
full['albi_good']      = (full['albi_bl'] <= -2.6).astype(int)
full['albi_grade_pre3'] = full['albi_pre3'].apply(
    lambda x: 1 if pd.notna(x) and x <= -2.6 else (2 if pd.notna(x) and x <= -1.39 else 3))
full['albi_grade_enc'] = full['albi_bl'].apply(
    lambda x: 1 if pd.notna(x) and x <= -2.6 else (2 if pd.notna(x) and x <= -1.39 else 3)).astype(float)
full['nlr_high']       = (full['nlr_bl'] >= 5).astype(int)

full['etiology_hbv'] = (full['etiology'] == 'HBV').astype(float) if 'etiology' in full.columns else 0.0
bclc_map = {'A': 0, 'B': 1, 'C': 2}
full['bclc_enc']      = full['bclc_stage'].map(bclc_map).fillna(1) if 'bclc_stage' in full.columns else 1
full['bclc_advanced'] = (full['bclc_enc'] >= 2).astype(int)
full['dbil'] = pd.to_numeric(full.get('dbil', pd.Series(dtype=float)), errors='coerce')
full['mono_bl'] = pd.to_numeric(full.get('mono', pd.Series(dtype=float)), errors='coerce')
full['hb_bl']  = pd.to_numeric(full.get('hb', pd.Series(dtype=float)), errors='coerce')
full['wbc_bl'] = pd.to_numeric(full.get('wbc', pd.Series(dtype=float)), errors='coerce')
full['cp_grade_enc'] = full['child_pugh_grade'].map({'A': 0, 'B': 1, 'C': 1}).fillna(0) if 'child_pugh_grade' in full.columns else 0
full['cp_score'] = pd.to_numeric(full.get('child_pugh_score', pd.Series(dtype=float)), errors='coerce')

neut_median = full['neut_bl'].median()
plr_median  = full['plr_bl'].median()
inr_median  = full['inr'].median() if 'inr' in full.columns else full.get('inr', pd.Series([1.0])).median()
mono_median = full['mono_bl'].median()
hb_median   = full['hb_bl'].median()
wbc_median  = full['wbc_bl'].median()

full['neut_high'] = (full['neut_bl'] >= neut_median).astype(int)
full['plr_high']  = (full['plr_bl'] >= plr_median).astype(int)
full['inr_high']  = (full['inr'] >= inr_median).astype(int) if 'inr' in full.columns else 0
full['mono_high'] = (full['mono_bl'] >= mono_median).astype(int)
full['hb_low']    = (full['hb_bl'] < hb_median).astype(int)
full['wbc_high']  = (full['wbc_bl'] >= wbc_median).astype(int)

print(f"有效分析患者: {len(full)}")
print(f"  免疫治疗组: {full['immune_added'].sum()}")
print(f"  HAIC单药组: {(full['immune_added']==0).sum()}")
print(f"有pre_haic_3数据: {full['afp_change_pre3'].notna().sum()} 患者")

# ============================================================================
# Phase 4: 非加权（占位列 ps / sw≡1，与旧 Cox 接口兼容）
# ============================================================================
print("\n--- Phase 4: 非加权分析（sw = 1）---")
full['ps'] = np.nan
full['sw'] = 1.0

# ============================================================================
# 全局绘图样式参数（Forest Plot 统一风格）
# ============================================================================
LM_ROW_H      = 0.24    # 每行高度（英寸）
LM_FIG_W      = 16.0    # 图宽（英寸）
LM_FS_HDR     = 12.0    # section header 字体
LM_FS_ROW     = 11.0    # 普通行标签字体
LM_FS_OVERALL = 12.0    # Overall 行字体
LM_FS_N       = 10.5    # n 列字体
LM_FS_TXT     = 10.5    # HR/RMST 数值字体
LM_FS_P       = 9.5     # P 值字体
LM_FS_AXIS    = 11.0    # 轴标签字体
LM_MS_OVERALL = 130     # Overall 散点大小
LM_MS_ROW     = 75      # 普通行散点大小
LM_LW         = 2.2     # CI 线宽

# ============================================================================
# 辅助函数
# ============================================================================

def calc_subgroup_hr(data, sw_col='sw'):
    """非加权 Cox PH（robust SE）。sw_col 仅保留签名兼容。"""
    _ = sw_col
    sub = data[data['os_months'] > 0].copy()
    try:
        cph = CoxPHFitter()
        cph.fit(sub[['os_months', 'death_status', 'immune_added']],
                duration_col='os_months', event_col='death_status', robust=True)
        hr = np.exp(cph.params_['immune_added'])
        ci = np.exp(cph.confidence_intervals_.loc['immune_added'])
        p  = cph.summary.loc['immune_added', 'p']
        return hr, ci.iloc[0], ci.iloc[1], p
    except Exception:
        return np.nan, np.nan, np.nan, np.nan


def _rmst_single(t, e, tau, w=None):
    """KM → RMST（trapz）；w=None 为非加权 KM"""
    k = KaplanMeierFitter()
    if w is not None:
        k.fit(t, e, weights=w)
    else:
        k.fit(t, e)
    sf = k.survival_function_
    tt = sf.index.values
    ss = sf.iloc[:, 0].values
    m  = tt <= tau
    if not m.any():
        return tau * 1.0
    tr = np.append(tt[m], tau)
    sr = np.append(ss[m], ss[m][-1])
    return float(np.trapz(sr, tr))


def calc_rmst_diff(data, tau=18, n_boot=500, sw_col='sw'):
    """两组 RMST 差（非加权 KM）。sw_col 仅保留签名兼容。"""
    _ = sw_col
    sub  = data[data['os_months'] > 0].copy()
    imm  = sub[sub['immune_added'] == 1]
    ctrl = sub[sub['immune_added'] == 0]
    if len(imm) < 1 or len(ctrl) < 1:
        return np.nan, np.nan, np.nan, np.nan

    w1, w0 = None, None
    r1   = _rmst_single(imm['os_months'],  imm['death_status'],  tau, w1)
    r0   = _rmst_single(ctrl['os_months'], ctrl['death_status'], tau, w0)
    diff = r1 - r0

    diffs = []
    rng = np.random.default_rng(42)
    for _ in range(n_boot):
        bi1 = rng.integers(0, len(imm),  size=len(imm))
        bi0 = rng.integers(0, len(ctrl), size=len(ctrl))
        try:
            d = (_rmst_single(imm['os_months'].values[bi1], imm['death_status'].values[bi1],
                              tau, None) -
                 _rmst_single(ctrl['os_months'].values[bi0], ctrl['death_status'].values[bi0],
                              tau, None))
            if not np.isnan(d):
                diffs.append(d)
        except Exception:
            continue

    if len(diffs) == 0:
        return diff, np.nan, np.nan, np.nan
    diffs = np.array(diffs)
    se = np.std(diffs, ddof=1)
    lo = np.percentile(diffs, 2.5)
    hi = np.percentile(diffs, 97.5)
    z  = diff / se if se > 0 else 0.0
    p  = float(2 * (1 - stats.norm.cdf(abs(z))))
    return diff, lo, hi, p


def calc_interaction_p_lrt(data, subgroup_col):
    """修饰因子 × treatment 交互项 LRT（1 df），非加权 Cox。"""
    cols = ['os_months', 'death_status', 'immune_added', subgroup_col]
    sub = data[cols].dropna().copy()
    sub = sub[sub['os_months'] > 0].copy()
    sub['interact'] = sub['immune_added'] * sub[subgroup_col]
    try:
        cph_r = CoxPHFitter()
        cph_f = CoxPHFitter()
        cph_r.fit(sub[['os_months', 'death_status', 'immune_added', subgroup_col]],
                  duration_col='os_months', event_col='death_status')
        cph_f.fit(sub[['os_months', 'death_status', 'immune_added',
                       subgroup_col, 'interact']],
                  duration_col='os_months', event_col='death_status')
        ll_r = cph_r.log_likelihood_
        ll_f = cph_f.log_likelihood_
        lrt_stat = max(2.0 * (ll_f - ll_r), 0.0)
        return float(stats.chi2.sf(lrt_stat, df=1))
    except Exception:
        return np.nan


full_pre3 = full[full['nlr_change_pre3'].notna()].copy()
print(f"有 nlr_change_pre3 数据: {len(full_pre3)} 患者（Landmark 亚组用）")

# ── Fig4 固定切分衍生列（临床预设）；缺失保持 NaN，不误归 0/1 ───────────
def _f4_bin_ge(series, thr):
    s = pd.to_numeric(series, errors='coerce')
    out = pd.Series(np.nan, index=s.index, dtype=float)
    m = s.notna()
    out.loc[m] = (s.loc[m] >= thr).astype(float)
    return out


def _f4_bin_gt(series, thr):
    s = pd.to_numeric(series, errors='coerce')
    out = pd.Series(np.nan, index=s.index, dtype=float)
    m = s.notna()
    out.loc[m] = (s.loc[m] > thr).astype(float)
    return out


def _f4_bin_lt(series, thr):
    s = pd.to_numeric(series, errors='coerce')
    out = pd.Series(np.nan, index=s.index, dtype=float)
    m = s.notna()
    out.loc[m] = (s.loc[m] < thr).astype(float)
    return out



# ============================================================================
# FIGURE 4: Comprehensive Forest Plot（扩展亚组）
# ============================================================================
print("\n--- Figure 4: Forest Plot (expanded subgroups) ---")

# label, col, val, is_calc, is_header, interact_col（第六项保留兼容，P_interaction 见下）
subgroups = [
    ('Overall', None, None, True, False, None),
    ('ALBI Grade 1 (≤ −2.6)', 'albi_good', 1, True, False, None),
    ('ALBI Grade 2–3 (> −2.6)', 'albi_good', 0, True, False, None),
    (f'INR low (< median {inr_median:.2f})', 'inr_high', 0, True, False, None),
    (f'INR high (≥ median {inr_median:.2f})', 'inr_high', 1, True, False, None),
    ('Single tumor', None, None, True, False, None),
    ('2–3 tumors', None, None, True, False, None),
    ('Multiple (>3) tumors', None, None, True, False, None),
    ('BCLC A/B', 'bclc_advanced', 0, True, False, None),
    ('BCLC C', 'bclc_advanced', 1, True, False, None),
    ('PVTT Vp3–4', 'pvtt_advanced', 1, True, False, None),
    ('PVTT absent/Vp1–2', 'pvtt_advanced', 0, True, False, None),
    ('HVTT absent', 'hvtt_binary', 0, True, False, None),
    ('HVTT present', 'hvtt_binary', 1, True, False, None),
    ('IVC/RA thrombus absent', 'ivc_ra_binary', 0, True, False, None),
    ('IVC/RA thrombus present', 'ivc_ra_binary', 1, True, False, None),
    ('No extrahepatic metastasis', 'metastasis_binary', 0, True, False, None),
    ('Extrahepatic metastasis', 'metastasis_binary', 1, True, False, None),
    ('No lymph node metastasis', 'lymph_node_binary', 0, True, False, None),
    ('Lymph node metastasis', 'lymph_node_binary', 1, True, False, None),
    ('No ascites', 'ascites_binary', 0, True, False, None),
    ('Ascites present', 'ascites_binary', 1, True, False, None),
    ('No esophagogastric varices', 'egv_binary', 0, True, False, None),
    ('Esophagogastric varices', 'egv_binary', 1, True, False, None),
]

# 肿瘤数目特殊处理（分层用 tumor_count_enc 三分类；交互 LRT 亦用该列，1df 线性趋势）
tumor_count_filter = {
    'Single tumor':         full['tumor_count_enc'] == 0,
    '2–3 tumors':           full['tumor_count_enc'] == 1,
    'Multiple (>3) tumors': full['tumor_count_enc'] == 2,
}


def _fig4_pint_key(lbl, col_key):
    """Fig4 每行对应的交互检验修饰因子列（与森林分层变量一致）。"""
    if lbl in tumor_count_filter:
        return 'tumor_count_enc'
    return col_key


def _fig4_pair_key(label, col_key):
    """同一分层因子下多行共用的配对键（用于找「ΔRMST 较小」那一行）。"""
    if label == 'Overall':
        return None
    if label in tumor_count_filter:
        return 'tumor_count_enc'
    return col_key


_LM_DUAL_COL_PAIRS = [
    ('opt_afp_bl_high', 'opt_afp_bl_low'),
    ('opt_nlr_bl_high', 'opt_nlr_bl_low'),
    ('opt_albi_bl_poor', 'opt_albi_bl_good'),
    ('opt_tumor_large', 'opt_tumor_small'),
    ('usr_afp_bl_hi', 'usr_afp_bl_lo'),
    ('usr_pivka_bl_hi', 'usr_pivka_bl_lo'),
    ('usr_nlr_bl_hi', 'usr_nlr_bl_lo'),
    ('usr_albi_bl_poor', 'usr_albi_bl_good'),
    ('usr_tumor_lg', 'usr_tumor_sm'),
    ('opt_afp_p3_high', 'opt_afp_p3_low'),
    ('opt_piv_p3_high', 'opt_piv_p3_low'),
    ('opt_nlr_p3_high', 'opt_nlr_p3_low'),
    ('sta_neut_high', 'sta_neut_low'),
    ('opt_albi_p3_poor', 'opt_albi_p3_good'),
    ('usr_afp_p3_hi', 'usr_afp_p3_lo'),
    ('usr_pivka_p3_hi', 'usr_pivka_p3_lo'),
    ('dyn_nlr_exempt', 'dyn_nlr_trigger'),
    ('dyn_afp_exempt', 'dyn_afp_trigger'),
    ('dyn_neut_exempt', 'dyn_neut_trigger'),
    ('dyn_plr_exempt', 'dyn_plr_trigger'),
    ('dyn_pivka_moderate', 'dyn_pivka_deep'),
    ('dyn_albi_stable', 'dyn_albi_improve'),
]
_LM_COL_TO_PAIRKEY = {}
for _i, (_a, _b) in enumerate(_LM_DUAL_COL_PAIRS):
    _k = f'lm_stratum_{_i}'
    _LM_COL_TO_PAIRKEY[_a] = _k
    _LM_COL_TO_PAIRKEY[_b] = _k


def _lm_pair_key(col):
    if col is None:
        return None
    return _LM_COL_TO_PAIRKEY.get(col, col)



def _apply_forest_gray_rules(results, is_overall_fn):
    """
    标记 forest_gray_gfx：亚组 HR 的 P>0.05，或同一配对键内 ΔRMST 较小者；
    该行 HR 与 ΔRMST 图形统一用灰色（由绘图侧读取）。
    """
    by_pair = defaultdict(list)
    for i, row in enumerate(results):
        if not row.get('is_calc') or is_overall_fn(row):
            continue
        pk = row.get('_pair_key')
        rm = row.get('rmst', np.nan)
        if pk is not None and pd.notna(rm):
            by_pair[pk].append((i, float(rm)))
    idx_small_rmst = set()
    for lst in by_pair.values():
        if len(lst) < 2:
            continue
        idx_small_rmst.add(min(lst, key=lambda t: t[1])[0])
    idx_hrp = set()
    for i, row in enumerate(results):
        if not row.get('is_calc') or is_overall_fn(row):
            continue
        hp = row.get('hr_p', np.nan)
        if pd.notna(hp) and hp > 0.05:
            idx_hrp.add(i)
    want = idx_small_rmst | idx_hrp
    for i, row in enumerate(results):
        row['forest_gray_gfx'] = i in want


print("\n--- Fig4 Forest: 非加权 Cox / KM ---")

results = []
# 所有可分层亚组：对 full 队列做 treatment × 修饰因子 Cox LRT（1 df）；同列多行共用同一 P
_f4_pint_keys = set()
for label, col, val, is_calc, is_header, interact_col in subgroups:
    if not is_calc:
        continue
    pk = _fig4_pint_key(label, col)
    if pk is not None:
        _f4_pint_keys.add(pk)

_f4_pint_by_key = {}
print(f"\n--- Fig4: P_interaction（各亚组修饰因子 × treatment，LRT 1 df），共 {len(_f4_pint_keys)} 个修饰变量 ---")
for k in sorted(_f4_pint_keys):
    if k not in full.columns:
        print(f"  P_interaction LRT [{k}] skipped (列不在 full)")
        _f4_pint_by_key[k] = np.nan
        continue
    p = calc_interaction_p_lrt(full, k)
    _f4_pint_by_key[k] = p
    print(f"  P_interaction LRT [{k}] = {format_p(p)}")

for label, col, val, is_calc, is_header, interact_col in subgroups:
    row = {'label': label, 'is_header': is_header, 'is_calc': is_calc, 'is_sep': label == 'SEP'}
    if is_calc:
        # 特殊处理肿瘤数目
        if label in tumor_count_filter:
            sub = full[tumor_count_filter[label]]
        elif col is None:
            sub = full
        else:
            sub = full[full[col] == val]

        nt = (sub['immune_added'] == 1).sum()
        nc = (sub['immune_added'] == 0).sum()
        row['nt'] = nt
        row['nc'] = nc

        pk = _fig4_pint_key(label, col)
        if pk is not None and pk in _f4_pint_by_key:
            row['p_interact'] = _f4_pint_by_key[pk]
        else:
            row['p_interact'] = np.nan

        hr, lo, hi, p = calc_subgroup_hr(sub)
        row.update({'hr': hr, 'hr_lo': lo, 'hr_hi': hi, 'hr_p': p})
        rd, rlo, rhi, rp = calc_rmst_diff(sub, tau=18, n_boot=300)
        row.update({'rmst': rd, 'rmst_lo': rlo, 'rmst_hi': rhi, 'rmst_p': rp})

        if pd.notna(row.get('hr')):
            p_int_str = format_p(row['p_interact']) if pd.notna(row.get('p_interact')) else 'N/A'
            print(f"  {label.strip()}: HR={row['hr']:.2f} ({row['hr_lo']:.2f}-{row['hr_hi']:.2f}), "
                  f"ΔRMST={row.get('rmst', np.nan):.1f}m, P={format_p(row['hr_p'])}, "
                  f"P_int={p_int_str}")
        row['_pair_key'] = _fig4_pair_key(label, col)
    else:
        row['_pair_key'] = None
    results.append(row)

_apply_forest_gray_rules(results, lambda r: r.get('label') == 'Overall')

# ── 绘制 Fig4 Forest Plot（与 SuppFig_lm_forest 统一样式）─────────────────────
# SEP 行压缩 y 坐标，普通行占 1.0 单位
F4_SEP_W  = 0.30
f4_y_pos  = []
f4_y_cur  = 0.0
for _r in reversed(results):
    f4_y_pos.insert(0, f4_y_cur)
    f4_y_cur += F4_SEP_W if _r.get('is_sep') else 1.0
f4_total_y = f4_y_cur

_f4_fig_h = max(10.0, f4_total_y * 0.30 + 2.8)
fig = plt.figure(figsize=(12.5, _f4_fig_h))
# V4: 新增第7列 P_interaction
gs_f = gridspec.GridSpec(1, 7, figure=fig,
                         width_ratios=[2.6, 0.55, 1.05, 1.95, 1.05, 2.10, 1.0], wspace=0.03)
ax_lab      = fig.add_subplot(gs_f[0])
ax_n        = fig.add_subplot(gs_f[1])
ax_hr       = fig.add_subplot(gs_f[2])
ax_hr_txt   = fig.add_subplot(gs_f[3])
ax_rmst     = fig.add_subplot(gs_f[4])
ax_rmst_txt = fig.add_subplot(gs_f[5])
ax_pint     = fig.add_subplot(gs_f[6])   # P_interaction 列

for ax in [ax_lab, ax_n, ax_hr_txt, ax_rmst_txt, ax_pint]:
    ax.set_xlim(0, 1)
    ax.set_ylim(-1.0, f4_total_y)
    ax.axis('off')
for ax in [ax_hr, ax_rmst]:
    ax.set_ylim(-1.0, f4_total_y)
    ax.grid(True, alpha=0.15, axis='x', linewidth=0.6)

ax_lab.text(0.0, f4_total_y + 0.5, 'Subgroup',         fontsize=LM_FS_HDR, fontweight='bold')
ax_n.text(0.5,   f4_total_y + 0.5, 'No. (I/H)',         fontsize=LM_FS_N,   fontweight='bold', ha='center')
ax_hr_txt.text(0.5, f4_total_y + 0.5, 'HR (95% CI)',    fontsize=LM_FS_N,   fontweight='bold', ha='center')
ax_rmst_txt.text(0.5, f4_total_y + 0.5, 'ΔRMST (95% CI)', fontsize=LM_FS_N, fontweight='bold', ha='center')
ax_pint.text(0.5, f4_total_y + 0.5, 'P-interaction',    fontsize=LM_FS_N,   fontweight='bold', ha='center')

for row, y in zip(results, f4_y_pos):
    if row['is_sep']:
        ax_hr.axhline(y=y, color='#e0e0e0', linewidth=0.4)
        ax_rmst.axhline(y=y, color='#e0e0e0', linewidth=0.4)
        ax_pint.axhline(y=y, color='#e0e0e0', linewidth=0.4)
        continue
    if row['is_header']:
        ax_lab.text(0.0, y, row['label'], fontsize=LM_FS_HDR, fontweight='bold',
                    va='center', color='#1565C0')
        continue
    if not row['is_calc']:
        continue

    is_overall = row['label'] == 'Overall'

    if is_overall:
        ax_lab.text(0.0, y, row['label'], fontsize=LM_FS_OVERALL, fontweight='bold', va='center')
    else:
        ax_lab.text(0.02, y, row['label'], fontsize=LM_FS_ROW, va='center', color='black')

    if 'nt' in row:
        ax_n.text(0.5, y, f"{row['nt']}/{row['nc']}", fontsize=LM_FS_N, va='center',
                  ha='center', color='#444')

    p_int = row.get('p_interact', np.nan)
    pint_mc = None
    if not is_overall and pd.notna(p_int):
        if p_int < 0.05:
            pint_mc = P_INTERACT_BRIGHT_RED
        elif p_int < 0.10:
            pint_mc = P_INTERACT_DARK_RED

    hr       = row.get('hr', np.nan)
    rmst_val = row.get('rmst', np.nan)
    gray_gfx = row.get('forest_gray_gfx', False)

    if pd.notna(hr):
        if is_overall:
            mc, ms, mk = 'black', LM_MS_OVERALL, 'D'
        elif gray_gfx:
            mc, ms, mk = CB_GRAY, LM_MS_ROW, 's'
        elif pint_mc is not None:
            mc, ms, mk = pint_mc, LM_MS_ROW, 's'
        elif hr < 0.70:
            mc, ms, mk = CB_GREEN, LM_MS_ROW, 's'
        elif hr > 0.90:
            mc, ms, mk = CB_GRAY, LM_MS_ROW - 5, 's'
        else:
            mc, ms, mk = CB_LBLUE, LM_MS_ROW, 's'

        ax_hr.scatter(hr, y, marker=mk, s=ms, color=mc, edgecolors='black', linewidth=0.7, zorder=5)
        ax_hr.plot([row['hr_lo'], row['hr_hi']], [y, y], color=mc, linewidth=LM_LW, zorder=4)

        fw = 'bold' if is_overall else 'normal'
        ax_hr_txt.text(0.02, y, f"{hr:.2f} ({row['hr_lo']:.2f}–{row['hr_hi']:.2f})",
                       fontsize=LM_FS_TXT, va='center', fontweight=fw)
        ax_hr_txt.text(0.98, y, f"P={format_p(row['hr_p'])}", fontsize=LM_FS_P, va='center',
                       ha='right', color='#222' if row['hr_p'] < 0.05 else '#999')

    if pd.notna(rmst_val):
        if is_overall:
            mc, ms, mk = 'black', LM_MS_OVERALL, 'D'
        elif gray_gfx:
            mc, ms, mk = CB_GRAY, LM_MS_ROW, 's'
        elif pint_mc is not None:
            mc, ms, mk = pint_mc, LM_MS_ROW, 's'
        elif rmst_val > 2.0:
            mc, ms, mk = CB_GREEN, LM_MS_ROW, 's'
        elif rmst_val < 0.5:
            mc, ms, mk = CB_GRAY, LM_MS_ROW - 5, 's'
        else:
            mc, ms, mk = CB_LBLUE, LM_MS_ROW, 's'

        ax_rmst.scatter(rmst_val, y, marker=mk, s=ms, color=mc, edgecolors='black',
                        linewidth=0.7, zorder=5)
        ax_rmst.plot([row['rmst_lo'], row['rmst_hi']], [y, y], color=mc, linewidth=LM_LW, zorder=4)

        fw = 'bold' if is_overall else 'normal'
        ax_rmst_txt.text(0.02, y,
                         f"{rmst_val:+.2f} ({row['rmst_lo']:+.2f}, {row['rmst_hi']:+.2f})",
                         fontsize=LM_FS_TXT, va='center', fontweight=fw)
        ax_rmst_txt.text(0.98, y, f"P={format_p(row['rmst_p'])}", fontsize=LM_FS_P,
                         va='center', ha='right',
                         color='#222' if row['rmst_p'] < 0.05 else '#999')

    # ── P_interaction 列：与 HR/ΔRMST 图形同色（鲜红 <0.05；暗红 0.05≤P<0.1）──
    if pd.notna(p_int):
        p_int_str = format_p(p_int)
        if p_int < 0.05:
            p_int_color = P_INTERACT_BRIGHT_RED
            fw_pi = 'bold'
        elif p_int < 0.10:
            p_int_color = P_INTERACT_DARK_RED
            fw_pi = 'bold'
        else:
            p_int_color = '#888888'
            fw_pi = 'normal'
        ax_pint.text(0.5, y, p_int_str, fontsize=LM_FS_P, va='center', ha='center',
                     color=p_int_color, fontweight=fw_pi)

ax_hr.axvline(x=1.0, color='#333', linestyle='-', linewidth=1.2)
ax_hr.set_xscale('linear')
ax_hr.set_xlim(0.0, 2.0)
ax_hr.set_xlabel(HR_AXIS_LABEL, fontsize=LM_FS_AXIS)
ax_hr.tick_params(axis='x', labelsize=LM_FS_AXIS - 1)
ax_hr.set_yticks([])

ax_rmst.axvline(x=0.0, color='#333', linestyle='-', linewidth=1.2)
ax_rmst.set_xlim(-5, 12)
ax_rmst.set_xlabel('ΔRMST (months)', fontsize=LM_FS_AXIS)
ax_rmst.tick_params(axis='x', labelsize=LM_FS_AXIS - 1)
ax_rmst.set_yticks([])

# Legend — keep only the 4 categorical colors (gray / 2 reds / favors-treatment)
f4_legend = [
    mpatches.Patch(color=CB_GREEN,             label='Favors immunotherapy (HR <0.70)'),
    mpatches.Patch(color=P_INTERACT_BRIGHT_RED, label='P-interaction <0.05'),
    mpatches.Patch(color=P_INTERACT_DARK_RED,   label='P-interaction 0.05–0.10'),
    mpatches.Patch(color=CB_GRAY,               label='No within-subgroup benefit'),
]
# Fixed bottom layout (in inches → figure fraction)
_bot_footnote_in = 0.85   # inches reserved for 3-line footnote
_bot_legend_in   = 0.35   # inches reserved for single-row legend
_bot_xlabel_in   = 0.55   # inches reserved for xlabel + ticks
_bot_total_in    = _bot_footnote_in + _bot_legend_in + _bot_xlabel_in
_top_in          = 0.95   # inches reserved for suptitle + column headers
fig.subplots_adjust(left=0.03, right=0.995,
                    bottom=_bot_total_in / _f4_fig_h,
                    top=(1 - _top_in / _f4_fig_h),
                    wspace=0.03)

# Legend sits in the band between xlabel and footnote
_legend_y_fig = (_bot_footnote_in + _bot_legend_in * 0.5) / _f4_fig_h
fig.legend(handles=f4_legend, fontsize=LM_FS_P, loc='center',
           ncol=4, framealpha=0.0,
           bbox_to_anchor=(0.5, _legend_y_fig),
           handlelength=1.2, handletextpad=0.4, columnspacing=2.0)

fig.suptitle('HAIC then immunotherapy vs HAIC alone — subgroup analysis (PSM-matched)',
             fontsize=13, fontweight='bold', y=(1 - 0.30 / _f4_fig_h))

# Footnote — wrap at ~140 chars to render as 3 readable lines within the figure width
_f4_foot = (
    'I = HAIC then immunotherapy; H = HAIC alone. '
    'HR from unweighted Cox model (robust SE); ΔRMST from unweighted Kaplan–Meier (τ = 18 months); '
    'P-interaction from stratifier × treatment LRT (exploratory; no multiplicity adjustment). '
    'Continuous covariates dichotomized at fixed clinical cutoffs (see Methods).'
)
import textwrap as _tw
fig.text(0.01, 0.005, _tw.fill(_f4_foot, width=140),
         fontsize=7.5, color='#444', style='italic', va='bottom')
plt.savefig(f'{OUT_DIR}/Fig4_forest_v4.pdf', bbox_inches=None)
plt.savefig(f'{OUT_DIR}/Fig4_forest_v4.png', dpi=600, bbox_inches=None)
plt.close()
print("  Figure 4 (v4) saved")

# ============================================================================
# Phase 7: Landmark 敏感性分析（42天，第3次HAIC前）
# 主分析（os_months）完全保留，本模块为独立追加，不修改任何现有代码
# ============================================================================
print("\n" + "=" * 70)
print("Phase 7: Landmark 敏感性分析 (day 42, pre_haic_3)")
print("=" * 70)

# ── 7a: 构建 landmark 专用数据集 ─────────────────────────────────────────────
full_pre3_lm = full_pre3.copy()
full_pre3_lm['os_lm'] = full_pre3_lm['os_months'] - LANDMARK_MONTHS
full_pre3_lm = full_pre3_lm[full_pre3_lm['os_lm'] > 0].copy()
print(f"Landmark 数据集: {len(full_pre3_lm)} 患者 "
      f"(排除 {len(full_pre3) - len(full_pre3_lm)} 例在第42天前截尾)")
print(f"  免疫组: {full_pre3_lm['immune_added'].sum()}, "
      f"HAIC组: {(full_pre3_lm['immune_added']==0).sum()}")

# ── 7e: Landmark Forest Plot（基线特征 + pre-haic-3 动态变量）────────────────
print("\n--- Phase 7e: Landmark Forest Plot (基线 + pre-haic-3 动态变量) ---")

def calc_subgroup_hr_lm(data, time_col='os_lm', sw_col='sw'):
    """Landmark 时间轴非加权 Cox HR。"""
    _ = sw_col
    sub = data[data[time_col] > 0].copy()
    try:
        cph = CoxPHFitter()
        cph.fit(sub[[time_col, 'death_status', 'immune_added']],
                duration_col=time_col, event_col='death_status', robust=True)
        hr = np.exp(cph.params_['immune_added'])
        ci = np.exp(cph.confidence_intervals_.loc['immune_added'])
        p  = cph.summary.loc['immune_added', 'p']
        return hr, ci.iloc[0], ci.iloc[1], p
    except Exception:
        return np.nan, np.nan, np.nan, np.nan

def calc_rmst_diff_lm(data, tau=None, n_boot=300, time_col='os_lm', sw_col='sw'):
    """Landmark 轴两组 RMST 差（非加权 KM）。"""
    _ = sw_col
    if tau is None:
        tau = TAU_LM
    sub = data[data[time_col] > 0].copy()
    imm  = sub[sub['immune_added'] == 1]
    ctrl = sub[sub['immune_added'] == 0]
    if len(imm) < 1 or len(ctrl) < 1:
        return np.nan, np.nan, np.nan, np.nan
    w1, w0 = None, None
    r1   = _rmst_single(imm[time_col], imm['death_status'], tau, w1)
    r0   = _rmst_single(ctrl[time_col], ctrl['death_status'], tau, w0)
    diff = r1 - r0
    diffs = []
    rng = np.random.default_rng(42)
    for _ in range(n_boot):
        i1 = rng.integers(0, len(imm),  size=len(imm))
        i0 = rng.integers(0, len(ctrl), size=len(ctrl))
        try:
            d = _rmst_single(imm[time_col].values[i1], imm['death_status'].values[i1],
                             tau, None) - \
                _rmst_single(ctrl[time_col].values[i0], ctrl['death_status'].values[i0],
                             tau, None)
            if not np.isnan(d):
                diffs.append(d)
        except Exception:
            continue
    if len(diffs) == 0:
        return diff, np.nan, np.nan, np.nan
    diffs = np.array(diffs)
    se = np.std(diffs, ddof=1)
    lo = np.percentile(diffs, 2.5)
    hi = np.percentile(diffs, 97.5)
    z  = diff / se if se > 0 else 0.0
    p  = float(2 * (1 - stats.norm.cdf(abs(z))))
    return diff, lo, hi, p

# 使用已有 full_pre3_lm（含 os_lm，已排除 42 天内截尾）
# 同时需要 pre-haic-3 动态变量，使用 full_pre3_lm（nlr_change_pre3 等已存在）
lm_fp = full_pre3_lm.copy()
# 从 full 传递 baseline 衍生列（已在 Fig4 阶段计算）
for _col in ['f4_neut_ge445']:
    if _col in full.columns and _col not in lm_fp.columns:
        lm_fp[_col] = full.loc[full.index.isin(lm_fp.index), _col]

# Landmark 森林图：连续变量固定切分（与主文 Fig4 一致，无数据驱动最优扫描）
# ── Landmark: 不再对连续变量做固定切分——亚组对齐主图 Fig4（仅分类变量）──

# ── 定义亚组：完全对齐主图 Fig4（仅分类变量；连续变量不再绘制）────────────────
_lm_inr_median = float(lm_fp['inr'].median()) if 'inr' in lm_fp.columns and lm_fp['inr'].notna().any() else float('nan')
lm_tumor_count_filter = {
    'Single tumor':         lm_fp['tumor_count_enc'] == 0,
    '2–3 tumors':           lm_fp['tumor_count_enc'] == 1,
    'Multiple (>3) tumors': lm_fp['tumor_count_enc'] == 2,
}
lm_subgroups = [
    # label, col, val, is_calc, is_header
    ('Overall (Landmark, day 42)', None, None, True, False),
    ('ALBI Grade 1 (≤ −2.6)',   'albi_good', 1, True, False),
    ('ALBI Grade 2–3 (> −2.6)', 'albi_good', 0, True, False),
    (f'INR low (< median {_lm_inr_median:.2f})',  'inr_high', 0, True, False),
    (f'INR high (≥ median {_lm_inr_median:.2f})', 'inr_high', 1, True, False),
    ('Single tumor',         None, None, True, False),
    ('2–3 tumors',           None, None, True, False),
    ('Multiple (>3) tumors', None, None, True, False),
    ('BCLC A/B',             'bclc_advanced', 0, True, False),
    ('BCLC C',               'bclc_advanced', 1, True, False),
    ('PVTT Vp3–4',           'pvtt_advanced', 1, True, False),
    ('PVTT absent/Vp1–2',    'pvtt_advanced', 0, True, False),
    ('HVTT absent',          'hvtt_binary', 0, True, False),
    ('HVTT present',         'hvtt_binary', 1, True, False),
    ('IVC/RA thrombus absent',  'ivc_ra_binary', 0, True, False),
    ('IVC/RA thrombus present', 'ivc_ra_binary', 1, True, False),
    ('No extrahepatic metastasis', 'metastasis_binary', 0, True, False),
    ('Extrahepatic metastasis',    'metastasis_binary', 1, True, False),
    ('No lymph node metastasis', 'lymph_node_binary', 0, True, False),
    ('Lymph node metastasis',    'lymph_node_binary', 1, True, False),
    ('No ascites',     'ascites_binary', 0, True, False),
    ('Ascites present','ascites_binary', 1, True, False),
    ('No esophagogastric varices', 'egv_binary', 0, True, False),
    ('Esophagogastric varices',    'egv_binary', 1, True, False),
]

print("\n--- Landmark Forest: 非加权 Cox / KM ---")

lm_results = []
for label, col, val, is_calc, is_header in lm_subgroups:
    row = {'label': label, 'is_header': is_header, 'is_calc': is_calc, 'is_sep': label == 'SEP'}
    if is_calc:
        if label in lm_tumor_count_filter:
            sub = lm_fp[lm_tumor_count_filter[label]]
        elif col is None:
            sub = lm_fp
        else:
            sub = lm_fp[lm_fp[col] == val]

        sub_valid = sub[sub['os_lm'].notna() & sub['death_status'].notna() & (sub['os_lm'] > 0)]
        nt = (sub_valid['immune_added'] == 1).sum()
        nc = (sub_valid['immune_added'] == 0).sum()
        row['nt'] = nt
        row['nc'] = nc

        try:
            hr, lo, hi, p = calc_subgroup_hr_lm(sub_valid)
            row.update({'hr': hr, 'hr_lo': lo, 'hr_hi': hi, 'hr_p': p})
        except Exception:
            row['hr'] = np.nan
        try:
            rd, rlo, rhi, rp = calc_rmst_diff_lm(sub_valid, tau=TAU_LM, n_boot=300)
            row.update({'rmst': rd, 'rmst_lo': rlo, 'rmst_hi': rhi, 'rmst_p': rp})
        except Exception:
            row['rmst'] = np.nan

        if pd.notna(row.get('hr')):
            print(f"  {label.strip()}: HR={row['hr']:.2f} ({row['hr_lo']:.2f}-{row['hr_hi']:.2f}), "
                  f"ΔRMST={row.get('rmst', np.nan):.2f}m, P={format_p(row['hr_p'])}")
        else:
            print(f"  {label.strip()}: n={nt}/{nc} — HR 不可估计（Cox 失败或数据不足）")
        row['_pair_key'] = _lm_pair_key(col)
    else:
        row['_pair_key'] = None
    lm_results.append(row)

_apply_forest_gray_rules(lm_results, lambda r: 'Overall' in str(r.get('label', '')))

# ── 绘制 Landmark Forest Plot ─────────────────────────────────────────────────
# ── 布局参数已在脚本顶部全局定义（LM_ROW_H / LM_FIG_W / LM_FS_* / LM_MS_* / LM_LW）──

# ── 计算 y 坐标：SEP 行只占 0.45 单位，普通行占 1.0 单位 ──────────────────────
SEP_WEIGHT = 0.30   # SEP 行相对高度权重（压缩分隔线占用空间）
y_positions_lm = []
y_cur = 0.0
for row in reversed(lm_results):
    y_positions_lm.insert(0, y_cur)
    y_cur += SEP_WEIGHT if row.get('is_sep') else 1.0
total_y_span = y_cur

n_rows_lm = len(lm_results)
fig_height = max(12.0, total_y_span * LM_ROW_H)
fig_lm = plt.figure(figsize=(LM_FIG_W, fig_height))
gs_lm = gridspec.GridSpec(1, 6, figure=fig_lm,
                           width_ratios=[3.5, 0.55, 1.3, 1.5, 1.3, 1.5], wspace=0.03)
ax_lm_lab      = fig_lm.add_subplot(gs_lm[0])
ax_lm_n        = fig_lm.add_subplot(gs_lm[1])
ax_lm_hr       = fig_lm.add_subplot(gs_lm[2])
ax_lm_hr_txt   = fig_lm.add_subplot(gs_lm[3])
ax_lm_rmst     = fig_lm.add_subplot(gs_lm[4])
ax_lm_rmst_txt = fig_lm.add_subplot(gs_lm[5])

for ax in [ax_lm_lab, ax_lm_n, ax_lm_hr_txt, ax_lm_rmst_txt]:
    ax.set_xlim(0, 1)
    ax.set_ylim(-1.0, total_y_span)
    ax.axis('off')
for ax in [ax_lm_hr, ax_lm_rmst]:
    ax.set_ylim(-1.0, total_y_span)
    ax.grid(True, alpha=0.15, axis='x', linewidth=0.6)

ax_lm_lab.text(0.0, total_y_span + 0.3, 'Subgroup', fontsize=LM_FS_HDR, fontweight='bold')
ax_lm_n.text(0.5, total_y_span + 0.3, 'I / H', fontsize=LM_FS_N, fontweight='bold', ha='center')
ax_lm_hr_txt.text(0.5, total_y_span + 0.3, 'HR (95% CI)', fontsize=LM_FS_N, fontweight='bold', ha='center')
ax_lm_rmst_txt.text(0.5, total_y_span + 0.3,
                    f'ΔRMST (95% CI)\n{RMST_LM_SUBCAP}',
                    fontsize=LM_FS_N - 0.5, fontweight='bold', ha='center')

# 颜色分区：动态变量豁免线 vs 触发线
EXEMPTION_COLOR = CB_GREEN
TRIGGER_COLOR   = CB_RED
BASELINE_COLOR  = CB_LBLUE

for row, y in zip(lm_results, y_positions_lm):
    if row['is_sep']:
        ax_lm_hr.axhline(y=y, color='#e0e0e0', linewidth=0.4)
        ax_lm_rmst.axhline(y=y, color='#e0e0e0', linewidth=0.4)
        continue
    if row['is_header']:
        # 区分三类 header：基线 / pre-haic-3 静态 / pre-haic-3 动态
        lbl_lower = row['label'].lower()
        if 'static' in lbl_lower or 'cycle 3' in lbl_lower:
            hdr_color = '#6A1B9A'   # 紫色：静态值
        elif 'dynamic' in lbl_lower or 'change rate' in lbl_lower:
            hdr_color = '#37474F'   # 深灰：动态变化率
        else:
            hdr_color = '#1565C0'   # 蓝色：基线
        ax_lm_lab.text(0.0, y, row['label'], fontsize=LM_FS_HDR, fontweight='bold',
                       va='center', color=hdr_color)
        continue
    if not row['is_calc']:
        continue

    is_overall    = 'Overall' in row['label']
    is_exemption  = 'exempt' in row['label'].lower()
    is_trigger    = ('trigger' in row['label'].lower() or
                     'worsening' in row['label'].lower() or
                     'moderate drop' in row['label'].lower() or
                     'stable/rise' in row['label'].lower())
    is_attenuated = ('attenuated' in row['label'].lower() or
                     'improve →' in row['label'].lower())
    is_static_pre3 = 'pre-haic-3' in row['label'].lower()   # pre-haic-3 静态值行
    is_dynamic    = row['label'].strip().startswith(('AFP', 'NLR', 'Neutrophil', 'PLR', 'PIVKA', 'ALBI'))

    # 标签
    if is_overall:
        ax_lm_lab.text(0.0, y, row['label'], fontsize=LM_FS_OVERALL, fontweight='bold', va='center')
    else:
        fc_lbl = (EXEMPTION_COLOR if is_exemption else
                  (TRIGGER_COLOR   if is_trigger   else
                   (CB_GRAY        if is_attenuated else
                    ('#6A1B9A'    if is_static_pre3 else 'black'))))
        ax_lm_lab.text(0.02, y, row['label'], fontsize=LM_FS_ROW, va='center', color=fc_lbl)

    if 'nt' in row:
        ax_lm_n.text(0.5, y, f"{row['nt']}/{row['nc']}", fontsize=LM_FS_N, va='center',
                     ha='center', color='#444')

    hr       = row.get('hr', np.nan)
    rmst_val = row.get('rmst', np.nan)
    gray_gfx = row.get('forest_gray_gfx', False)

    STATIC_PRE3_COLOR = '#6A1B9A'   # 紫色：pre-haic-3 静态值

    if pd.notna(hr):
        if is_overall:
            mc, ms, mk = 'black', LM_MS_OVERALL, 'D'
        elif gray_gfx:
            mc, ms, mk = CB_GRAY, LM_MS_ROW, 's'
        elif is_exemption:
            mc, ms, mk = EXEMPTION_COLOR, LM_MS_ROW, 's'
        elif is_trigger:
            mc, ms, mk = TRIGGER_COLOR, LM_MS_ROW, 's'
        elif is_attenuated:
            mc, ms, mk = CB_GRAY, LM_MS_ROW - 5, 'o'
        elif is_static_pre3:
            mc, ms, mk = STATIC_PRE3_COLOR, LM_MS_ROW, 'D'
        elif hr < 0.70:
            mc, ms, mk = CB_GREEN, LM_MS_ROW, 's'
        elif hr > 0.90:
            mc, ms, mk = CB_GRAY, LM_MS_ROW - 5, 's'
        else:
            mc, ms, mk = BASELINE_COLOR, LM_MS_ROW, 's'

        ax_lm_hr.scatter(hr, y, marker=mk, s=ms, color=mc, edgecolors='black',
                         linewidth=0.7, zorder=5)
        ax_lm_hr.plot([row['hr_lo'], row['hr_hi']], [y, y], color=mc, linewidth=LM_LW, zorder=4)

        fw = 'bold' if is_overall else 'normal'
        ax_lm_hr_txt.text(0.02, y, f"{hr:.2f} ({row['hr_lo']:.2f}–{row['hr_hi']:.2f})",
                          fontsize=LM_FS_TXT, va='center', fontweight=fw)
        ax_lm_hr_txt.text(0.98, y, f"P={format_p(row['hr_p'])}", fontsize=LM_FS_P, va='center',
                          ha='right', color='#222' if row['hr_p'] < 0.05 else '#999')

    if pd.notna(rmst_val):
        if is_overall:
            mc, ms, mk = 'black', LM_MS_OVERALL, 'D'
        elif gray_gfx:
            mc, ms, mk = CB_GRAY, LM_MS_ROW, 's'
        elif is_exemption:
            mc, ms, mk = EXEMPTION_COLOR, LM_MS_ROW, 's'
        elif is_trigger:
            mc, ms, mk = TRIGGER_COLOR, LM_MS_ROW, 's'
        elif is_attenuated:
            mc, ms, mk = CB_GRAY, LM_MS_ROW - 5, 'o'
        elif is_static_pre3:
            mc, ms, mk = STATIC_PRE3_COLOR, LM_MS_ROW, 'D'
        elif rmst_val > 1.5:
            mc, ms, mk = CB_GREEN, LM_MS_ROW, 's'
        elif rmst_val < 0.3:
            mc, ms, mk = CB_GRAY, LM_MS_ROW - 5, 's'
        else:
            mc, ms, mk = BASELINE_COLOR, LM_MS_ROW, 's'

        ax_lm_rmst.scatter(rmst_val, y, marker=mk, s=ms, color=mc, edgecolors='black',
                            linewidth=0.7, zorder=5)
        ax_lm_rmst.plot([row['rmst_lo'], row['rmst_hi']], [y, y], color=mc, linewidth=LM_LW, zorder=4)

        fw = 'bold' if is_overall else 'normal'
        ax_lm_rmst_txt.text(0.02, y,
                             f"{rmst_val:+.2f} ({row['rmst_lo']:+.2f}, {row['rmst_hi']:+.2f})",
                             fontsize=LM_FS_TXT, va='center', fontweight=fw)
        ax_lm_rmst_txt.text(0.98, y, f"P={format_p(row['rmst_p'])}", fontsize=LM_FS_P,
                             va='center', ha='right',
                             color='#222' if row['rmst_p'] < 0.05 else '#999')

# 轴设置
ax_lm_hr.axvline(x=1.0, color='#333', linestyle='-', linewidth=1.2)
ax_lm_hr.set_xscale('linear')
ax_lm_hr.set_xlim(0.0, 2.0)
ax_lm_hr.set_xlabel(HR_AXIS_LABEL, fontsize=LM_FS_AXIS)
ax_lm_hr.tick_params(axis='x', labelsize=LM_FS_AXIS - 1)
ax_lm_hr.set_yticks([])
ax_lm_hr.text(0.5, 1.02, HR_AXIS_LABEL, fontsize=LM_FS_AXIS + 0.5, fontweight='bold',
              ha='center', transform=ax_lm_hr.transAxes, va='bottom')
ax_lm_hr.text(0.02, -0.05, '← Favors\nimmunotherapy', fontsize=LM_FS_AXIS - 1.5, color=CB_GREEN,
              transform=ax_lm_hr.transAxes, ha='left', va='top')
ax_lm_hr.text(0.98, -0.05, 'Favors\nHAIC alone →', fontsize=LM_FS_AXIS - 1.5, color=CB_PINK,
              transform=ax_lm_hr.transAxes, ha='right', va='top')

ax_lm_rmst.axvline(x=0.0, color='#333', linestyle='-', linewidth=1.2)
ax_lm_rmst.set_xlim(-7, 12)
ax_lm_rmst.set_xlabel('ΔRMST (months)', fontsize=LM_FS_AXIS)
ax_lm_rmst.tick_params(axis='x', labelsize=LM_FS_AXIS - 1)
ax_lm_rmst.set_yticks([])
_lm_rmst_title = f'ΔRMST at {TAU_LM}m from landmark (unweighted)'
ax_lm_rmst.text(0.5, 1.02, _lm_rmst_title,
                fontsize=LM_FS_AXIS + 0.5, fontweight='bold',
                ha='center', transform=ax_lm_rmst.transAxes, va='bottom')

# 图例
legend_handles = [
    mpatches.Patch(color=EXEMPTION_COLOR, label='Exempt (deep drop → HAIC alone)'),
    mpatches.Patch(color=TRIGGER_COLOR,   label='Trigger (add immunotherapy)'),
    mpatches.Patch(color=CB_GRAY,         label='Attenuated benefit'),
    mpatches.Patch(color='#6A1B9A',       label='Pre-HAIC-3 static value'),
    mpatches.Patch(color=BASELINE_COLOR,  label='Baseline subgroup'),
]
ax_lm_hr.legend(handles=legend_handles, fontsize=LM_FS_P - 1, loc='upper left',
                framealpha=0.95, edgecolor='#ccc', borderpad=0.4,
                handlelength=1.2, handletextpad=0.4, labelspacing=0.3)

_lm_st_line = (
    f'Unweighted HR (Sandwich SE) and ΔRMST (τ = {TAU_LM} months); '
    f'PSM matched, TIDY cohort'
)
fig_lm.suptitle(
    f'Landmark Sensitivity Analysis (day 42): Forest Plot\n'
    f'Subgroups aligned with main Fig4 (categorical variables only;'
    f' continuous biomarkers are not plotted)\n'
    f'{_lm_st_line}',
    fontsize=12, fontweight='bold', y=1.01
)
plt.tight_layout(rect=[0, 0, 1, 0.99])
plt.savefig(f'{OUT_DIR}/{OUT_PREFIX}_SuppFig_lm_forest.pdf', bbox_inches='tight')
plt.savefig(f'{OUT_DIR}/{OUT_PREFIX}_SuppFig_lm_forest.png', dpi=300, bbox_inches='tight')
plt.close()
print("  SuppFig_lm_forest saved")

# 保存 Landmark Forest Plot 数据表
lm_fp_rows = []
for row in lm_results:
    if row.get('is_sep') or row.get('is_header') or not row.get('is_calc'):
        continue
    if pd.isna(row.get('hr', np.nan)):
        continue
    lm_fp_rows.append({
        'Subgroup': row['label'].strip(),
        'I (n)': row.get('nt', ''),
        'H (n)': row.get('nc', ''),
        HR_COL_CSV: f"{row['hr']:.2f} ({row['hr_lo']:.2f}–{row['hr_hi']:.2f})",
        'P (HR)': format_p(row['hr_p']),
        f'ΔRMST τ={TAU_LM}m (95% CI)': (
            f"{row.get('rmst', np.nan):+.2f} "
            f"({row.get('rmst_lo', np.nan):+.2f}, {row.get('rmst_hi', np.nan):+.2f})"
        ) if pd.notna(row.get('rmst')) else 'N/A',
        'P (RMST)': format_p(row.get('rmst_p', 1)) if pd.notna(row.get('rmst_p')) else 'N/A',
    })
pd.DataFrame(lm_fp_rows).to_csv(f'{OUT_DIR}/{OUT_PREFIX}_Table_lm_forest.csv', index=False)
print("  Table_lm_forest.csv saved")

print("\n" + "=" * 70)
print("完成: Fig4_forest_v4, SuppFig_lm_forest")
print(f"表: {OUT_PREFIX}_Table_lm_forest.csv")
print(f"输出目录: {OUT_DIR}/")
print("=" * 70)

