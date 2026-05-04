#!/usr/bin/env python3
"""
PSM 队列构建：基于 PSM 匹配 ID 筛选患者，合并 TIDY 基线与纵向数据，生成6个 cohort CSV。

PSM ID 来源:
  HAIC_NO_TACE_4_TIDY/update_group_7/results/psm_balance_tables_complete/matched_ids_0X_*.csv

TIDY 数据来源:
  data/immunotherapy_decision_analysis_副本/true_pre_haic_3/nonlinear_analysis/rms_rcs_relaxed_dual/
    - HAIC_NO_TACE_4_TIDY_baseline_imputed.csv
    - HAIC_NO_TACE_4_TIDY_longitudinal.csv

输出（本脚本同目录）:
  cohort_01_HAIC_alone_vs_HAIC_I_conc.csv
  cohort_02_HAIC_alone_vs_HAIC_then_I.csv
  cohort_03_HAIC_alone_vs_HAIC_T_conc.csv
  cohort_04_HAIC_alone_vs_HAIC_then_T.csv
  cohort_05_HAIC_alone_vs_HAIC_IT_conc.csv
  cohort_06_HAIC_alone_vs_HAIC_then_IT.csv

下游: RCS_PSM_dual_timescale.R
"""

import os
import warnings

warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# ── 数据路径 ──────────────────────────────────────────────────────────────────
WORKSPACE = os.path.normpath(os.path.join(SCRIPT_DIR, "..", "..", ".."))
PSM_DIR = os.path.join(
    WORKSPACE,
    "HAIC_NO_TACE_4_TIDY",
    "update_group_7",
    "results",
    "psm_balance_tables_complete",
)
TIDY_DIR = os.path.join(
    SCRIPT_DIR, "..", "..", "data", "tidy_data"
)
PATH_BASELINE = os.path.join(TIDY_DIR, "HAIC_NO_TACE_4_TIDY_baseline_imputed.csv")
PATH_LONG = os.path.join(TIDY_DIR, "HAIC_NO_TACE_4_TIDY_longitudinal.csv")

# ── 6个对子配置 ───────────────────────────────────────────────────────────────
PAIRS = [
    {
        "id": "01",
        "psm_file": "matched_ids_01_HAIC_alone_vs_HAIC+I_concurrent.csv",
        "out_file": "cohort_01_HAIC_alone_vs_HAIC_I_conc.csv",
        "group1": "HAIC_alone",
        "group2": "HAIC+I_concurrent",
        "label": "HAIC_alone vs HAIC+I_concurrent",
    },
    {
        "id": "02",
        "psm_file": "matched_ids_02_HAIC_alone_vs_HAIC_then_I.csv",
        "out_file": "cohort_02_HAIC_alone_vs_HAIC_then_I.csv",
        "group1": "HAIC_alone",
        "group2": "HAIC_then_I",
        "label": "HAIC_alone vs HAIC_then_I",
    },
    {
        "id": "03",
        "psm_file": "matched_ids_03_HAIC_alone_vs_HAIC+T_concurrent.csv",
        "out_file": "cohort_03_HAIC_alone_vs_HAIC_T_conc.csv",
        "group1": "HAIC_alone",
        "group2": "HAIC+T_concurrent",
        "label": "HAIC_alone vs HAIC+T_concurrent",
    },
    {
        "id": "04",
        "psm_file": "matched_ids_04_HAIC_alone_vs_HAIC_then_T.csv",
        "out_file": "cohort_04_HAIC_alone_vs_HAIC_then_T.csv",
        "group1": "HAIC_alone",
        "group2": "HAIC_then_T",
        "label": "HAIC_alone vs HAIC_then_T",
    },
    {
        "id": "05",
        "psm_file": "matched_ids_05_HAIC_alone_vs_HAIC+I+T_concurrent.csv",
        "out_file": "cohort_05_HAIC_alone_vs_HAIC_IT_conc.csv",
        "group1": "HAIC_alone",
        "group2": "HAIC+I+T_concurrent",
        "label": "HAIC_alone vs HAIC+I+T_concurrent",
    },
    {
        "id": "06",
        "psm_file": "matched_ids_06_HAIC_alone_vs_HAIC_then_I+T.csv",
        "out_file": "cohort_06_HAIC_alone_vs_HAIC_then_IT.csv",
        "group1": "HAIC_alone",
        "group2": "HAIC_then_I+T",
        "label": "HAIC_alone vs HAIC_then_I+T",
    },
]

LANDMARK_MONTHS = 42 / 30.44

# ── 辅助函数 ──────────────────────────────────────────────────────────────────
def death_to_event(s: pd.Series) -> pd.Series:
    if s.dtype == object:
        m = s.astype(str).str.strip().str.lower()
        out = np.where(m.isin(("yes", "1", "true", "是")), 1.0, np.nan)
        out = np.where(m.isin(("no", "0", "false", "否")), 0.0, out)
        return pd.Series(out, index=s.index)
    return pd.to_numeric(s, errors="coerce")


def pct_change(new_col, old_col, df):
    new = df[new_col]
    old = df[old_col]
    return np.where(
        (old.notna()) & (old != 0) & (new.notna()),
        (new - old) / old.abs() * 100,
        np.nan,
    )


# ── 加载 TIDY 数据 ────────────────────────────────────────────────────────────
print("=" * 70)
print("build_cohort_psm.py — PSM 队列构建")
print("=" * 70)

print(f"\n[1] 加载基线数据: {PATH_BASELINE}")
baseline = pd.read_csv(PATH_BASELINE)
print(f"    rows={len(baseline)}")

print(f"\n[2] 加载纵向数据: {PATH_LONG}")
long_df = pd.read_csv(PATH_LONG)
print(f"    rows={len(long_df)}")

# ── 纵向：pre_haic_1 (baseline) ───────────────────────────────────────────────
lon_cols = ["patient_id", "afp", "pivka", "nlr", "neut", "lymph", "plt", "mono", "tbil", "alb", "alt", "ast"]
pre1 = long_df.loc[long_df["timepoint_type"] == "baseline", lon_cols].copy()
pre1.columns = [
    "patient_id", "afp_bl_lon", "pivka_bl_lon", "nlr_bl_lon",
    "neut_bl_lon", "lymph_bl_lon", "plt_bl_lon", "mono_bl_lon",
    "tbil_bl_lon", "alb_bl_lon", "alt_bl_lon", "ast_bl_lon",
]
pre1["plr_bl_lon"] = pre1["plt_bl_lon"] / pre1["lymph_bl_lon"].replace(0, np.nan)
pre1["albi_bl_lon"] = 0.66 * np.log10(pre1["tbil_bl_lon"].clip(lower=0.1)) - 0.085 * pre1["alb_bl_lon"]
pre1["sii_bl_lon"] = pre1["plt_bl_lon"] * pre1["nlr_bl_lon"]
pre1["piv_bl_lon"] = pre1["mono_bl_lon"] * pre1["plt_bl_lon"] * pre1["nlr_bl_lon"]

# ── 纵向：pre_haic_3 ──────────────────────────────────────────────────────────
pre3 = long_df.loc[long_df["timepoint_type"] == "pre_haic_3", lon_cols].copy()
pre3.columns = [
    "patient_id", "afp_pre3", "pivka_pre3", "nlr_pre3",
    "neut_pre3", "lymph_pre3", "plt_pre3", "mono_pre3",
    "tbil_pre3", "alb_pre3", "alt_pre3", "ast_pre3",
]
pre3["plr_pre3"] = pre3["plt_pre3"] / pre3["lymph_pre3"].replace(0, np.nan)
pre3["albi_pre3"] = 0.66 * np.log10(pre3["tbil_pre3"].clip(lower=0.1)) - 0.085 * pre3["alb_pre3"]

dyn = pre1.merge(pre3, on="patient_id", how="inner")
print(f"\n    baseline+pre_haic_3 纵向配对: {len(dyn)} 患者")

# 计算变化率
dyn["afp_change_pre3"]   = pct_change("afp_pre3",  "afp_bl_lon",  dyn)
dyn["pivka_change_pre3"] = pct_change("pivka_pre3", "pivka_bl_lon", dyn)
dyn["nlr_change_pre3"]   = pct_change("nlr_pre3",  "nlr_bl_lon",  dyn)
dyn["neut_change_pre3"]  = pct_change("neut_pre3", "neut_bl_lon", dyn)
dyn["albi_change_pre3"]  = pct_change("albi_pre3", "albi_bl_lon", dyn)
dyn["plr_change_pre3"]   = pct_change("plr_pre3",  "plr_bl_lon",  dyn)
dyn["lymph_change_pre3"] = pct_change("lymph_pre3","lymph_bl_lon",dyn)
dyn["mono_change_pre3"]  = pct_change("mono_pre3", "mono_bl_lon", dyn)
dyn["plt_change_pre3"]   = pct_change("plt_pre3",  "plt_bl_lon",  dyn)
dyn["alb_change_pre3"]   = pct_change("alb_pre3",  "alb_bl_lon",  dyn)
dyn["tbil_change_pre3"]  = pct_change("tbil_pre3", "tbil_bl_lon", dyn)
dyn["alt_change_pre3"]   = pct_change("alt_pre3",  "alt_bl_lon",  dyn)
dyn["ast_change_pre3"]   = pct_change("ast_pre3",  "ast_bl_lon",  dyn)
dyn["sii_pre3"]          = dyn["plt_pre3"] * dyn["nlr_pre3"]
dyn["piv_pre3"]          = dyn["mono_pre3"] * dyn["plt_pre3"] * dyn["nlr_pre3"]
dyn["sii_change_pre3"]   = pct_change("sii_pre3",  "sii_bl_lon",  dyn)
dyn["piv_change_pre3"]   = pct_change("piv_pre3",  "piv_bl_lon",  dyn)

DYN_COLS = [
    "patient_id", "afp_pre3", "pivka_pre3", "nlr_pre3", "neut_pre3",
    "lymph_pre3", "albi_pre3", "plr_pre3", "mono_pre3", "plt_pre3",
    "tbil_pre3", "alb_pre3", "alt_pre3", "ast_pre3",
    "sii_pre3", "piv_pre3",
    "afp_change_pre3", "pivka_change_pre3", "nlr_change_pre3",
    "neut_change_pre3", "albi_change_pre3", "plr_change_pre3",
    "lymph_change_pre3", "mono_change_pre3", "plt_change_pre3",
    "sii_change_pre3", "piv_change_pre3", "alb_change_pre3", "tbil_change_pre3",
    "alt_change_pre3", "ast_change_pre3",
]

# ── 基线特征工程（与原 build_cohort.py 完全一致）────────────────────────────
def engineer_baseline(bl: pd.DataFrame) -> pd.DataFrame:
    bl = bl.copy()
    bl["death_status"] = death_to_event(bl["death_status"])
    bl["os_months"] = pd.to_numeric(bl["os_days"], errors="coerce") / 30.44
    bl = bl[(bl["os_months"] > 0) & bl["death_status"].notna()].copy()

    bl["log_afp_bl"]   = np.log1p(pd.to_numeric(bl["afp"],   errors="coerce").clip(lower=0))
    bl["log_pivka_bl"] = np.log1p(pd.to_numeric(bl["pivka"], errors="coerce").clip(lower=0))

    pvtt_raw = bl["pvtt_classification"].astype(str).str.strip().str.lower()
    pvtt_map = {"absent": 0, "vp1/2": 1, "vp3/4": 2, "无": 0}
    bl["pvtt_grade"] = pvtt_raw.map(lambda x: pvtt_map.get(x, 0))

    bl["hvtt_binary"]   = (bl["hvtt"].astype(str).str.strip().str.lower().isin(("yes", "有"))).astype(float)
    bl["ivc_ra_binary"] = (bl["ivc_or_ra_thrombus"].astype(str).str.strip().str.lower().isin(("yes", "有"))).astype(float)
    bl["metastasis_binary"] = (bl["distant_metastasis"].astype(str).str.strip().str.lower().isin(("yes", "是"))).astype(float)

    ascites_map = {
        "absent": 0, "mild": 1, "moderate-severe": 2, "moderate": 2, "severe": 2,
        "无": 0, "少量": 1, "中-大量": 2,
    }
    bl["ascites_score_enc"] = bl["ascites"].map(lambda v: ascites_map.get(str(v).strip().lower().replace(" ", "-"), 0))

    tumor_count_map = {"solitary": 0, "2-3": 1, ">3": 2, "单发 (1个)": 0, "2-3个": 1, "多发 (>3个)": 2}
    bl["tumor_count_enc"] = bl["tumor_count_category"].astype(str).str.strip().str.lower().map(tumor_count_map).fillna(0).astype(int)

    bl["alb_bl"]  = pd.to_numeric(bl["alb"],  errors="coerce")
    bl["tbil_bl"] = pd.to_numeric(bl["tbil"], errors="coerce")
    bl["alt_bl"]  = pd.to_numeric(bl["alt"],  errors="coerce")
    bl["ast_bl"]  = pd.to_numeric(bl["ast"],  errors="coerce")
    _tb = bl["tbil_bl"].clip(lower=0.1)
    bl["albi_bl"] = 0.66 * np.log10(_tb) - 0.085 * bl["alb_bl"]

    bl["nlr_bl"]   = pd.to_numeric(bl["nlr"],   errors="coerce")
    bl["neut_bl"]  = pd.to_numeric(bl["neut"],  errors="coerce")
    bl["lymph_bl"] = pd.to_numeric(bl["lymph"], errors="coerce")
    bl["plr_bl"]   = pd.to_numeric(bl["plt"],   errors="coerce") / bl["lymph_bl"].replace(0, np.nan)
    bl["mono_bl"]  = pd.to_numeric(bl["mono"],  errors="coerce")

    bl["lymph_node_binary"] = (bl["lymph_node_metastasis"].astype(str).str.strip().str.lower().isin(("yes", "是"))).astype(float)
    bl["egv_binary"] = (bl["varices"].astype(str).str.strip().str.lower().isin(("yes", "有"))).astype(float)

    return bl


OUT_COLS = [
    "patient_id", "pair_id", "group1", "group2", "trt",
    "main_group", "os_months", "death_status",
    "afp", "pivka",
    "nlr_bl", "neut_bl", "albi_bl", "alb_bl", "tbil_bl", "alt_bl", "ast_bl",
    "plr_bl", "mono_bl", "sii_bl", "piv_bl",
    "pvtt_grade", "hvtt_binary", "metastasis_binary", "tumor_count_enc",
    "afp_pre3", "pivka_pre3", "nlr_pre3", "neut_pre3", "albi_pre3",
    "plr_pre3", "mono_pre3", "plt_pre3", "tbil_pre3", "alb_pre3",
    "alt_pre3", "ast_pre3",
    "sii_pre3", "piv_pre3",
    "afp_change_pre3", "pivka_change_pre3", "nlr_change_pre3",
    "neut_change_pre3", "albi_change_pre3", "plr_change_pre3",
    "lymph_change_pre3", "mono_change_pre3", "plt_change_pre3",
    "sii_change_pre3", "piv_change_pre3", "alb_change_pre3", "tbil_change_pre3",
    "alt_change_pre3", "ast_change_pre3",
    "tumor_max_diameter_cm", "inr", "plt", "creatinine",
    "ivc_ra_binary", "ascites_score_enc",
    "log_afp_bl", "log_pivka_bl",
    "lymph_node_binary", "lymph_bl", "egv_binary", "lymph_pre3",
]

# ── 对所有基线数据预先做特征工程 ─────────────────────────────────────────────
print("\n[3] 基线特征工程...")
bl_all = engineer_baseline(baseline)
bl_all["sii_bl"] = pd.to_numeric(bl_all["plt"], errors="coerce") * pd.to_numeric(bl_all["nlr_bl"], errors="coerce")
bl_all["piv_bl"] = (
    pd.to_numeric(bl_all["mono_bl"], errors="coerce")
    * pd.to_numeric(bl_all["plt"], errors="coerce")
    * pd.to_numeric(bl_all["nlr_bl"], errors="coerce")
)
# 合并纵向
bl_all = bl_all.merge(dyn[DYN_COLS], on="patient_id", how="left")
print(f"    全队列工程后: {len(bl_all)} 行")

# ── 逐对子构建队列 ────────────────────────────────────────────────────────────
print("\n[4] 按 PSM 匹配 ID 筛选并输出...")
for pair in PAIRS:
    psm_path = os.path.join(PSM_DIR, pair["psm_file"])
    out_path  = os.path.join(SCRIPT_DIR, pair["out_file"])

    psm = pd.read_csv(psm_path)
    # treatment 列：0=group1(HAIC_alone), 1=group2
    matched_ids = psm[["patient_id", "treatment"]].drop_duplicates("patient_id")

    cohort = bl_all.merge(matched_ids, on="patient_id", how="inner")
    cohort["trt"]    = cohort["treatment"].astype(int)
    cohort["pair_id"] = pair["id"]
    cohort["group1"] = pair["group1"]
    cohort["group2"] = pair["group2"]

    available = [c for c in OUT_COLS if c in cohort.columns]
    missing   = [c for c in OUT_COLS if c not in cohort.columns]
    if missing:
        print(f"    [{pair['id']}] 缺少列 (跳过): {missing}")

    cohort[available].to_csv(out_path, index=False)
    n0 = (cohort["trt"] == 0).sum()
    n1 = (cohort["trt"] == 1).sum()
    lm_n = ((cohort["os_months"] - LANDMARK_MONTHS) > 0).sum()
    print(
        f"    [{pair['id']}] {pair['label']}: "
        f"total={len(cohort)} | {pair['group1']}={n0} | {pair['group2']}={n1} | "
        f"landmark(os_lm>0)={lm_n}"
    )

print("\n=== build_cohort_psm.py 完成 ===")
