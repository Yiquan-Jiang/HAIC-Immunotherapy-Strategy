#!/usr/bin/env python3
"""
构建 HAIC_alone 与其余 6 个治疗组的两两配对队列 CSV。

数据源（脚本同目录）:
  - HAIC_NO_TACE_4_TIDY_baseline_imputed.csv
  - HAIC_NO_TACE_4_TIDY_longitudinal.csv  （timepoint_type=baseline 对应原 pre-HAIC-1）

输出（每个配对一个 CSV，同目录）:
  - cohort_HAIC_alone_vs_HAIC_then_I.csv
  - cohort_HAIC_alone_vs_HAIC_then_IT.csv
  - cohort_HAIC_alone_vs_HAIC_then_T.csv
  - cohort_HAIC_alone_vs_HAIC_I_conc.csv
  - cohort_HAIC_alone_vs_HAIC_IT_conc.csv
  - cohort_HAIC_alone_vs_HAIC_T_conc.csv

下游: rcs_all_pairs_dual_timescale.R
"""

import os
import warnings

warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PATH_BASELINE = os.path.join(SCRIPT_DIR, "HAIC_NO_TACE_4_TIDY_baseline_imputed.csv")
PATH_LONG = os.path.join(SCRIPT_DIR, "HAIC_NO_TACE_4_TIDY_longitudinal.csv")

# ── 6 个配对定义 ──────────────────────────────────────────────────────────────
# (对照组名称（main_group 中的值）, 治疗变量列名, 输出 CSV 文件名)
PAIR_CONFIGS = [
    ("HAIC_then_I",         "trt_compare", "cohort_HAIC_alone_vs_HAIC_then_I.csv"),
    ("HAIC_then_I+T",       "trt_compare", "cohort_HAIC_alone_vs_HAIC_then_IT.csv"),
    ("HAIC_then_T",         "trt_compare", "cohort_HAIC_alone_vs_HAIC_then_T.csv"),
    ("HAIC+I_concurrent",   "trt_compare", "cohort_HAIC_alone_vs_HAIC_I_conc.csv"),
    ("HAIC+I+T_concurrent", "trt_compare", "cohort_HAIC_alone_vs_HAIC_IT_conc.csv"),
    ("HAIC+T_concurrent",   "trt_compare", "cohort_HAIC_alone_vs_HAIC_T_conc.csv"),
]

LANDMARK_MONTHS = 42 / 30.44


# ── 工具函数 ──────────────────────────────────────────────────────────────────
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


def map_ascites(val):
    ascites_map = {
        "absent": 0, "mild": 1, "moderate-severe": 2, "moderate": 2, "severe": 2,
        "无": 0, "少量": 1, "中-大量": 2,
    }
    x = str(val).strip().lower().replace(" ", "-")
    return ascites_map.get(x, 0)


# ── 加载原始数据 ──────────────────────────────────────────────────────────────
print("=" * 70)
print("build_all_pairs_cohorts — TIDY baseline (imputed) + longitudinal")
print("=" * 70)

baseline = pd.read_csv(PATH_BASELINE)
print(f"\nLoaded baseline imputed: {len(baseline)} rows")
print("main_group distribution:")
print(baseline["main_group"].value_counts().to_string())

long_df = pd.read_csv(PATH_LONG)
print(f"\nLoaded longitudinal: {len(long_df)} rows")

# ── 构建纵向衍生列（baseline → pre_haic_1；pre3 → pre_haic_3）─────────────────
lon_cols = ["patient_id", "afp", "pivka", "nlr", "neut", "lymph", "plt", "mono", "tbil", "alb", "alt", "ast"]

pre1 = long_df.loc[long_df["timepoint_type"] == "baseline", lon_cols].copy()
pre1.columns = [
    "patient_id", "afp_bl_lon", "pivka_bl_lon", "nlr_bl_lon", "neut_bl_lon",
    "lymph_bl_lon", "plt_bl_lon", "mono_bl_lon", "tbil_bl_lon", "alb_bl_lon",
    "alt_bl_lon", "ast_bl_lon",
]
pre1["plr_bl_lon"] = pre1["plt_bl_lon"] / pre1["lymph_bl_lon"].replace(0, np.nan)
pre1["albi_bl_lon"] = (
    0.66 * np.log10(pre1["tbil_bl_lon"].clip(lower=0.1)) - 0.085 * pre1["alb_bl_lon"]
)
pre1["sii_bl_lon"] = pre1["plt_bl_lon"] * pre1["nlr_bl_lon"]
pre1["piv_bl_lon"] = pre1["mono_bl_lon"] * pre1["plt_bl_lon"] * pre1["nlr_bl_lon"]

# ── 三级回退提取 pre3 等效数据 ──────────────────────────────────────────────
# 优先级: 1) pre_haic_3  2) 离 pre_haic_3 时间最近的 post_haic  3) pre_haic_2
long_df["haic_date"] = pd.to_datetime(long_df["haic_date"])

_pre3_raw = long_df.loc[long_df["timepoint_type"] == "pre_haic_3", ["patient_id", "haic_date"] + lon_cols[1:]].copy()
_pre2_raw = long_df.loc[long_df["timepoint_type"] == "pre_haic_2", ["patient_id", "haic_date"] + lon_cols[1:]].copy()
_post_raw = long_df.loc[long_df["timepoint_type"] == "post_haic", ["patient_id", "haic_date"] + lon_cols[1:]].copy()

# 所有有 pre_haic_3 日期的患者（用作 post_haic 距离计算基准）
pre3_dates = _pre3_raw[["patient_id", "haic_date"]].rename(columns={"haic_date": "pre3_date"})
# 没有 pre_haic_3 的患者，用 pre_haic_2 + 中位周期间隔估算 pre3 日期
_both_dates = _pre2_raw[["patient_id", "haic_date"]].rename(columns={"haic_date": "pre2_date"}).merge(
    pre3_dates, on="patient_id", how="inner"
)
median_cycle_gap = (_both_dates["pre3_date"] - _both_dates["pre2_date"]).dt.days.median()
print(f"\nMedian pre_haic_2 → pre_haic_3 gap: {median_cycle_gap:.0f} days")

# 第 1 级: pre_haic_3
pre3_tier1 = _pre3_raw.drop(columns=["haic_date"]).copy()
tier1_pids = set(pre3_tier1["patient_id"])
print(f"Tier 1 (pre_haic_3): {len(tier1_pids)} patients")

# 第 2 级: 对没有 pre_haic_3 的患者，找离 pre_haic_3 预期日期最近的 post_haic
# 2a: 有 pre_haic_2 的患者 → 用 pre2 + 中位间隔估算 pre3 日期
need_tier2a = _pre2_raw[~_pre2_raw["patient_id"].isin(tier1_pids)][["patient_id", "haic_date"]].rename(
    columns={"haic_date": "pre2_date"}
)
need_tier2a["expected_pre3_date"] = need_tier2a["pre2_date"] + pd.Timedelta(days=median_cycle_gap)

# 2b: 既没有 pre3 也没有 pre2 的患者 → 用 baseline + 2倍中位间隔估算 pre3 日期
_bl_raw = long_df.loc[long_df["timepoint_type"] == "baseline", ["patient_id", "haic_date"]].copy()
pids_no_pre2_no_pre3 = set(_bl_raw["patient_id"]) - tier1_pids - set(need_tier2a["patient_id"])
need_tier2b = _bl_raw[_bl_raw["patient_id"].isin(pids_no_pre2_no_pre3)][["patient_id", "haic_date"]].rename(
    columns={"haic_date": "bl_date"}
)
need_tier2b["expected_pre3_date"] = need_tier2b["bl_date"] + pd.Timedelta(days=median_cycle_gap * 2)

# 合并 2a + 2b 的 expected_pre3_date
need_tier2 = pd.concat([
    need_tier2a[["patient_id", "expected_pre3_date"]],
    need_tier2b[["patient_id", "expected_pre3_date"]],
], ignore_index=True)

post_for_tier2 = _post_raw[_post_raw["patient_id"].isin(need_tier2["patient_id"])].copy()
post_for_tier2 = post_for_tier2.merge(
    need_tier2[["patient_id", "expected_pre3_date"]], on="patient_id"
)
post_for_tier2["abs_gap"] = (post_for_tier2["haic_date"] - post_for_tier2["expected_pre3_date"]).dt.days.abs()
# 取距离最近的一条
post_for_tier2 = post_for_tier2.sort_values("abs_gap").groupby("patient_id").first().reset_index()
pre3_tier2 = post_for_tier2.drop(columns=["haic_date", "expected_pre3_date", "abs_gap"]).copy()
tier2_pids = set(pre3_tier2["patient_id"])
print(f"Tier 2 (nearest post_haic): {len(tier2_pids)} patients "
      f"(from pre2-based={len(need_tier2a)}, from bl-based={len(need_tier2b)})")

# 第 3 级: pre_haic_2（有 pre2 但 Tier2 没匹配到 post_haic 的）
covered = tier1_pids | tier2_pids
pre3_tier3 = _pre2_raw[~_pre2_raw["patient_id"].isin(covered)].drop(columns=["haic_date"]).copy()
tier3_pids = set(pre3_tier3["patient_id"])
print(f"Tier 3 (pre_haic_2): {len(tier3_pids)} patients")

# 合并三级
pre3 = pd.concat([pre3_tier1, pre3_tier2, pre3_tier3], ignore_index=True)
pre3.columns = [
    "patient_id", "afp_pre3", "pivka_pre3", "nlr_pre3", "neut_pre3",
    "lymph_pre3", "plt_pre3", "mono_pre3", "tbil_pre3", "alb_pre3",
    "alt_pre3", "ast_pre3",
]
pre3["plr_pre3"] = pre3["plt_pre3"] / pre3["lymph_pre3"].replace(0, np.nan)
pre3["albi_pre3"] = (
    0.66 * np.log10(pre3["tbil_pre3"].clip(lower=0.1)) - 0.085 * pre3["alb_pre3"]
)
print(f"Total pre3-equivalent: {len(pre3)} patients "
      f"(Tier1={len(tier1_pids)}, Tier2={len(tier2_pids)}, Tier3={len(tier3_pids)})")

dyn = pre1.merge(pre3, on="patient_id", how="inner")
print(f"Patients with baseline + pre3-equivalent in longitudinal: {len(dyn)}")

dyn["afp_change_pre3"]   = pct_change("afp_pre3",   "afp_bl_lon",   dyn)
dyn["pivka_change_pre3"] = pct_change("pivka_pre3",  "pivka_bl_lon", dyn)
dyn["nlr_change_pre3"]   = pct_change("nlr_pre3",    "nlr_bl_lon",   dyn)
dyn["neut_change_pre3"]  = pct_change("neut_pre3",   "neut_bl_lon",  dyn)
dyn["albi_change_pre3"]  = pct_change("albi_pre3",   "albi_bl_lon",  dyn)
dyn["plr_change_pre3"]   = pct_change("plr_pre3",    "plr_bl_lon",   dyn)
dyn["lymph_change_pre3"] = pct_change("lymph_pre3",  "lymph_bl_lon", dyn)
dyn["mono_change_pre3"]  = pct_change("mono_pre3",   "mono_bl_lon",  dyn)
dyn["plt_change_pre3"]   = pct_change("plt_pre3",    "plt_bl_lon",   dyn)
dyn["alb_change_pre3"]   = pct_change("alb_pre3",    "alb_bl_lon",   dyn)
dyn["tbil_change_pre3"]  = pct_change("tbil_pre3",   "tbil_bl_lon",  dyn)
dyn["alt_change_pre3"]   = pct_change("alt_pre3",    "alt_bl_lon",   dyn)
dyn["ast_change_pre3"]   = pct_change("ast_pre3",    "ast_bl_lon",   dyn)
dyn["sii_pre3"]          = dyn["plt_pre3"] * dyn["nlr_pre3"]
dyn["piv_pre3"]          = dyn["mono_pre3"] * dyn["plt_pre3"] * dyn["nlr_pre3"]
dyn["sii_change_pre3"]   = pct_change("sii_pre3", "sii_bl_lon", dyn)
dyn["piv_change_pre3"]   = pct_change("piv_pre3", "piv_bl_lon", dyn)

dyn_cols = [
    "patient_id",
    "afp_pre3", "pivka_pre3", "nlr_pre3", "neut_pre3", "lymph_pre3",
    "albi_pre3", "plr_pre3", "mono_pre3", "plt_pre3", "tbil_pre3", "alb_pre3",
    "alt_pre3", "ast_pre3",
    "sii_pre3", "piv_pre3",
    "afp_change_pre3", "pivka_change_pre3", "nlr_change_pre3", "neut_change_pre3",
    "albi_change_pre3", "plr_change_pre3", "lymph_change_pre3", "mono_change_pre3",
    "plt_change_pre3", "sii_change_pre3", "piv_change_pre3",
    "alb_change_pre3", "tbil_change_pre3",
    "alt_change_pre3", "ast_change_pre3",
]


# ── 构建基线衍生列（一次性处理全量 baseline）────────────────────────────────
def build_baseline_features(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["death_status"] = death_to_event(df["death_status"])
    df["os_months"] = pd.to_numeric(df["os_days"], errors="coerce") / 30.44
    df = df[(df["os_months"] > 0) & df["death_status"].notna()].copy()

    df["sex_binary"] = (
        df["sex"].astype(str).str.strip().str.lower().isin(("male", "男", "m"))
    ).astype(float)
    df["log_afp_bl"]   = np.log1p(pd.to_numeric(df["afp"],   errors="coerce").clip(lower=0))
    df["log_pivka_bl"] = np.log1p(pd.to_numeric(df["pivka"], errors="coerce").clip(lower=0))

    pvtt_raw = df["pvtt_classification"].astype(str).str.strip().str.lower()
    pvtt_map = {"absent": 0, "vp1/2": 1, "vp3/4": 2, "无": 0}
    df["pvtt_grade"] = pvtt_raw.map(lambda x: pvtt_map.get(x, 0))

    df["hvtt_binary"] = (
        df["hvtt"].astype(str).str.strip().str.lower().isin(("yes", "有"))
    ).astype(float)
    df["ivc_ra_binary"] = (
        df["ivc_or_ra_thrombus"].astype(str).str.strip().str.lower().isin(("yes", "有"))
    ).astype(float)
    df["metastasis_binary"] = (
        df["distant_metastasis"].astype(str).str.strip().str.lower().isin(("yes", "是"))
    ).astype(float)
    df["ascites_score_enc"] = df["ascites"].map(map_ascites)

    tumor_count_map = {
        "solitary": 0, "2-3": 1, ">3": 2,
        "单发 (1个)": 0, "2-3个": 1, "多发 (>3个)": 2,
    }
    tc = df["tumor_count_category"].astype(str).str.strip().str.lower()
    df["tumor_count_enc"] = tc.map(tumor_count_map).fillna(0).astype(int)

    df["alb_bl"]  = pd.to_numeric(df["alb"],  errors="coerce")
    df["tbil_bl"] = pd.to_numeric(df["tbil"], errors="coerce")
    df["alt_bl"]  = pd.to_numeric(df["alt"],  errors="coerce")
    df["ast_bl"]  = pd.to_numeric(df["ast"],  errors="coerce")
    # ALBI: log10(TBIL μmol/L)×0.66 − ALB(g/L)×0.085
    _tb = df["tbil_bl"].clip(lower=0.1)
    df["albi_bl"]       = 0.66 * np.log10(_tb) - 0.085 * df["alb_bl"]
    df["albi_grade_enc"] = pd.to_numeric(df["albi_grade"], errors="coerce")
    df["nlr_bl"]  = pd.to_numeric(df["nlr"],  errors="coerce")
    df["neut_bl"] = pd.to_numeric(df["neut"], errors="coerce")
    df["lymph_bl"] = pd.to_numeric(df["lymph"], errors="coerce")
    df["plr_bl"]  = pd.to_numeric(df["plt"], errors="coerce") / df["lymph_bl"].replace(0, np.nan)
    df["mono_bl"] = pd.to_numeric(df["mono"], errors="coerce")
    df["sii_bl"]  = pd.to_numeric(df["plt"], errors="coerce") * df["nlr_bl"]
    df["piv_bl"]  = df["mono_bl"] * pd.to_numeric(df["plt"], errors="coerce") * df["nlr_bl"]

    df["lymph_node_binary"] = (
        df["lymph_node_metastasis"].astype(str).str.strip().str.lower().isin(("yes", "是"))
    ).astype(float)
    df["egv_binary"] = (
        df["varices"].astype(str).str.strip().str.lower().isin(("yes", "有"))
    ).astype(float)

    cp_map = {"a": 0, "b": 1, "c": 1}
    df["cp_grade_enc"] = (
        df["child_pugh_grade"].astype(str).str.strip().str.lower().map(cp_map).fillna(0)
    )
    df["cp_score"] = pd.to_numeric(df["child_pugh_score"], errors="coerce")

    bclc_map = {"a": 0, "b": 1, "c": 2}
    df["bclc_enc"] = (
        df["bclc_stage"].astype(str).str.strip().str.upper().map(bclc_map).fillna(1)
    )
    return df


baseline_feat = build_baseline_features(baseline)
print(f"\nBaseline after feature engineering: {len(baseline_feat)} rows (os>0, death not NA)")

# ── 输出列清单 ────────────────────────────────────────────────────────────────
OUT_COLS = [
    "patient_id", "main_group", "trt_compare",
    "os_months", "death_status",
    "afp", "pivka",
    "nlr_bl", "neut_bl", "albi_bl", "alb_bl", "tbil_bl", "alt_bl", "ast_bl",
    "plr_bl", "mono_bl", "sii_bl", "piv_bl",
    "pvtt_grade", "hvtt_binary", "metastasis_binary", "tumor_count_enc",
    "afp_pre3", "pivka_pre3", "nlr_pre3", "neut_pre3", "albi_pre3",
    "plr_pre3", "mono_pre3", "plt_pre3", "tbil_pre3", "alb_pre3",
    "alt_pre3", "ast_pre3",
    "sii_pre3", "piv_pre3",
    "afp_change_pre3", "pivka_change_pre3", "nlr_change_pre3", "neut_change_pre3",
    "albi_change_pre3", "plr_change_pre3", "lymph_change_pre3", "mono_change_pre3",
    "plt_change_pre3", "sii_change_pre3", "piv_change_pre3",
    "alb_change_pre3", "tbil_change_pre3",
    "alt_change_pre3", "ast_change_pre3",
    "tumor_max_diameter_cm", "inr", "plt", "creatinine",
    "ivc_ra_binary", "ascites_score_enc",
    "log_afp_bl", "log_pivka_bl",
    "lymph_node_binary", "lymph_bl", "egv_binary", "lymph_pre3",
]


# ── 逐配对构建并保存 ──────────────────────────────────────────────────────────
print("\n" + "=" * 70)
for compare_group, trt_col, out_fname in PAIR_CONFIGS:
    print(f"\n>>> 配对: HAIC_alone vs {compare_group}")

    mask = baseline_feat["main_group"].isin(("HAIC_alone", compare_group))
    bl = baseline_feat[mask].copy()
    bl["trt_compare"] = (bl["main_group"] == compare_group).astype(int)

    full = bl.merge(dyn[dyn_cols], on="patient_id", how="left")

    n_chg = full["nlr_change_pre3"].notna().sum()
    print(
        f"    Merged: {len(full)} rows | "
        f"HAIC_alone: {(full['trt_compare']==0).sum()} | "
        f"{compare_group}: {full['trt_compare'].sum()} | "
        f"nlr_change_pre3 non-NA: {n_chg}"
    )

    available = [c for c in OUT_COLS if c in full.columns]
    missing   = [c for c in OUT_COLS if c not in full.columns]
    if missing:
        print(f"    警告: 缺少列 (跳过): {missing}")

    out_path = os.path.join(SCRIPT_DIR, out_fname)
    full[available].to_csv(out_path, index=False)

    # 打印 landmark 统计
    os_lm = full["os_months"] - LANDMARK_MONTHS
    lm = full[os_lm > 0]
    print(
        f"    Saved: {out_fname}  rows={len(full)}  cols={len(available)}"
    )
    print(
        f"    Landmark os_lm>0: {len(lm)} | "
        f"alone: {(lm['trt_compare']==0).sum()} | "
        f"{compare_group}: {lm['trt_compare'].sum()}"
    )

print("\n" + "=" * 70)
print("Done. 所有配对队列 CSV 已保存至:", SCRIPT_DIR)
