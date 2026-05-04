#!/usr/bin/env python3
"""
RCS 分析队列：main_group 为 HAIC_alone 与 HAIC+I+T_concurrent（TIDY 插补基线 + TIDY 纵向）。

数据源（脚本同目录）:
  - HAIC_NO_TACE_4_TIDY_baseline_imputed.csv
  - HAIC_NO_TACE_4_TIDY_longitudinal.csv  （纵向中 timepoint_type=baseline 对应原 pre-HAIC-1）

输出: data/RCS_INT_HAIC_ALONE_AND_I_T_CONC_cohort.csv

下游: RCS_INT_HAIC_ALONE_AND_I_T_CONC_dual_timescale.R
"""

import os
import warnings

warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, "data")
PATH_BASELINE = os.path.join(DATA_DIR, "HAIC_NO_TACE_4_TIDY_baseline_imputed.csv")
PATH_LONG = os.path.join(DATA_DIR, "HAIC_NO_TACE_4_TIDY_longitudinal.csv")
OUT_PATH = os.path.join(DATA_DIR, "RCS_INT_HAIC_ALONE_AND_I_T_CONC_cohort.csv")

MAIN_GROUPS = ("HAIC_alone", "HAIC+I+T_concurrent")
TRT_COL = "trt_haic_i_t_conc"
TRT_VAL = "HAIC+I+T_concurrent"


def death_to_event(s: pd.Series) -> pd.Series:
    if s.dtype == object:
        m = s.astype(str).str.strip().str.lower()
        out = np.where(m.isin(("yes", "1", "true", "是")), 1.0, np.nan)
        out = np.where(m.isin(("no", "0", "false", "否")), 0.0, out)
        return pd.Series(out, index=s.index)
    return pd.to_numeric(s, errors="coerce")


print("=" * 70)
print("RCS_INT_HAIC_ALONE_AND_I_T_CONC — cohort from TIDY baseline (imputed) + longitudinal")
print("=" * 70)

baseline = pd.read_csv(PATH_BASELINE)
print(f"\nLoaded baseline imputed: {len(baseline)} rows")

bl = baseline[baseline["main_group"].isin(MAIN_GROUPS)].copy()
bl[TRT_COL] = (bl["main_group"] == TRT_VAL).astype(int)

bl["death_status"] = death_to_event(bl["death_status"])
bl["os_months"] = pd.to_numeric(bl["os_days"], errors="coerce") / 30.44
bl = bl[(bl["os_months"] > 0) & bl["death_status"].notna()].copy()
print(f"After main_group filter + os>0 + death: {len(bl)}")

# ── 与旧 trace 对齐的派生列（TIDY 为英文分类）────────────────────────────
bl["sex_binary"] = (bl["sex"].astype(str).str.strip().str.lower().isin(("male", "男", "m"))).astype(float)
bl["log_afp_bl"] = np.log1p(pd.to_numeric(bl["afp"], errors="coerce").clip(lower=0))
bl["log_pivka_bl"] = np.log1p(pd.to_numeric(bl["pivka"], errors="coerce").clip(lower=0))

pvtt_raw = bl["pvtt_classification"].astype(str).str.strip().str.lower()
pvtt_map = {"absent": 0, "vp1/2": 1, "vp3/4": 2, "无": 0}
bl["pvtt_grade"] = pvtt_raw.map(lambda x: pvtt_map.get(x, 0))
bl["hvtt_binary"] = (bl["hvtt"].astype(str).str.strip().str.lower().isin(("yes", "有"))).astype(float)
bl["ivc_ra_binary"] = (
    bl["ivc_or_ra_thrombus"].astype(str).str.strip().str.lower().isin(("yes", "有"))
).astype(float)
bl["metastasis_binary"] = (
    bl["distant_metastasis"].astype(str).str.strip().str.lower().isin(("yes", "是"))
).astype(float)

ascites_map = {
    "absent": 0, "mild": 1, "moderate-severe": 2, "moderate": 2, "severe": 2,
    "无": 0, "少量": 1, "中-大量": 2,
}


def map_ascites(val):
    x = str(val).strip().lower().replace(" ", "-")
    return ascites_map.get(x, 0)


bl["ascites_score_enc"] = bl["ascites"].map(map_ascites)

tumor_count_map = {"solitary": 0, "2-3": 1, ">3": 2, "单发 (1个)": 0, "2-3个": 1, "多发 (>3个)": 2}
tc = bl["tumor_count_category"].astype(str).str.strip().str.lower()
bl["tumor_count_enc"] = tc.map(tumor_count_map).fillna(0).astype(int)

bl["alb_bl"] = pd.to_numeric(bl["alb"], errors="coerce")
bl["tbil_bl"] = pd.to_numeric(bl["tbil"], errors="coerce")
_tb = bl["tbil_bl"].clip(lower=0.1)
bl["albi_bl"] = 0.66 * np.log10(_tb) - 0.085 * bl["alb_bl"]
bl["albi_grade_enc"] = pd.to_numeric(bl["albi_grade"], errors="coerce")
bl["nlr_bl"] = pd.to_numeric(bl["nlr"], errors="coerce")
bl["neut_bl"] = pd.to_numeric(bl["neut"], errors="coerce")
bl["lymph_bl"] = pd.to_numeric(bl["lymph"], errors="coerce")
bl["plr_bl"] = pd.to_numeric(bl["plt"], errors="coerce") / bl["lymph_bl"].replace(0, np.nan)

bl["lymph_node_binary"] = (
    bl["lymph_node_metastasis"].astype(str).str.strip().str.lower().isin(("yes", "是"))
).astype(float)
bl["egv_binary"] = (bl["varices"].astype(str).str.strip().str.lower().isin(("yes", "有"))).astype(float)
bl["mono_bl"] = pd.to_numeric(bl["mono"], errors="coerce")

cp_map = {"a": 0, "b": 1, "c": 1}
bl["cp_grade_enc"] = bl["child_pugh_grade"].astype(str).str.strip().str.lower().map(cp_map).fillna(0)
bl["cp_score"] = pd.to_numeric(bl["child_pugh_score"], errors="coerce")

bclc_map = {"a": 0, "b": 1, "c": 2}
bl["bclc_enc"] = bl["bclc_stage"].astype(str).str.strip().str.upper().map(bclc_map).fillna(1)

# ── 纵向：baseline = 原 pre_haic_1 ───────────────────────────────────────
print(f"\nLoaded longitudinal: {PATH_LONG}")
long_df = pd.read_csv(PATH_LONG)

lon_cols = ["patient_id", "afp", "pivka", "nlr", "neut", "lymph", "plt", "mono", "tbil", "alb"]
pre1 = long_df.loc[long_df["timepoint_type"] == "baseline", lon_cols].copy()
pre1.columns = [
    "patient_id",
    "afp_bl_lon",
    "pivka_bl_lon",
    "nlr_bl_lon",
    "neut_bl_lon",
    "lymph_bl_lon",
    "plt_bl_lon",
    "mono_bl_lon",
    "tbil_bl_lon",
    "alb_bl_lon",
]
pre1["plr_bl_lon"] = pre1["plt_bl_lon"] / pre1["lymph_bl_lon"].replace(0, np.nan)
pre1["albi_bl_lon"] = 0.66 * np.log10(pre1["tbil_bl_lon"].clip(lower=0.1)) - 0.085 * pre1["alb_bl_lon"]
pre1["sii_bl_lon"] = pre1["plt_bl_lon"] * pre1["nlr_bl_lon"]
pre1["piv_bl_lon"] = pre1["mono_bl_lon"] * pre1["plt_bl_lon"] * pre1["nlr_bl_lon"]

pre3 = long_df.loc[long_df["timepoint_type"] == "pre_haic_3", lon_cols].copy()
pre3.columns = [
    "patient_id",
    "afp_pre3",
    "pivka_pre3",
    "nlr_pre3",
    "neut_pre3",
    "lymph_pre3",
    "plt_pre3",
    "mono_pre3",
    "tbil_pre3",
    "alb_pre3",
]
pre3["plr_pre3"] = pre3["plt_pre3"] / pre3["lymph_pre3"].replace(0, np.nan)
pre3["albi_pre3"] = 0.66 * np.log10(pre3["tbil_pre3"].clip(lower=0.1)) - 0.085 * pre3["alb_pre3"]

dyn = pre1.merge(pre3, on="patient_id", how="inner")
print(f"Patients with baseline+pre_haic_3 in longitudinal: {len(dyn)}")


def pct_change(new_col, old_col, df):
    new = df[new_col]
    old = df[old_col]
    return np.where((old.notna()) & (old != 0) & (new.notna()), (new - old) / old.abs() * 100, np.nan)


dyn["afp_change_pre3"] = pct_change("afp_pre3", "afp_bl_lon", dyn)
dyn["pivka_change_pre3"] = pct_change("pivka_pre3", "pivka_bl_lon", dyn)
dyn["nlr_change_pre3"] = pct_change("nlr_pre3", "nlr_bl_lon", dyn)
dyn["neut_change_pre3"] = pct_change("neut_pre3", "neut_bl_lon", dyn)
dyn["albi_change_pre3"] = pct_change("albi_pre3", "albi_bl_lon", dyn)
dyn["plr_change_pre3"] = pct_change("plr_pre3", "plr_bl_lon", dyn)
dyn["lymph_change_pre3"] = pct_change("lymph_pre3", "lymph_bl_lon", dyn)
dyn["mono_change_pre3"] = pct_change("mono_pre3", "mono_bl_lon", dyn)
dyn["plt_change_pre3"] = pct_change("plt_pre3", "plt_bl_lon", dyn)
dyn["alb_change_pre3"] = pct_change("alb_pre3", "alb_bl_lon", dyn)
dyn["tbil_change_pre3"] = pct_change("tbil_pre3", "tbil_bl_lon", dyn)
dyn["sii_pre3"] = dyn["plt_pre3"] * dyn["nlr_pre3"]
dyn["piv_pre3"] = dyn["mono_pre3"] * dyn["plt_pre3"] * dyn["nlr_pre3"]
dyn["sii_change_pre3"] = pct_change("sii_pre3", "sii_bl_lon", dyn)
dyn["piv_change_pre3"] = pct_change("piv_pre3", "piv_bl_lon", dyn)

dyn_cols = [
    "patient_id",
    "afp_pre3",
    "pivka_pre3",
    "nlr_pre3",
    "neut_pre3",
    "lymph_pre3",
    "albi_pre3",
    "plr_pre3",
    "mono_pre3",
    "plt_pre3",
    "tbil_pre3",
    "alb_pre3",
    "sii_pre3",
    "piv_pre3",
    "afp_change_pre3",
    "pivka_change_pre3",
    "nlr_change_pre3",
    "neut_change_pre3",
    "albi_change_pre3",
    "plr_change_pre3",
    "lymph_change_pre3",
    "mono_change_pre3",
    "plt_change_pre3",
    "sii_change_pre3",
    "piv_change_pre3",
    "alb_change_pre3",
    "tbil_change_pre3",
]

full = bl.merge(dyn[dyn_cols], on="patient_id", how="left")
full["sii_bl"] = pd.to_numeric(full["plt"], errors="coerce") * pd.to_numeric(full["nlr_bl"], errors="coerce")
full["piv_bl"] = (
    pd.to_numeric(full["mono_bl"], errors="coerce")
    * pd.to_numeric(full["plt"], errors="coerce")
    * pd.to_numeric(full["nlr_bl"], errors="coerce")
)

n_chg = full["nlr_change_pre3"].notna().sum()
print(f"Merged cohort: {len(full)} | nlr_change_pre3 non-NA: {n_chg}")

out_cols = [
    "patient_id",
    "main_group",
    TRT_COL,
    "os_months",
    "death_status",
    "afp",
    "pivka",
    "nlr_bl",
    "neut_bl",
    "albi_bl",
    "alb_bl",
    "tbil_bl",
    "plr_bl",
    "mono_bl",
    "sii_bl",
    "piv_bl",
    "pvtt_grade",
    "hvtt_binary",
    "metastasis_binary",
    "tumor_count_enc",
    "afp_pre3",
    "pivka_pre3",
    "nlr_pre3",
    "neut_pre3",
    "albi_pre3",
    "plr_pre3",
    "mono_pre3",
    "plt_pre3",
    "tbil_pre3",
    "alb_pre3",
    "sii_pre3",
    "piv_pre3",
    "afp_change_pre3",
    "pivka_change_pre3",
    "nlr_change_pre3",
    "neut_change_pre3",
    "albi_change_pre3",
    "plr_change_pre3",
    "lymph_change_pre3",
    "mono_change_pre3",
    "plt_change_pre3",
    "sii_change_pre3",
    "piv_change_pre3",
    "alb_change_pre3",
    "tbil_change_pre3",
    "tumor_max_diameter_cm",
    "inr",
    "plt",
    "creatinine",
    "ivc_ra_binary",
    "ascites_score_enc",
    "log_afp_bl",
    "log_pivka_bl",
    "lymph_node_binary",
    "lymph_bl",
    "egv_binary",
    "lymph_pre3",
]

available = [c for c in out_cols if c in full.columns]
missing = [c for c in out_cols if c not in full.columns]
if missing:
    print("  警告: 缺少列 (跳过):", missing)

full[available].to_csv(OUT_PATH, index=False)
print(f"\nSaved: {OUT_PATH}  rows={len(full)}  cols={len(available)}")
print(f"  HAIC_alone: {(full[TRT_COL] == 0).sum()}  |  HAIC+I+T_concurrent: {full[TRT_COL].sum()}")

LANDMARK_MONTHS = 42 / 30.44
full["os_lm"] = full["os_months"] - LANDMARK_MONTHS
lm = full[full["os_lm"] > 0]
print(f"\nLandmark os_lm>0: {len(lm)}  |  alone: {(lm[TRT_COL]==0).sum()}  |  I_T_conc: {lm[TRT_COL].sum()}")
print("Done.")