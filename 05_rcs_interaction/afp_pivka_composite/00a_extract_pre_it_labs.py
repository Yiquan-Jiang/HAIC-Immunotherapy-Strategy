#!/usr/bin/env python3
"""
提取 pre-IT（免疫/靶向治疗前最近一次）实验室检查数据。

从原始 GBK 编码检验 CSV 中提取肝功能、炎症、血细胞指标，
匹配到距离首次免疫/靶向治疗开始日期最近的时间点（30天内）。

规则：
  - 序贯组（THEN_*）：取 first_it_date 前 30 天内最近的检查
  - 并发组（*_CONC）：取 first_haic_date 前 30 天内最近的检查（= baseline）
  - HAIC_alone 对照组：用所有时间点的平均值

包含 AFP / PIVKA 的 pre-IT 提取及变化率计算。

参考代码：scripts_clean/01_data_prep/10_add_post_haic_timepoints.py (load_lab_data)
下游：01_rcs_afp_pivka_composite.R
"""

import os
import warnings

warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.normpath(os.path.join(SCRIPT_DIR, "..", "..", ".."))
BASE_PATH = os.path.normpath(os.path.join(PROJECT_ROOT, ".."))

SWIMMER_CSV = os.path.join(
    PROJECT_ROOT, "data", "publication_cohorts", "00_swimmer_plot_events.csv"
)

RAW_BIOCHEM = os.path.join(BASE_PATH, "基线数据.基线_生化检查_gbk.csv")
RAW_HEMA = os.path.join(BASE_PATH, "基线数据.基线_血常规_gbk.csv")
RAW_COAG = os.path.join(BASE_PATH, "基线数据.基线_凝血检查_gbk.csv")
RAW_TUMOR = os.path.join(BASE_PATH, "基线数据.基线_肿瘤标记物检查_gbk.csv")
RAW_PIVKA = os.path.join(BASE_PATH, "基线数据.基线_PIVKA_gbk.csv")

COHORT_FILES = {
    "THEN_IT": {"csv": "composite_THEN_IT_cohort.csv", "type": "sequential"},
    "THEN_I": {"csv": "composite_THEN_I_cohort.csv", "type": "sequential"},
    "THEN_T": {"csv": "composite_THEN_T_cohort.csv", "type": "sequential"},
    "T_CONC": {"csv": "composite_T_CONC_cohort.csv", "type": "concurrent"},
    "I_CONC": {"csv": "composite_I_CONC_cohort.csv", "type": "concurrent"},
    "IT_CONC": {"csv": "composite_IT_CONC_cohort.csv", "type": "concurrent"},
}

PRE_IT_WINDOW_DAYS = 30

# ── 原始数据加载 ─────────────────────────────────────────────────────────────


def load_lab_data(file_path, item_keywords, col_names, target_patients=None):
    """从 GBK 原始 CSV 加载实验室数据（复用 10_add_post_haic_timepoints.py 的模式）"""
    if not os.path.exists(file_path):
        print(f"  WARNING: {file_path} not found")
        return {}

    df = pd.read_csv(file_path, encoding="gbk", skiprows=2, low_memory=False)
    df = df[df["item_name"].notna()].copy()
    if target_patients is not None:
        df = df[df["patient_sn"].isin(target_patients)].copy()
    df["test_time"] = pd.to_datetime(df["test_time"], errors="coerce")
    df = df[df["test_time"].notna()].copy()

    results = {}
    for keyword, col_name in zip(item_keywords, col_names):
        item_df = df[df["item_name"].str.contains(keyword, na=False)].copy()
        item_df["value"] = pd.to_numeric(item_df["test_result"], errors="coerce")
        item_df = item_df[item_df["value"].notna()].copy()
        item_df = item_df.rename(
            columns={"patient_sn": "patient_id", "test_time": "date"}
        )
        results[col_name] = item_df[["patient_id", "date", "value"]].copy()
    return results


def load_all_lab_data(target_patients):
    """加载所有需要的原始检验数据"""
    print("\n[加载原始检验数据]")
    all_lab = {}

    print("  - 肝功能 (生化检查)...")
    biochem = load_lab_data(
        RAW_BIOCHEM,
        ["总胆红素", "白蛋白", "丙氨酸氨基转移酶", "天门冬氨酸氨基转移酶"],
        ["tbil", "alb", "alt", "ast"],
        target_patients,
    )
    all_lab.update(biochem)

    print("  - 血常规...")
    hema = load_lab_data(
        RAW_HEMA,
        ["中性粒细胞计数", "淋巴细胞计数", "单核细胞计数", "血小板计数"],
        ["neut", "lymph", "mono", "plt"],
        target_patients,
    )
    all_lab.update(hema)

    print("  - 凝血检查...")
    coag = load_lab_data(
        RAW_COAG, ["凝血酶原时间", "INR"], ["pt", "inr"], target_patients
    )
    all_lab.update(coag)

    print("  - 肿瘤标志物 (AFP)...")
    tumor = load_lab_data(
        RAW_TUMOR, ["甲胎蛋白"], ["afp"], target_patients
    )
    all_lab.update(tumor)

    print("  - 肿瘤标志物 (PIVKA)...")
    pivka = load_lab_data(
        RAW_PIVKA, ["异常凝血酶原"], ["pivka"], target_patients
    )
    all_lab.update(pivka)

    for name, lab_df in all_lab.items():
        print(f"    {name}: {len(lab_df)} records, {lab_df['patient_id'].nunique()} patients")

    return all_lab


# ── 匹配函数 ─────────────────────────────────────────────────────────────────


def find_nearest_before(patient_id, target_date, lab_df, max_days=30):
    """在 target_date 前 max_days 天内找最近的检查值。

    返回 (value, lab_date) 或 (NaN, NaT)
    """
    if lab_df is None or len(lab_df) == 0:
        return np.nan, pd.NaT

    pdata = lab_df[lab_df["patient_id"] == patient_id]
    if len(pdata) == 0:
        return np.nan, pd.NaT

    target_dt = pd.to_datetime(target_date)
    pdata = pdata.copy()
    pdata["date"] = pd.to_datetime(pdata["date"]).dt.normalize()
    target_dt_norm = target_dt.normalize()

    # 仅取 target_date 前 max_days 天到 target_date 当天的检查
    pdata["days_before"] = (target_dt_norm - pdata["date"]).dt.days
    valid = pdata[(pdata["days_before"] >= 0) & (pdata["days_before"] <= max_days)]

    if len(valid) == 0:
        return np.nan, pd.NaT

    nearest = valid.loc[valid["days_before"].idxmin()]
    return nearest["value"], nearest["date"]


def extract_pre_it_for_patient(patient_id, target_date, all_lab, lab_cols):
    """为单个患者在 target_date 前提取所有检验指标"""
    row = {"patient_id": patient_id}
    dates = []

    for col in lab_cols:
        lab_df = all_lab.get(col)
        val, dt = find_nearest_before(patient_id, target_date, lab_df, PRE_IT_WINDOW_DAYS)
        row[f"{col}_pre_it"] = val
        if pd.notna(dt):
            dates.append(dt)

    # 取所有匹配检查中最近的日期作为 pre_it_lab_date
    if dates:
        row["pre_it_lab_date"] = max(dates)
        row["pre_it_days_before_it"] = (
            pd.to_datetime(target_date).normalize() - pd.to_datetime(max(dates))
        ).days
    else:
        row["pre_it_lab_date"] = pd.NaT
        row["pre_it_days_before_it"] = np.nan

    return row


def compute_patient_mean(patient_id, all_lab, lab_cols):
    """为 HAIC_alone 患者计算所有时间点的平均值。"""
    row = {"patient_id": patient_id}

    for col in lab_cols:
        lab_df = all_lab.get(col)
        if lab_df is None or len(lab_df) == 0:
            row[f"{col}_pre_it"] = np.nan
            continue

        pdata = lab_df[lab_df["patient_id"] == patient_id]
        if len(pdata) == 0:
            row[f"{col}_pre_it"] = np.nan
            continue

        vals = pd.to_numeric(pdata["value"], errors="coerce").dropna()
        row[f"{col}_pre_it"] = vals.mean() if len(vals) > 0 else np.nan

    row["pre_it_lab_date"] = pd.NaT
    row["pre_it_days_before_it"] = np.nan
    return row


# ── 派生指标 ─────────────────────────────────────────────────────────────────


def compute_derived(df):
    """计算 NLR, PLR, SII, PIV, ALBI 的 pre_it 版本"""
    neut = df["neut_pre_it"]
    lymph = df["lymph_pre_it"]
    mono = df["mono_pre_it"]
    plt_ = df["plt_pre_it"]
    tbil = df["tbil_pre_it"]
    alb = df["alb_pre_it"]

    df["nlr_pre_it"] = neut / lymph.replace(0, np.nan)
    df["plr_pre_it"] = plt_ / lymph.replace(0, np.nan)
    df["sii_pre_it"] = plt_ * neut / lymph.replace(0, np.nan)
    df["piv_pre_it"] = mono * plt_ * neut / lymph.replace(0, np.nan)
    # ALBI: 0.66 * log10(tbil_umol_L) - 0.085 * alb_g_L
    df["albi_pre_it"] = np.where(
        tbil.notna() & alb.notna() & (tbil > 0),
        0.66 * np.log10(tbil) - 0.085 * alb,
        np.nan,
    )
    return df


def add_pre_it_composites(df):
    """为 pre-IT 的 AFP/PIVKA 添加 4 个组合指标，
    与 00_build_composite_cohorts.py 中 baseline / pre3 的 add_static_composites 逻辑对齐。
    输出列: log_afp_pivka_ratio_pre_it, log_afp_pivka_sum_pre_it,
            afp_pivka_pc1_pre_it,       afp_pivka_pc2_pre_it
    """
    afp_vals = pd.to_numeric(df.get("afp_pre_it"), errors="coerce")
    pivka_vals = pd.to_numeric(df.get("pivka_pre_it"), errors="coerce")

    log_afp = np.log1p(afp_vals.clip(lower=0))
    log_pivka = np.log1p(pivka_vals.clip(lower=0))

    df["log_afp_pivka_ratio_pre_it"] = log_afp - log_pivka
    df["log_afp_pivka_sum_pre_it"] = log_afp + log_pivka

    mask = log_afp.notna() & log_pivka.notna()
    X = np.column_stack([log_afp.values, log_pivka.values])
    valid = mask.values & np.isfinite(X).all(axis=1)

    pc1 = np.full(len(df), np.nan)
    pc2 = np.full(len(df), np.nan)

    if valid.sum() >= 10:
        scaler = StandardScaler()
        X_valid = scaler.fit_transform(X[valid])
        pca = PCA(n_components=2)
        scores = pca.fit_transform(X_valid)
        pc1[valid] = scores[:, 0]
        pc2[valid] = scores[:, 1]

        loadings = pca.components_
        explained = pca.explained_variance_ratio_
        print(f"  PCA_pre_it: loadings={np.round(loadings, 3).tolist()}, "
              f"var_explained={np.round(explained, 3).tolist()}")

    df["afp_pivka_pc1_pre_it"] = pc1
    df["afp_pivka_pc2_pre_it"] = pc2
    return df


def add_pre_it_dynamic_composites(df):
    """为 pre-IT 的 AFP/PIVKA 变化率添加 4 个动态组合指标，
    与 00_build_composite_cohorts.py 中 pre3 动态组合（afp_pivka_change_diff / sum /
    pc1_dyn / pc2_dyn）的 add_dynamic_composites 逻辑对齐。
    输出列: afp_pivka_change_diff_pre_it, afp_pivka_change_sum_pre_it,
            afp_pivka_pc1_dyn_pre_it,    afp_pivka_pc2_dyn_pre_it
    """
    afp_chg = pd.to_numeric(df.get("afp_change_pre_it"), errors="coerce")
    pivka_chg = pd.to_numeric(df.get("pivka_change_pre_it"), errors="coerce")

    df["afp_pivka_change_diff_pre_it"] = afp_chg - pivka_chg
    df["afp_pivka_change_sum_pre_it"] = afp_chg + pivka_chg

    mask = afp_chg.notna() & pivka_chg.notna()
    X = np.column_stack([afp_chg.values, pivka_chg.values])
    valid = mask.values & np.isfinite(X).all(axis=1)

    pc1 = np.full(len(df), np.nan)
    pc2 = np.full(len(df), np.nan)

    if valid.sum() >= 10:
        scaler = StandardScaler()
        X_valid = scaler.fit_transform(X[valid])
        pca = PCA(n_components=2)
        scores = pca.fit_transform(X_valid)
        pc1[valid] = scores[:, 0]
        pc2[valid] = scores[:, 1]

        loadings = pca.components_
        explained = pca.explained_variance_ratio_
        print(f"  PCA_dyn_pre_it: loadings={np.round(loadings, 3).tolist()}, "
              f"var_explained={np.round(explained, 3).tolist()}")

    df["afp_pivka_pc1_dyn_pre_it"] = pc1
    df["afp_pivka_pc2_dyn_pre_it"] = pc2
    return df


def compute_change_rates(df):
    """计算 pre_it 相对于 baseline 的变化率（%）"""
    # baseline 列名映射：composite CSV 中 baseline 列名
    bl_map = {
        "alb": "alb_bl",
        "tbil": "tbil_bl",
        "alt": "alt_bl",
        "ast": "ast_bl",
        "albi": "albi_bl",
        "nlr": "nlr_bl",
        "plr": "plr_bl",
        "sii": "sii_bl",
        "piv": "piv_bl",
        "neut": "neut_bl",
        "lymph": "lymph_bl",
        "mono": "mono_bl",
        "plt": "plt",  # baseline plt 列名是 "plt"
        "afp": "afp",  # baseline afp 列名是 "afp"
        "pivka": "pivka",  # baseline pivka 列名是 "pivka"
    }

    for var, bl_col in bl_map.items():
        pre_it_col = f"{var}_pre_it"
        change_col = f"{var}_change_pre_it"
        if pre_it_col in df.columns and bl_col in df.columns:
            bl_vals = pd.to_numeric(df[bl_col], errors="coerce")
            pre_it_vals = pd.to_numeric(df[pre_it_col], errors="coerce")
            df[change_col] = np.where(
                bl_vals.notna() & pre_it_vals.notna() & (bl_vals != 0),
                (pre_it_vals - bl_vals) / bl_vals.abs() * 100,
                np.nan,
            )
        else:
            df[change_col] = np.nan
    return df


# ── 主流程 ───────────────────────────────────────────────────────────────────


def main():
    print("=" * 70)
    print("Pre-IT Lab Value Extraction")
    print("=" * 70)

    # 1. 加载 swimmer events，提取首次 I/T 日期
    print("\n[1] 加载 swimmer events，提取首次 I/T 日期...")
    sw = pd.read_csv(SWIMMER_CSV)
    sw["start_date"] = pd.to_datetime(sw["start_date"])
    it_events = sw[sw["treatment_category"].isin(["Immunotherapy", "Targeted Therapy"])]
    first_it = (
        it_events.groupby("patient_id")["start_date"].min().reset_index()
    )
    first_it.columns = ["patient_id", "first_it_date"]
    print(f"  共 {len(first_it)} 名患者有 I/T 记录")

    # 获取首次 HAIC 日期（用于并发组 baseline 时间点）
    haic_events = sw[sw["treatment_category"].isin(["HAIC", "HAIC+TACE"])]
    first_haic = (
        haic_events.groupby("patient_id")["start_date"].min().reset_index()
    )
    first_haic.columns = ["patient_id", "first_haic_date"]
    print(f"  共 {len(first_haic)} 名患者有 HAIC 记录")

    # 2. 收集所有需要查询的患者
    all_patients = set()
    for key, info in COHORT_FILES.items():
        csv_path = os.path.join(SCRIPT_DIR, info["csv"])
        if os.path.exists(csv_path):
            df = pd.read_csv(csv_path, usecols=["patient_id"])
            all_patients.update(df["patient_id"].unique())
    print(f"  所有队列共 {len(all_patients)} 名唯一患者")

    # 3. 加载原始检验数据
    raw_lab_cols = ["tbil", "alb", "alt", "ast", "neut", "lymph", "mono", "plt", "pt", "inr", "afp", "pivka"]
    all_lab = load_all_lab_data(all_patients)

    # 4. 处理每个队列
    for key, info in COHORT_FILES.items():
        csv_path = os.path.join(SCRIPT_DIR, info["csv"])
        if not os.path.exists(csv_path):
            print(f"\n  SKIP: {csv_path} not found")
            continue

        print(f"\n{'='*70}")
        print(f"[{key}] {info['csv']} (type={info['type']})")
        print(f"{'='*70}")

        df = pd.read_csv(csv_path)
        n_total = len(df)
        is_trt = df["main_group"] != "HAIC_alone"

        # 删除上一轮可能残留的 pre_it 列，避免 merge 冲突
        drop_cols = [c for c in df.columns if c.endswith("_pre_it") or c in
                     ["first_it_date", "first_haic_date", "pre_it_source",
                      "pre_it_lab_date", "pre_it_days_before_it"]]
        df = df.drop(columns=drop_cols, errors="ignore")

        # 合并 I/T 和 HAIC 日期
        df = df.merge(first_it, on="patient_id", how="left")
        df = df.merge(first_haic, on="patient_id", how="left")

        # 分三种情况处理
        pre_it_rows = []
        source_counts = {"pre_it_matched": 0, "baseline_concurrent": 0,
                         "haic_alone_mean": 0, "no_match": 0}

        for idx, row in df.iterrows():
            pid = row["patient_id"]
            is_treatment = is_trt.iloc[idx]

            if not is_treatment:
                # HAIC_alone：用所有时间点的平均值
                pre_it_row = compute_patient_mean(pid, all_lab, raw_lab_cols)
                pre_it_row["pre_it_source"] = "haic_alone_mean"
                pre_it_row["first_it_date"] = pd.NaT
                source_counts["haic_alone_mean"] += 1

            elif info["type"] == "concurrent":
                # 并发组：用 first_haic_date 前的 baseline
                target_date = row.get("first_haic_date")
                if pd.notna(target_date):
                    pre_it_row = extract_pre_it_for_patient(
                        pid, target_date, all_lab, raw_lab_cols
                    )
                    pre_it_row["pre_it_source"] = "baseline_concurrent"
                    pre_it_row["first_it_date"] = row.get("first_it_date")
                    source_counts["baseline_concurrent"] += 1
                else:
                    pre_it_row = {"patient_id": pid, "pre_it_source": "no_match"}
                    for col in raw_lab_cols:
                        pre_it_row[f"{col}_pre_it"] = np.nan
                    pre_it_row["first_it_date"] = pd.NaT
                    pre_it_row["pre_it_lab_date"] = pd.NaT
                    pre_it_row["pre_it_days_before_it"] = np.nan
                    source_counts["no_match"] += 1

            else:
                # 序贯组：用 first_it_date 前 30 天
                target_date = row.get("first_it_date")
                if pd.notna(target_date):
                    pre_it_row = extract_pre_it_for_patient(
                        pid, target_date, all_lab, raw_lab_cols
                    )
                    pre_it_row["pre_it_source"] = "pre_it_matched"
                    pre_it_row["first_it_date"] = target_date
                    # 检查是否实际匹配到数据
                    has_any = any(
                        pd.notna(pre_it_row.get(f"{c}_pre_it"))
                        for c in raw_lab_cols
                    )
                    if not has_any:
                        pre_it_row["pre_it_source"] = "no_match"
                        source_counts["no_match"] += 1
                    else:
                        source_counts["pre_it_matched"] += 1
                else:
                    pre_it_row = {"patient_id": pid, "pre_it_source": "no_match"}
                    for col in raw_lab_cols:
                        pre_it_row[f"{col}_pre_it"] = np.nan
                    pre_it_row["first_it_date"] = pd.NaT
                    pre_it_row["pre_it_lab_date"] = pd.NaT
                    pre_it_row["pre_it_days_before_it"] = np.nan
                    source_counts["no_match"] += 1

            pre_it_rows.append(pre_it_row)

        pre_it_df = pd.DataFrame(pre_it_rows)

        print(f"\n  Source distribution:")
        for src, cnt in source_counts.items():
            print(f"    {src}: {cnt}")

        # 合并回原始数据
        # 先删除 df 中可能已有的 pre_it 列和临时合并列
        drop_cols = [c for c in df.columns if c.endswith("_pre_it") or c in
                     ["first_it_date", "first_haic_date", "pre_it_source",
                      "pre_it_lab_date", "pre_it_days_before_it"]]
        df = df.drop(columns=drop_cols, errors="ignore")

        df = df.merge(pre_it_df, on="patient_id", how="left")

        # 计算派生指标
        df = compute_derived(df)

        # 计算变化率
        df = compute_change_rates(df)

        # 计算 pre-IT AFP-PIVKA 组合指标（ratio, sum, PC1, PC2）
        df = add_pre_it_composites(df)

        # 计算 pre-IT AFP-PIVKA 动态组合指标（change_diff / sum / PC1_dyn / PC2_dyn）
        df = add_pre_it_dynamic_composites(df)

        # 打印统计
        pre_it_cols = [c for c in df.columns if c.endswith("_pre_it") or c.endswith("_change_pre_it")]
        print(f"\n  New pre-IT columns ({len(pre_it_cols)}):")
        for c in sorted(pre_it_cols):
            if c in ["pre_it_source", "pre_it_lab_date", "pre_it_days_before_it", "first_it_date"]:
                continue
            n_valid = pd.to_numeric(df[c], errors="coerce").notna().sum()
            print(f"    {c}: n_valid={n_valid}/{n_total}")

        # 时间窗统计
        days_col = df["pre_it_days_before_it"].dropna()
        if len(days_col) > 0:
            print(f"\n  pre_it_days_before_it: median={days_col.median():.0f}, "
                  f"mean={days_col.mean():.1f}, min={days_col.min():.0f}, max={days_col.max():.0f}")

        # 保存
        df.to_csv(csv_path, index=False)
        print(f"  Saved: {csv_path} ({len(df)} rows, {len(df.columns)} cols)")

    print("\n\n=== Pre-IT extraction completed ===")
    print("Ready for 01_rcs_afp_pivka_composite.R")


if __name__ == "__main__":
    main()
