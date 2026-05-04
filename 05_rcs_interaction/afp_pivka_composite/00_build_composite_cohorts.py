#!/usr/bin/env python3
"""
AFP-PIVKA 组合指标构建脚本。

读取上层目录已有的 4 个 cohort CSV，为每个队列新增组合指标后保存到本目录。
不修改原始 cohort 文件。

组合指标体系
==========
静态（baseline & pre3）:
  1. log_afp_pivka_ratio  = log(AFP+1) - log(PIVKA+1)       → 相对平衡
  2. log_afp_pivka_sum    = log(AFP+1) + log(PIVKA+1)       → 双标志物总负荷
  3. afp_pivka_pc1 / pc2  = PCA 正交分解                    → 共同趋势 / 分离方向

动态（变化率）:
  4. afp_pivka_change_diff    = AFP_change - PIVKA_change   → 变化分离度
  5. afp_pivka_change_sum     = AFP_change + PIVKA_change   → 整体反应方向
  6. afp_pivka_pc1_dyn / pc2_dyn = PCA 正交分解（动态）

下游: 01_rcs_afp_pivka_composite.R
"""

import os
import warnings

warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PARENT_DIR = os.path.dirname(SCRIPT_DIR)

COHORT_FILES = {
    "THEN_IT": {
        "csv": "cohort_HAIC_alone_vs_HAIC_then_IT.csv",
        "trt_col": "trt_haic_then_it",
        "src_trt_col": "trt_compare",
    },
    "THEN_I": {
        "csv": "cohort_HAIC_alone_vs_HAIC_then_I.csv",
        "trt_col": "trt_haic_then_i",
        "src_trt_col": "trt_compare",
    },
    "THEN_T": {
        "csv": "cohort_HAIC_alone_vs_HAIC_then_T.csv",
        "trt_col": "trt_haic_then_t",
        "src_trt_col": "trt_compare",
    },
    "T_CONC": {
        "csv": "cohort_HAIC_alone_vs_HAIC_T_conc.csv",
        "trt_col": "trt_haic_t_conc",
        "src_trt_col": "trt_compare",
    },
    "I_CONC": {
        "csv": "cohort_HAIC_alone_vs_HAIC_I_conc.csv",
        "trt_col": "trt_haic_i_conc",
        "src_trt_col": "trt_compare",
    },
    "IT_CONC": {
        "csv": "cohort_HAIC_alone_vs_HAIC_IT_conc.csv",
        "trt_col": "trt_haic_it_conc",
        "src_trt_col": "trt_compare",
    },
}


def add_static_composites(df: pd.DataFrame, suffix: str) -> pd.DataFrame:
    """为 baseline 或 pre3 的 AFP/PIVKA 添加组合指标。

    suffix: '_bl' 或 '_pre3'
    """
    afp_col = "afp" if suffix == "_bl" else "afp_pre3"
    pivka_col = "pivka" if suffix == "_bl" else "pivka_pre3"

    afp_vals = pd.to_numeric(df[afp_col], errors="coerce")
    pivka_vals = pd.to_numeric(df[pivka_col], errors="coerce")

    log_afp = np.log1p(afp_vals.clip(lower=0))
    log_pivka = np.log1p(pivka_vals.clip(lower=0))

    df[f"log_afp_pivka_ratio{suffix}"] = log_afp - log_pivka
    df[f"log_afp_pivka_sum{suffix}"] = log_afp + log_pivka

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
        print(f"  PCA{suffix}: loadings={np.round(loadings, 3).tolist()}, "
              f"var_explained={np.round(explained, 3).tolist()}")

    df[f"afp_pivka_pc1{suffix}"] = pc1
    df[f"afp_pivka_pc2{suffix}"] = pc2

    return df


def add_dynamic_composites(df: pd.DataFrame) -> pd.DataFrame:
    """为变化率（AFP_change_pre3, PIVKA_change_pre3）添加组合指标。"""
    afp_chg = pd.to_numeric(df["afp_change_pre3"], errors="coerce")
    pivka_chg = pd.to_numeric(df["pivka_change_pre3"], errors="coerce")

    df["afp_pivka_change_diff"] = afp_chg - pivka_chg
    df["afp_pivka_change_sum"] = afp_chg + pivka_chg

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
        print(f"  PCA_dyn: loadings={np.round(loadings, 3).tolist()}, "
              f"var_explained={np.round(explained, 3).tolist()}")

    df["afp_pivka_pc1_dyn"] = pc1
    df["afp_pivka_pc2_dyn"] = pc2

    return df


def process_one_cohort(key: str, info: dict):
    src = os.path.join(PARENT_DIR, info["csv"])
    if not os.path.exists(src):
        print(f"  SKIP: {src} not found")
        return

    df = pd.read_csv(src)
    # 重命名 trt 列（源文件统一为 trt_compare，下游 R 脚本需要各自的列名）
    src_trt = info.get("src_trt_col", "")
    dst_trt = info["trt_col"]
    if src_trt and src_trt in df.columns and src_trt != dst_trt:
        df = df.rename(columns={src_trt: dst_trt})
    print(f"\n{'='*60}")
    print(f"[{key}] {info['csv']}  n={len(df)}  trt_col={dst_trt}")
    print(f"{'='*60}")

    df = add_static_composites(df, "_bl")
    df = add_static_composites(df, "_pre3")
    df = add_dynamic_composites(df)

    new_cols = [c for c in df.columns if any(
        c.startswith(p) for p in [
            "log_afp_pivka_", "afp_pivka_pc", "afp_pivka_change_",
        ]
    )]
    print(f"\n  New composite columns ({len(new_cols)}):")
    for c in new_cols:
        n_valid = df[c].notna().sum()
        stats = df[c].describe()
        print(f"    {c}: n_valid={n_valid}, "
              f"mean={stats['mean']:.3f}, std={stats['std']:.3f}, "
              f"min={stats['min']:.3f}, max={stats['max']:.3f}")

    out_path = os.path.join(SCRIPT_DIR, f"composite_{key}_cohort.csv")
    df.to_csv(out_path, index=False)
    print(f"  Saved: {out_path}")

    concordance_patterns(df)


def concordance_patterns(df: pd.DataFrame):
    afp_chg = df["afp_change_pre3"]
    pivka_chg = df["pivka_change_pre3"]
    mask = afp_chg.notna() & pivka_chg.notna()
    n = mask.sum()
    if n == 0:
        return
    afp_down = afp_chg[mask] < 0
    pivka_down = pivka_chg[mask] < 0
    print(f"\n  Concordance/Discordance (n={n}):")
    print(f"    Both down:      {(afp_down & pivka_down).sum()} "
          f"({(afp_down & pivka_down).mean()*100:.1f}%)")
    print(f"    Both up:        {(~afp_down & ~pivka_down).sum()} "
          f"({(~afp_down & ~pivka_down).mean()*100:.1f}%)")
    print(f"    AFP down/PIVKA up:  {(afp_down & ~pivka_down).sum()} "
          f"({(afp_down & ~pivka_down).mean()*100:.1f}%)")
    print(f"    AFP up/PIVKA down:  {(~afp_down & pivka_down).sum()} "
          f"({(~afp_down & pivka_down).mean()*100:.1f}%)")


if __name__ == "__main__":
    print("=" * 60)
    print("AFP-PIVKA Composite Indicator Builder")
    print("=" * 60)

    for key, info in COHORT_FILES.items():
        process_one_cohort(key, info)

    print("\n\nAll cohorts processed. Ready for 01_rcs_afp_pivka_composite.R")
