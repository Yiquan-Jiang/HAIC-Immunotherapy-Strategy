#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
数据准备脚本 — update_group_7
=============================
将 baseline_imputed.csv 与 patient_treatment_sequence_labels.csv 合并，
用 main_group 替换旧分组，排除 grey_zone 和 before_haic，
输出 analysis_ready.csv 供后续 PSM / KM / 森林图 / 表格使用。
"""

import pandas as pd
import os

BASE = "/Users/yqj/Nutstore Files/我的坚果云/Liver_tumor_big_data/FIRST_LINE_HAIC_2025-12-30/HAIC_NO_TACE_4_TIDY"
DATA_SRC = os.path.join(BASE, "data", "haic_no_tace_4")
DATA_DST = os.path.join(BASE, "update_group_7", "data")

EXCLUDE_GROUPS = {'grey_zone', 'before_haic'}

print("=" * 60)
print("  数据准备 — update_group_7")
print("=" * 60)

# 1. 读取
bl = pd.read_csv(os.path.join(DATA_SRC, "HAIC_NO_TACE_4_TIDY_baseline_imputed.csv"))
seq = pd.read_csv(os.path.join(DATA_SRC, "patient_treatment_sequence_labels.csv"))
print(f"\nbaseline_imputed: {len(bl)} 行")
print(f"sequence_labels:  {len(seq)} 行")

# 2. 合并 main_group
merged = bl.merge(seq[['patient_id', 'main_group']], on='patient_id', how='left')
print(f"\n合并后: {len(merged)} 行")
print(f"main_group 匹配成功: {merged['main_group'].notna().sum()}")
print(f"main_group 匹配失败: {merged['main_group'].isna().sum()}")

# 3. 排除 grey_zone 和 before_haic
before = len(merged)
df = merged[~merged['main_group'].isin(EXCLUDE_GROUPS) & merged['main_group'].notna()].copy()
print(f"\n排除 {EXCLUDE_GROUPS} 后: {len(df)} 行（排除 {before - len(df)} 例）")

# 4. 删除旧的 treatment_pattern 列（如果存在）
if 'treatment_pattern' in df.columns:
    df.drop(columns=['treatment_pattern'], inplace=True)

# 5. 分组分布
print(f"\n最终分组分布（7组）:")
print(df['main_group'].value_counts().to_string())
print(f"\n总计: {len(df)} 例")

# 6. 保存
out_file = os.path.join(DATA_DST, "analysis_ready.csv")
df.to_csv(out_file, index=False)
print(f"\n已保存: {out_file}")
print(f"  行数: {len(df)}, 列数: {len(df.columns)}")
