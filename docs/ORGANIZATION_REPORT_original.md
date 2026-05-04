# 研究项目整理报告

## 整理时间
2026-04-06

## 项目概述

两个研究项目已整理到 `HAIC_NO_TACE_4_TIDY/update_group_7/scripts/` 目录下，确保可独立运行。

---

## 1. interaction_analysis

**路径**: `HAIC_NO_TACE_4_TIDY/update_group_7/scripts/interaction_analysis/`

### 功能
PSM匹配队列的森林图分析，比较 HAIC alone vs HAIC then I/I+T 的治疗效果。

### 目录结构
```
interaction_analysis/
├── README.md                    # 项目说明
├── run_analysis.sh              # 运行脚本 (可执行)
├── 01_publication_figures.py    # PSM02分析 (HAIC alone vs HAIC then I)
├── 02_publication_figures_ids06_IplusT.py # PSM06分析 (HAIC alone vs HAIC then I+T)
├── data/                        # 数据文件 (6个)
│   ├── HAIC_NO_TACE_4_TIDY_baseline.csv (1.8MB)
│   ├── HAIC_NO_TACE_4_TIDY_longitudinal.csv (6.1MB)
│   ├── matched_ids_02_HAIC_alone_vs_HAIC_then_I.csv (30KB)
│   ├── matched_ids_06_HAIC_alone_vs_HAIC_then_I+T.csv (46KB)
│   ├── 00_swimmer_plot_events.csv (9.0MB)
│   └── trace_planB_relaxed.csv (485KB)
└── output/                      # 输出文件
    ├── Fig4_forest_v4.pdf/png
    ├── psm02_tidymatched_unweighted_* (森林图+表格)
    ├── planB_* (Plan B分析结果)
    └── psm06_IplusT/ (PSM06 I+T分析)
```

### 运行方式
```bash
cd HAIC_NO_TACE_4_TIDY/update_group_7/scripts/interaction_analysis
./run_analysis.sh
```

---

## 2. rcs_non_linear_int_investigation

**路径**: `HAIC_NO_TACE_4_TIDY/update_group_7/scripts/rcs_non_linear_int_investigation/`

### 功能
RCS非线性交互分析，探索治疗效应与生物标志物的非线性关系。

### 目录结构
```
rcs_non_linear_int_investigation/
├── README.md                    # 项目说明
├── run_analysis.sh              # 运行脚本 (可执行)
├── RCS_INT_HAIC_ALONE_AND_THEN_I_build_cohort.py  # THEN_I队列构建
├── RCS_INT_HAIC_ALONE_AND_THEN_I_dual_timescale.R # THEN_I RCS分析
├── RCS_INT_HAIC_ALONE_AND_THEN_IT_build_cohort.py # THEN_I+T队列构建
├── RCS_INT_HAIC_ALONE_AND_THEN_IT_dual_timescale.R # THEN_I+T RCS分析
├── data/                        # 数据文件 (5个)
│   ├── HAIC_NO_TACE_4_TIDY_baseline.csv (1.8MB)
│   ├── HAIC_NO_TACE_4_TIDY_baseline_imputed.csv (1.6MB)
│   ├── HAIC_NO_TACE_4_TIDY_longitudinal.csv (6.1MB)
│   ├── RCS_INT_HAIC_ALONE_AND_THEN_I_cohort.csv (523KB)
│   └── RCS_INT_HAIC_ALONE_AND_THEN_IT_cohort.csv (550KB)
├── output/                      # 主分析输出
│   ├── RCS_INT_HAIC_ALONE_AND_THEN_I/
│   │   ├── landmark/ (动态/静态RCS图+ANOVA汇总)
│   │   └── total_os/
│   └── RCS_INT_HAIC_ALONE_AND_THEN_IT/
│       ├── landmark/
│       └── total_os/
└── afp_pivka_composite/         # AFP-PIVKA组合指标分析
    ├── 00_build_composite_cohorts.py
    ├── 01_rcs_afp_pivka_composite.R
    ├── plot_afp_negative_trajectory.R
    ├── plot_baseline_distribution.R
    ├── composite_*_cohort.csv (4个队列文件)
    └── output/
        ├── afp_negative_dynamic_trajectory.pdf/png
        ├── afp_pivka_baseline_distribution.pdf/png
        ├── THEN_I/landmark/total_os/
        ├── THEN_IT/landmark/total_os/
        ├── THEN_T/landmark/total_os/
        └── T_CONC/landmark/total_os/
```

### 运行方式
```bash
cd HAIC_NO_TACE_4_TIDY/update_group_7/scripts/rcs_non_linear_int_investigation
./run_analysis.sh
```

---

## 脚本修改说明

### interaction_analysis
- 修改路径配置为相对路径，使用 `data/` 子目录
- 添加 `SWIMMER_EVENTS_PATH` 全局变量
- 所有数据文件自包含，无需外部依赖

### rcs_non_linear_int_investigation
- R脚本已使用环境变量或相对路径
- Python脚本使用 `SCRIPT_DIR` 相对路径
- 数据文件自包含

---

## 依赖要求

### Python
```bash
pip install pandas numpy matplotlib lifelines scipy
```

### R
```r
install.packages(c("survival", "rms", "ggplot2", "dplyr", "gridExtra", "glmnet"))
```

---

## 数据完整性检查

| 项目 | 数据文件 | 状态 |
|------|----------|------|
| interaction_analysis | 6个数据文件 | ✓ 完整 |
| rcs_non_linear_int_investigation | 5个主数据 + 4个队列 | ✓ 完整 |
| afp_pivka_composite | 4个组合队列 | ✓ 完整 |

---

## 输出文件统计

### interaction_analysis
- 森林图: 10个PDF + 10个PNG
- 数据表: 5个CSV

### rcs_non_linear_int_investigation
- RCS曲线图: 24个PDF + 24个PNG (主分析)
- AFP-PIVKA图: 10个PDF + 10个PNG
- ANOVA汇总: 8个CSV
- HR交叉点: 4个CSV