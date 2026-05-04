# Interaction Analysis - PSM Matched Cohort Forest Plots

## 项目概述

本项目分析 HAIC alone vs HAIC then I (或 HAIC then I+T) 的治疗策略比较，使用PSM匹配后的队列进行森林图绘制和亚组交互分析。

## 目录结构

```
interaction_analysis/
├── README.md                   # 本说明文件
├── 01_publication_figures.py   # 主分析脚本 - PSM02 (HAIC alone vs HAIC then I)
├── 02_publication_figures_ids06_IplusT.py # PSM06分析脚本 (HAIC alone vs HAIC then I+T)
├── data/                       # 数据目录
│   ├── HAIC_NO_TACE_4_TIDY_baseline.csv      # 基线数据
│   ├── HAIC_NO_TACE_4_TIDY_longitudinal.csv  # 纵向数据
│   ├── matched_ids_06_HAIC_alone_vs_HAIC_then_I+T.csv # PSM06匹配ID
│   ├── 00_swimmer_plot_events.csv            # Swimmer plot事件数据
│   └── trace_planB_relaxed.csv               # Trace plan B数据
└── output/                     # 输出目录
    ├── Fig4_forest_v4.pdf/png                 # Figure 4森林图
    ├── psm02_tidymatched_unweighted_*         # PSM02分析结果
    ├── planB_decision_haic3_relaxed_*         # Plan B决策分析
    ├── planB_psm02_matched_unweighted_*       # Plan B PSM02结果
    └── psm06_IplusT/                          # PSM06 I+T分析结果
        ├── psm06_tidymatched_unweighted_IplusT_Fig4_forest_v4.*
        ├── psm06_tidymatched_unweighted_IplusT_SuppFig_lm_forest.*
        └── psm06_tidymatched_unweighted_IplusT_Table_lm_forest.csv
```

## 数据依赖

- **基线数据**: `HAIC_NO_TACE_4_TIDY_baseline.csv` - 包含患者基线特征
- **纵向数据**: `HAIC_NO_TACE_4_TIDY_longitudinal.csv` - 包含时间序列生物标志物数据
- **匹配ID**: 来自 `HAIC_NO_TACE_4_TIDY/update_group_7/results/psm_balance_tables_complete/`

## 运行步骤

### 1. 确保数据文件存在

```bash
# 检查数据文件
ls data/
```

### 2. 运行PSM02分析 (HAIC alone vs HAIC then I)

```bash
cd HAIC_NO_TACE_4_TIDY/update_group_7/scripts/interaction_analysis
python 01_publication_figures.py
```

### 3. 运行PSM06分析 (HAIC alone vs HAIC then I+T)

```bash
python 02_publication_figures_ids06_IplusT.py
```

## Python依赖

```python
# 必需包
pandas
numpy
matplotlib
lifelines  # KaplanMeierFitter, CoxPHFitter
scipy
```

安装依赖:
```bash
pip install pandas numpy matplotlib lifelines scipy
```

## 输出说明

### Figure 4 森林图
- 展示各亚组的HR和95% CI
- 包含P_interaction标注
- PSM02 vs PSM06 分别输出

### Landmark敏感性分析
- Landmark时间: 42天 (约1.38个月)
- τ = 16个月 (从landmark起)

### NPG配色标准
- 对照组(HAIC alone): `#4DBBD5` 青色
- 治疗组(HAIC then I/I+T): `#E64B35` 红橙
- P_interaction显著: `#FF1744` 鲜红 (P<0.05)

## 分析队列

| 对比组 | PSM编号 | 说明 |
|--------|---------|------|
| HAIC alone vs HAIC then I | PSM02 | 序贯免疫治疗 |
| HAIC alone vs HAIC then I+T | PSM06 | 序贯免疫+靶向治疗 |

## 注意事项

1. 脚本使用相对路径，需在脚本所在目录运行
2. 输出PDF+PNG双格式（300 DPI）
3. 所有分析使用非加权方法（unweighted）
4. 遵循Nature期刊配色标准