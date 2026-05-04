# RCS Non-linear Interaction Investigation

## 项目概述

本项目使用RMS库的受限立方样条（Restricted Cubic Splines, RCS）方法，探索治疗效应与生物标志物之间的非线性交互关系。采用双时间尺度（landmark + total OS）分析策略。

## 分析队列（6组）

| 队列ID | 对照组 | 治疗组 | 说明 |
|--------|--------|--------|------|
| THEN_I | HAIC_alone | HAIC_then_I | HAIC后序贯免疫治疗 |
| THEN_IT | HAIC_alone | HAIC_then_I+T | HAIC后序贯免疫+靶向治疗 |
| THEN_T | HAIC_alone | HAIC_then_T | HAIC后序贯靶向治疗 |
| I_CONC | HAIC_alone | HAIC+I_concurrent | HAIC与免疫同时进行 |
| T_CONC | HAIC_alone | HAIC+T_concurrent | HAIC与靶向同时进行 |
| I_T_CONC | HAIC_alone | HAIC+I+T_concurrent | HAIC与免疫+靶向同时进行 |

## 目录结构

```
rcs_non_linear_int_investigation/
├── README.md                                        # 本说明文件
├── run_analysis.sh                                  # 一键运行脚本
│
├── RCS_INT_HAIC_ALONE_AND_THEN_I_*.py/R             # THEN_I队列
├── RCS_INT_HAIC_ALONE_AND_THEN_IT_*.py/R            # THEN_I+T队列
├── RCS_INT_HAIC_ALONE_AND_THEN_T_*.py/R             # THEN_T队列
├── RCS_INT_HAIC_ALONE_AND_I_CONC_*.py/R             # I_CONC队列
├── RCS_INT_HAIC_ALONE_AND_T_CONC_*.py/R             # T_CONC队列
├── RCS_INT_HAIC_ALONE_AND_I_T_CONC_*.py/R           # I_T_CONC队列
│
├── data/                                            # 数据目录
│   ├── HAIC_NO_TACE_4_TIDY_baseline_imputed.csv     # 插补后基线数据
│   ├── HAIC_NO_TACE_4_TIDY_longitudinal.csv         # 纵向数据
│   └── RCS_INT_*_cohort.csv                         # 各队列分析数据
│
├── output/                                          # 主分析输出
│   ├── RCS_INT_HAIC_ALONE_AND_THEN_I/
│   │   ├── landmark/                                # 42天landmark分析
│   │   └── total_os/                                # 全生存期分析
│   ├── RCS_INT_HAIC_ALONE_AND_THEN_IT/
│   ├── RCS_INT_HAIC_ALONE_AND_THEN_T/
│   ├── RCS_INT_HAIC_ALONE_AND_I_CONC/
│   ├── RCS_INT_HAIC_ALONE_AND_T_CONC/
│   └── RCS_INT_HAIC_ALONE_AND_I_T_CONC/
│
└── afp_pivka_composite/                             # AFP-PIVKA组合指标分析
    ├── 00_build_composite_cohorts.py                # 组合队列构建
    ├── 01_rcs_afp_pivka_composite.R                  # 组合指标RCS分析
    └── output/
        ├── THEN_I/landmark/total_os/
        ├── THEN_IT/landmark/total_os/
        ├── THEN_T/landmark/total_os/
        ├── T_CONC/landmark/total_os/
        ├── I_CONC/landmark/total_os/
        └── I_T_CONC/landmark/total_os/
```

## 分析方法

### RCS模型

模型形式: `Surv(time, death_status) ~ trt * rcs(rcsx, nk=3)`

- **暴露项**: 治疗组 × RCS(生物标志物)
- **nk=3**: 3个节点（默认）
- **IPTW加权**: 使用glmnet ridge倾向评分 + 稳定权重

### 双时间尺度

| 时间尺度 | 说明 |
|----------|------|
| landmark | 42天landmark后生存分析 |
| total_os | 全生存期OS分析 |

## 运行步骤

### 方式一：一键运行（推荐）

```bash
cd HAIC_NO_TACE_4_TIDY/update_group_7/scripts/rcs_non_linear_int_investigation
./run_analysis.sh
```

### 方式二：单独运行

#### 1. 运行THEN_I分析

```bash
# Step 1: 构建队列
python RCS_INT_HAIC_ALONE_AND_THEN_I_build_cohort.py

# Step 2: RCS分析
Rscript RCS_INT_HAIC_ALONE_AND_THEN_I_dual_timescale.R
```

#### 2. 运行THEN_I+T分析

```bash
python RCS_INT_HAIC_ALONE_AND_THEN_IT_build_cohort.py
Rscript RCS_INT_HAIC_ALONE_AND_THEN_IT_dual_timescale.R
```

#### 3. 运行THEN_T分析

```bash
python RCS_INT_HAIC_ALONE_AND_THEN_T_build_cohort.py
Rscript RCS_INT_HAIC_ALONE_AND_THEN_T_dual_timescale.R
```

#### 4. 运行I_CONC分析

```bash
python RCS_INT_HAIC_ALONE_AND_I_CONC_build_cohort.py
Rscript RCS_INT_HAIC_ALONE_AND_I_CONC_dual_timescale.R
```

#### 5. 运行T_CONC分析

```bash
python RCS_INT_HAIC_ALONE_AND_T_CONC_build_cohort.py
Rscript RCS_INT_HAIC_ALONE_AND_T_CONC_dual_timescale.R
```

#### 6. 运行I_T_CONC分析

```bash
python RCS_INT_HAIC_ALONE_AND_I_T_CONC_build_cohort.py
Rscript RCS_INT_HAIC_ALONE_AND_I_T_CONC_dual_timescale.R
```

#### 7. AFP-PIVKA组合分析

```bash
cd afp_pivka_composite

# 构建组合队列
python 00_build_composite_cohorts.py

# RCS分析 (可选择队列)
Rscript 01_rcs_afp_pivka_composite.R THEN_I
Rscript 01_rcs_afp_pivka_composite.R THEN_IT
Rscript 01_rcs_afp_pivka_composite.R THEN_T
Rscript 01_rcs_afp_pivka_composite.R T_CONC
Rscript 01_rcs_afp_pivka_composite.R I_CONC
Rscript 01_rcs_afp_pivka_composite.R I_T_CONC
Rscript 01_rcs_afp_pivka_composite.R ALL  # 全部运行
```

## Python依赖

```python
pandas
numpy
scikit-learn  # 用于PCA（组合指标）
```

## R依赖

```r
library(survival)
library(rms)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(glmnet)
library(grid)
```

安装R依赖:
```r
install.packages(c("survival", "rms", "ggplot2", "dplyr", "gridExtra", "glmnet"))
```

## 输出说明

### RCS曲线图

- **dynamic**: 动态RCS（使用变化率）
- **static_baseline**: 静态RCS（使用基线值）
- **static_pre3**: 静态RCS（使用pre-HAIC-3值）

### ANOVA汇总表

- P值检验非线性交互效应
- HR=0.7, HR=0.8交叉点定位

### 组合指标分析

- AFP-PIVKA组合变量
- 6种治疗策略对比
- 基线分布 + 动态轨迹

## 环境变量（可选）

```bash
# 自定义参数
export RMS_RCS_NK=3         # RCS节点数
export RMS_RCS_N_BOOT=200   # Bootstrap次数
export RMS_RCS_DATA_CSV=/path/to/data.csv
export RMS_RCS_OUT_DIR=/path/to/output
export AFP_PIVKA_COHORT=THEN_I  # 指定队列
```

## 注意事项

1. Python脚本构建队列，R脚本执行分析
2. 队列CSV需先运行Python脚本生成
3. IPTW权重排除静态交互变量
4. 所有分析输出PDF+PNG双格式
5. Bootstrap默认200次（可通过环境变量调整）
6. 每个队列对应PSM匹配名单在 `../../../results/psm_balance_tables_complete/`