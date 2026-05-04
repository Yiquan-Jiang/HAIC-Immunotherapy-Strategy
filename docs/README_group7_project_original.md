# update_group_7 — HAIC 7组治疗方案对比分析

## 研究概述

基于"前4次 HAIC 未联合 TACE"入组标准的患者队列（N=4,233），将原 5 分组方案更新为 **7 组精细化治疗分组**，排除 `grey_zone` 和 `before_haic` 后最终纳入 **3,885 例**患者，完成全部 **21 组两两 PSM 匹配**及生存分析。

## 7 个治疗组

| 组别 | 英文缩写 | 样本量 | 说明 |
|---|---|---|---|
| HAIC 单药 | `HAIC_alone` | 1,119 | 仅接受 HAIC 治疗 |
| HAIC + 免疫（同步） | `HAIC+I_concurrent` | 274 | HAIC 同步联合免疫治疗 |
| HAIC → 免疫（序贯） | `HAIC_then_I` | 152 | HAIC 后序贯免疫治疗 |
| HAIC + 靶向（同步） | `HAIC+T_concurrent` | 372 | HAIC 同步联合靶向治疗 |
| HAIC → 靶向（序贯） | `HAIC_then_T` | 247 | HAIC 后序贯靶向治疗 |
| HAIC + 免疫 + 靶向（同步） | `HAIC+I+T_concurrent` | 1,500 | HAIC 同步联合免疫+靶向 |
| HAIC → 免疫 + 靶向（序贯） | `HAIC_then_I+T` | 221 | HAIC 后序贯免疫+靶向 |

> 排除组：`grey_zone`（178例）、`before_haic`（170例）

## 目录结构

```
update_group_7/
├── README.md                          ← 本文件
├── scripts/
│   ├── run_all.sh                     一键复现脚本（含 Step7）
│   ├── step0_prepare_data.py          数据准备（合并分组、排除）
│   ├── step3_psm_analysis.R           PSM 匹配 + 生存分析（21组）
│   ├── step4_km_curves.py             KM 生存曲线（Python，发表级）
│   ├── step5_forest_plot.py           森林图（HR 对比，PSM 前后）
│   ├── step6_tables_and_loveplots.R   Table 1 + 平衡表 + Love Plot
│   ├── step7_subgroup_ow.R            ★ 高危亚组 OW 分析（主分析）
│   ├── step7_subgroup_plots.py        ★ 亚组森林图 + KM 曲线（Python）
│   ├── step7_ow_balance_table.R       ★ OW 加权后基线平衡表（R）
│   ├── step7_ow_balance_figure.py     ★ 基线平衡表可视化（Python）
│   └── swimmer_plot_7groups.R       7 组 Swimmer（Fig2d 模板：geom_rect + inset 图例 + patchwork）
├── data/
│   ├── analysis_ready.csv             ★ 分析用数据（3,885例 × 67列）
│   ├── HAIC_NO_TACE_4_TIDY_baseline_imputed.csv  插补后基线数据
│   ├── patient_treatment_sequence_labels.csv       新分组标签来源
│   └── ...（其他原始数据副本）
├── results/
│   ├── psm_balance_tables_complete/
│   │   ├── survival_analysis_final.csv             21组 HR/CI/P 汇总
│   │   └── matched_ids_XX_*.csv                    21组匹配患者ID
│   ├── tables/
│   │   ├── table1_overall_baseline.docx            Table 1（7组基线）
│   │   └── tableXX_compXX_*_psm_balance.docx       21组 PSM 平衡表
│   └── subgroup_analysis/                          ★ Step7 输出目录
│       ├── ow_subgroup_results.csv                 亚组 OW HR/CI/P/等价性
│       ├── ow_smd_balance.csv                      各亚组 SMD（OW前后）
│       ├── ow_interaction_tests.csv                交互作用检验 P 值
│       └── ow_weighted_ids_*.csv                   各亚组加权患者数据（供 KM）
├── figures/
│   ├── km/                            7组整体 + 21组两两 KM 曲线
│   ├── psm_pub_quality/               森林图
│   ├── loveplots/                     21组 Love Plot
│   ├── swimmer_7groups/               7组分别：swimmer_<group>.pdf/png
│   └── subgroup/                      ★ Step7 输出目录
│       ├── ow_love_*.pdf/png          各亚组 Love Plot
│       ├── ow_subgroup_forest_plot.pdf/png   亚组森林图
│       ├── ow_subgroup_km_curves.pdf/png     亚组 KM 曲线（4面板）
│       └── ow_balance_table_*.pdf/png        OW 加权后基线平衡表
└── logs/
    ├── psm_analysis.log
    └── tables_loveplots.log
```

## 快速复现

### 环境要求

**Python 3.8+**
```bash
pip install pandas numpy lifelines matplotlib scikit-learn
```

**R 4.0+**
```r
install.packages(c("tidyverse", "MatchIt", "survival", "survminer",
                    "cobalt", "gtsummary", "flextable", "officer",
                    "WeightIt", "survey"))
```

### 一键运行

```bash
cd update_group_7/scripts

# 运行全部步骤（约 25-35 分钟，含 Step7）
bash run_all.sh

# 从 Step3 开始（跳过数据准备）
bash run_all.sh 3

# 只运行亚组分析（Step7，约 5-8 分钟）
bash run_all.sh 7

# 只运行特定步骤
bash run_all.sh 4 5    # 只跑 KM 曲线和森林图
bash run_all.sh 7a 7b 7c 7d  # 只跑 Step7 各子步骤
```

### 分步运行

```bash
# Step 0: 数据准备（~10秒）
python3 step0_prepare_data.py

# Step 3: PSM 匹配 + 生存分析（~5-10分钟）
Rscript step3_psm_analysis.R

# Step 4: KM 生存曲线（~2-3分钟）
python3 step4_km_curves.py

# Step 5: 森林图（~1分钟）
python3 step5_forest_plot.py

# Step 6: Table 1 + 平衡表 + Love Plot（~10-15分钟）
Rscript step6_tables_and_loveplots.R

# ── Step 7: 高危亚组分析（OW）──────────────────────────────
# Step 7a: OW 分析（R）— 生成 HR/CI/P/等价性 + Love Plot（~2分钟）
Rscript step7_subgroup_ow.R

# Step 7b: 亚组森林图 + KM 曲线（Python）（~1分钟）
python3 step7_subgroup_plots.py

# Step 7c: OW 加权后基线平衡表（R）（~2分钟）
Rscript step7_ow_balance_table.R

# Step 7d: 基线平衡表可视化（Python）（~1分钟）
python3 step7_ow_balance_figure.py
```

## 分析流程说明

### Step 0: 数据准备 (`step0_prepare_data.py`)

- **输入**: `HAIC_NO_TACE_4_TIDY_baseline_imputed.csv`（4,233例）+ `patient_treatment_sequence_labels.csv`
- **操作**: 合并 `main_group` 列 → 排除 `grey_zone`/`before_haic` → 删除旧 `treatment_pattern` 列
- **输出**: `analysis_ready.csv`（3,885例 × 67列）

### Step 3: PSM 匹配 (`step3_psm_analysis.R`)

- **方法**: 1:1 最近邻匹配（logistic regression 倾向评分）
- **匹配变量**（20个）:
  - 肿瘤标志物: AFP 分类、PIVKA-II 分类、PIVKA-II 标准化值
  - 肿瘤特征: 直径>10cm、多发、PVTT 分级、PVTT/HVTT/IVC-RA 有无
  - 转移: 远处转移、淋巴结转移
  - 肝功能: ALBI 分级、胆红素标准化、白蛋白标准化、血小板标准化
  - 其他: 腹水、静脉曲张、年龄标准化、肿瘤直径标准化、NLR 标准化
- **自适应卡钳**: 小组(<150) → 0.25; 中组(150-300) → 0.15; 大组(>300) → 0.10
- **输出**: 21组匹配患者ID + 生存分析汇总（HR、95%CI、P值、中位OS）

### Step 4: KM 生存曲线 (`step4_km_curves.py`)

- 7组整体 KM 曲线（PSM 前）
- 21组两两 KM 曲线（左：PSM 前，右：PSM 后）
- 使用 `lifelines` 库，配合 R 的匹配结果
- 输出 PDF（矢量）+ PNG（300 DPI）

### Step 5: 森林图 (`step5_forest_plot.py`)

- 21组两两对比的 HR 汇总森林图
- 每组显示两行：灰色=PSM前，蓝色=PSM后，橙红=P<0.05
- 对数坐标轴，菱形标记 HR 点估计

### Step 5b: 参照组森林图 (`step5b_forest_vs_IT_concurrent.py`)

- **仅 PSM 后**：其余 6 组各自与 `HAIC+I+T_concurrent` 对比，**同一幅森林图**（HR = 行内组 / 参照组）
- R 中 Cox 参照为 `Group1`，故当 `Group2` 为参照时需对 HR 与 CI **取倒数**
- 输出：`figures/psm_pub_quality/HR_forest_vs_IT_concurrent_psm_after.pdf/png` 与 `forest_vs_IT_concurrent_psm_after.csv`

### Step 6: 表格与 Love Plot (`step6_tables_and_loveplots.R`)

- Table 1: 7组整体基线特征（Kruskal-Wallis / Chi-squared 检验）
- 21组 PSM 前后平衡表（Word 格式，`gtsummary` 生成）
- 21组 Love Plot（标准化均值差，`cobalt` 生成）

---

## Step 7: 高危亚组分析（OW）

> **研究问题**：在高危 HCC 亚组中，`HAIC+I_concurrent`（同步免疫）vs `HAIC_then_I`（延迟免疫）的 OS 是否存在显著差异？
>
> **方法**：Overlap Weighting (OW) + 加权 Cox 回归（robust sandwich SE）
>
> **等价性判定**：90% CI 落入 [0.60, 1.67] 即认为等价

### 亚组定义

| 亚组 | 定义 | 样本量（估计） |
|---|---|---|
| Composite high-risk | IVC/RA 或 肿瘤数>3 或 PVTT Vp3/4 或 远处转移 或 直径>10cm | ~200 |
| Tumor count >3 | 肿瘤数目 >3 | ~150 |
| Tumor diameter >10 cm | 最大肿瘤直径 >10 cm | ~120 |
| PVTT Vp3/4 | 门静脉癌栓 Vp3 或 Vp4 | ~130 |
| Extrahepatic metastasis | 远处转移阳性 | ~102 |

互补亚组（Non-*）同步分析，用于交互作用检验。

### 方法学说明

**为何选择 OW 而非 PSM？**
- 亚组内样本量有限（最小 n≈102），PSM 1:1 后每组仅约 30-40 例，等价性检验效能严重不足
- OW 权重有界（最大值 = PS×(1-PS) ≤ 0.25），小样本下天然稳定，无极端权重风险
- OW 推断目标（ATO）聚焦于"两种治疗均可选"的重叠人群，与本研究临床问题高度匹配
- 参考文献：Li, Morgan & Zaslavsky (2018) *JASA*; Mao et al. (2019) *Statistics in Medicine*

### Step 7a: OW 分析 (`step7_subgroup_ow.R`)

- **输入**: `data/analysis_ready.csv`
- **PS 模型**: logistic regression，15 个核心协变量（亚组定义变量自动排除）
- **权重**: `w_treat = 1 - PS`，`w_ctrl = PS`（estimand = "ATO"）
- **效应估计**: `survey::svycoxph`，robust sandwich SE
- **平衡性**: `cobalt::bal.tab`，SMD 阈值 < 0.10
- **等价性**: 90% CI = HR ± 1.645 × robust SE
- **输出**:
  - `results/subgroup_analysis/ow_subgroup_results.csv`：HR/CI/P/等价性/ESS/SMD
  - `results/subgroup_analysis/ow_smd_balance.csv`：各变量 SMD（OW前后）
  - `results/subgroup_analysis/ow_interaction_tests.csv`：交互作用 P 值
  - `results/subgroup_analysis/ow_weighted_ids_*.csv`：加权患者数据（供 Python KM）
  - `figures/subgroup/ow_love_*.pdf/png`：各亚组 Love Plot

### Step 7b: 亚组可视化 (`step7_subgroup_plots.py`)

- **输入**: `ow_subgroup_results.csv`，`ow_interaction_tests.csv`，`ow_weighted_ids_*.csv`
- **输出**:
  - `figures/subgroup/ow_subgroup_forest_plot.pdf/png`：亚组森林图（含等价区间 [0.60, 1.67]）
  - `figures/subgroup/ow_subgroup_km_curves.pdf/png`：4面板 KM 曲线（主要高危亚组）

### Step 7c: OW 加权后基线平衡表 (`step7_ow_balance_table.R`)

- **输入**: `data/analysis_ready.csv`
- **输出**:
  - `results/subgroup_analysis/ow_balance_table_full.csv`：长格式（变量×亚组×统计量）
  - `results/subgroup_analysis/ow_balance_table_wide.csv`：宽格式（供可视化）

### Step 7d: 基线平衡表可视化 (`step7_ow_balance_figure.py`)

- **输入**: `ow_balance_table_wide.csv`
- **输出**:
  - `figures/subgroup/ow_balance_table_*.pdf/png`：各亚组独立平衡表图
  - `figures/subgroup/ow_balance_table_combined.pdf/png`：合并平衡表（所有亚组）

## 关键结果摘要（PSM 后）

| # | 对比 | N (匹配后) | HR | 95% CI | P |
|---|---|---|---|---|---|
| 01 | HAIC alone vs HAIC+I concurrent | 271 vs 271 | 0.62 | 0.48–0.79 | **<0.001** |
| 02 | HAIC alone vs HAIC→I | 150 vs 150 | 0.72 | 0.53–0.98 | **0.036** |
| 03 | HAIC alone vs HAIC+T concurrent | 364 vs 364 | 0.75 | 0.61–0.92 | **0.005** |
| 04 | HAIC alone vs HAIC→T | 247 vs 247 | 1.03 | 0.82–1.28 | 0.802 |
| 05 | HAIC alone vs HAIC+I+T concurrent | 958 vs 958 | 0.74 | 0.64–0.84 | **<0.001** |
| 06 | HAIC alone vs HAIC→I+T | 219 vs 219 | 0.84 | 0.64–1.10 | 0.197 |
| 12 | HAIC→I vs HAIC+T concurrent | 142 vs 142 | 1.54 | 1.11–2.14 | **0.010** |
| 13 | HAIC→I vs HAIC→T | 134 vs 134 | 1.44 | 1.06–1.96 | **0.020** |

> 完整 21 组结果见 `results/psm_balance_tables_complete/survival_analysis_final.csv`

## 与旧 5 分组的关系

新 7 分组是对旧 5 分组的精细化拆分：
- 旧 `HAIC+I_early` → 新 `HAIC+I_concurrent`（+部分归入 `HAIC+I+T_concurrent`）
- 旧 `HAIC_then_I` → 拆分为新 `HAIC_then_I` + `HAIC+T_concurrent` + `HAIC_then_T`
- 旧 `Other` → 大量拆分到 `HAIC+T_concurrent`、`HAIC_then_T` 等具体组
- 新增 `HAIC+T_concurrent` 和 `HAIC_then_T` 两个靶向相关组

## 注意事项

1. **路径**: 脚本中使用绝对路径，如需迁移请全局替换 `BASE_DIR`
2. **R 包版本**: `gtsummary` ≥ 2.0 可能有 API 变更，建议使用 1.7.x
3. **运行顺序**: Step3 必须在 Step4/5 之前运行（生成匹配 ID 和生存分析结果）
4. **内存**: Step6 的 21 组 PSM 重新匹配较耗内存，建议 ≥8GB RAM
