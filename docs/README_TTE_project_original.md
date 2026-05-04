# HAIC 按需加入免疫/靶向治疗的决策策略研究

## 一、研究概述

本项目探索**肝动脉灌注化疗（HAIC）**治疗肝细胞癌（HCC）过程中，**按需加入免疫治疗或靶向治疗**的最优策略。

### 核心研究问题

> 在 HAIC 治疗过程中，哪些患者应该加入免疫/靶向治疗？何时加入？依据什么指标决策？

### 研究路线（三步走）

```
Step 1: RCS 交互分析 → 发现连续变量的非线性效应修饰作用，确定切割点
    ↓
Step 2: 交互森林图 → 验证亚组异质性，确认哪些亚组从联合治疗中获益
    ↓
Step 3: 目标试验模拟（TTE）→ 用因果推断框架正式比较"动态决策策略"vs"早期联合策略"
```

---

## 二、分析流程详解

### Step 1: RCS 交互分析（Restricted Cubic Splines × Treatment Interaction）

**目的**：探索连续生物标志物是否以非线性方式修饰治疗效应（HAIC alone vs 各联合方案），并找到 HR 曲线与参考线的交叉点（即潜在切割阈值）。

**方法**：

- 以 HAIC alone 为参照，与其余 6 种方案做两两配对
- 对每个连续修饰变量拟合：`Surv(time, event) ~ trt × rcs(x, nk=3)`
- IPTW 加权（relaxed dual）或 PSM 匹配后未加权（PSM 版）
- Bootstrap（200次）得到 HR 的 95% 置信带
- 双时间尺度：42天 Landmark 后残存 OS + 总 OS

**关键输出**：

- HR 随连续变量变化的曲线图（含 CI 带）
- ANOVA 交互 P 值（整体交互 + 非线性交互）
- HR 与参考线（1.0, 0.85, 0.7）的交叉点 → **用于确定后续分析的切割阈值**

**脚本**：


| 文件                                                       | 说明                            |
| -------------------------------------------------------- | ----------------------------- |
| `scripts/rcs_interaction/build_all_pairs_cohorts.py`     | 构建 6 个两两配对队列（IPTW 版）          |
| `scripts/rcs_interaction/rcs_all_pairs_dual_timescale.R` | **主分析**：RCS×治疗交互 Cox（IPTW 加权） |
| `scripts/rcs_interaction/build_cohort_psm.py`            | 构建 6 个 PSM 匹配队列               |
| `scripts/rcs_interaction/RCS_PSM_dual_timescale.R`       | **敏感性分析**：PSM 队列上的 RCS 交互     |
| `scripts/rcs_interaction/afp_pivka_composite/`           | AFP-PIVKA 复合指标的 RCS 交互分析      |


### Step 2: 交互森林图（Interaction Forest Plot）

**目的**：将 Step 1 发现的连续变量切割点应用于亚组分析，生成发表级交互森林图，展示各亚组的治疗效应异质性。

**方法**：

- 基于 PSM 匹配队列（HAIC alone vs HAIC then I / HAIC then I+T）
- 对每个亚组拟合 Cox 模型，计算 HR + 95% CI
- 计算 ΔRMST（受限平均生存时间差）+ Bootstrap CI
- 似然比检验（LRT）得到 P_interaction

**关键输出**：

- Fig4 森林图：亚组 HR + ΔRMST + P_interaction
- Landmark 补充森林图
- 亚组分析汇总表

**脚本**：


| 文件                                                                   | 说明                                   |
| -------------------------------------------------------------------- | ------------------------------------ |
| `scripts/publication_figures/01_publication_figures.py`              | PSM02 队列：HAIC alone vs HAIC then I   |
| `scripts/publication_figures/02_publication_figures_ids06_IplusT.py` | PSM06 队列：HAIC alone vs HAIC then I+T |


### Step 3: 目标试验模拟（Target Trial Emulation, TTE）

**目的**：用因果推断框架正式比较两种治疗策略：

- **策略 A（动态策略）**：根据 NLR/PIV 等指标的动态变化，在满足触发条件时才加入免疫治疗
- **策略 B（早期联合）**：HAIC 后 14 天内即加入免疫治疗

**方法**：

- **Clone-Censor-Weight (CCW)** 框架：每人克隆为两条观测，按策略规则人工删失
- **稳定化 IPCW**：对人为删失用逆概率加权（person-period pooled logistic）
- **加权 Cox** + 稳健三明治 SE + PH 检验
- **RMST**：加权 KM 阶梯积分 + Bootstrap CI（500次，每次重估 IPCW）
- **敏感性分析**：权重截断、IPTW×IPCW、不同 grace period、E-value

**动态策略触发规则（NLR 版）**：

```
绝对触发（pre-HAIC-3 时评估）：
  Rule 1: 基线 PVTT Vp3/4 / 远处转移 / HVTT / IVC-RA癌栓 / 淋巴结转移
  Rule 2: 基线或 pre-HAIC-3 NLR ≥ 2.5
  Rule 3: 基线最大肿瘤径 > 13 cm
  Rule 4: 基线 PIVKA > 12000 mAU/mL
  Rule 5: 基线 AFP < 20 ng/mL

豁免规则（cycle ≥ 4 时评估）：
  ExRule 2: NLR 持续 < 2.5 → 继续观察
  ExRule 3: AFP 较基线下降 > 50% → 继续观察
```

**PIV 版**：用 PIV（PLT×单核×中性粒/淋巴）替代 NLR 作为核心触发指标。

**关键输出**：

- R 导出 CSV（HR、RMST、敏感性、克隆数据集、KM 曲线、风险表）
- Python 生成发表级图：KM 曲线、ΔRMST 森林图、流程图、权重诊断图、敏感性森林图

**脚本**：


| 文件                                                      | 说明                       |
| ------------------------------------------------------- | ------------------------ |
| `scripts/tte_core/tte_nlr_R_core_cohort_7group_psm02.R` | TTE 核心分析（NLR 规则版）        |
| `scripts/tte_core/tte_piv_R_core_cohort_7group_psm02.R` | TTE 核心分析（PIV 规则版）        |
| `scripts/tte_core/tte_nlr_R_figures.py`                 | TTE 结果可视化（通用，读 R 输出 CSV） |
| `scripts/runners/run_tte_cohort_7group_psm02.sh`        | 一键运行 NLR 版（R + Python）   |
| `scripts/runners/run_tte_cohort_7group_psm02_piv.sh`    | 一键运行 PIV 版               |


### 辅助：策略示意图


| 文件                                                             | 说明                   |
| -------------------------------------------------------------- | -------------------- |
| `scripts/schematic_figures/plot_tte_schematic_v3.py`           | TTE 策略示意图（Nature 风格） |
| `scripts/schematic_figures/plot_tte_nlr_strategy_flowchart.py` | NLR 策略流程图            |


---

## 三、数据说明

### data/tidy_data/（核心数据）


| 文件                                         | 说明                        | 行数     |
| ------------------------------------------ | ------------------------- | ------ |
| `HAIC_NO_TACE_4_TIDY_baseline.csv`         | 全队列基线数据（无 TACE 的 HAIC 患者） | ~4235  |
| `HAIC_NO_TACE_4_TIDY_baseline_imputed.csv` | 基线数据（缺失值已插补）              | ~4235  |
| `HAIC_NO_TACE_4_TIDY_longitudinal.csv`     | 纵向随访数据（各 HAIC 周期前检查）      | ~43672 |
| `analysis_ready.csv`                       | 7 组分类后的分析就绪数据             | ~3887  |


> **注**：原 `HAIC_IMMUNE_baseline.xlsx` 和 `HAIC_IMMUNE_longitudinal.xlsx` 已被 TIDY CSV 替代。
> TIDY CSV 是 xlsx 的超集（4234 vs 1614 人，68 vs 59 列），所有脚本已适配。
> 列名映射：`days_haic_to_immune_y` → `days_haic_to_immune`，`albi_score_calculated` → `albi_score`。

### data/psm_matched_ids/（PSM 匹配 ID）

6 个文件，对应 HAIC alone 与其余 6 种联合方案的两两 PSM 匹配结果。

### data/rcs_cohorts/（RCS 分析队列）

- `relaxed_dual/`：IPTW 版 6 个配对队列 + AFP-PIVKA 复合队列
- `psm/`：PSM 版 6 个配对队列

### data/publication_cohorts/（发表图用数据）


| 文件                                               | 说明                           |
| ------------------------------------------------ | ---------------------------- |
| `00_swimmer_plot_events.csv`                     | Swimmer plot 事件数据（含 TACE 信息） |
| `matched_ids_06_HAIC_alone_vs_HAIC_then_I+T.csv` | PSM06 匹配 ID                  |


---

## 四、7 组治疗方案


| 编号  | 方案                    | 说明             |
| --- | --------------------- | -------------- |
| 1   | HAIC alone            | 单纯 HAIC        |
| 2   | HAIC + I concurrent   | HAIC 同期联合免疫    |
| 3   | HAIC then I           | HAIC 后序贯免疫     |
| 4   | HAIC + T concurrent   | HAIC 同期联合靶向    |
| 5   | HAIC then T           | HAIC 后序贯靶向     |
| 6   | HAIC + I+T concurrent | HAIC 同期联合免疫+靶向 |
| 7   | HAIC then I+T         | HAIC 后序贯免疫+靶向  |


---

## 五、复现指南

### 环境依赖

**R 包**：

```r
install.packages(c("dplyr", "tidyr", "stringr",
                    "survival", "survey", "cobalt", "boot",
                    "rms", "glmnet", "ggplot2"))
```

**Python 包**：

```bash
pip install numpy pandas matplotlib lifelines scipy
```

### 执行顺序

> 所有脚本的输出统一写入项目根目录下的 `output/`，按分析步骤分子目录组织。

#### A. RCS 交互分析（Step 1）

```bash
# 1. 构建 IPTW 版配对队列
cd scripts/rcs_interaction/
python3 build_all_pairs_cohorts.py

# 2. 运行 RCS 交互分析（IPTW 版）→ 输出至 output/step1_rcs_interaction/iptw/
Rscript rcs_all_pairs_dual_timescale.R

# 3. 构建 PSM 版配对队列
python3 build_cohort_psm.py

# 4. 运行 RCS 交互分析（PSM 版）→ 输出至 output/step1_rcs_interaction/psm/
Rscript RCS_PSM_dual_timescale.R
```

> **可选覆盖**：`RCS_OUT_ROOT` 环境变量可自定义 IPTW 版输出根目录。

#### B. 交互森林图（Step 2）

```bash
cd scripts/publication_figures/
python3 01_publication_figures.py      # → output/step2_interaction_forest/psm02_HAIC_then_I/
python3 02_publication_figures_ids06_IplusT.py  # → output/step2_interaction_forest/psm06_HAIC_then_IplusT/
```

#### C. 目标试验模拟（Step 3）

```bash
# NLR 规则版（推荐使用一键脚本）→ output/step3_tte/NLR_BASED_RULES_R/cohort_7group_psm02/
bash scripts/runners/run_tte_cohort_7group_psm02.sh

# PIV 规则版 → output/step3_tte/PIV_BASED_RULES_R/cohort_7group_psm02/
bash scripts/runners/run_tte_cohort_7group_psm02_piv.sh
```

> 也可手动分步执行：
>
> ```bash
> Rscript scripts/tte_core/tte_piv_R_core_cohort_7group_psm02.R data/tidy_data
> python3 scripts/tte_core/tte_nlr_R_figures.py output/step3_tte/PIV_BASED_RULES_R/cohort_7group_psm02
> ```

---

## 六、关键统计方法总结


| 方法                          | 用途       | 实现                              |
| --------------------------- | -------- | ------------------------------- |
| RCS × Treatment Interaction | 非线性效应修饰  | `rms::cph` + `rcs()`            |
| IPTW（Ridge PS）              | 治疗组平衡    | `glmnet` α=0, 5-fold CV         |
| PSM                         | 治疗组匹配    | 外部已完成（`matched_ids_*.csv`）      |
| Cox PH Model                | 生存分析     | `survival::coxph` / `lifelines` |
| RMST                        | 受限平均生存时间 | 加权 KM 阶梯积分                      |
| Clone-Censor-Weight         | 目标试验模拟   | 手写 CCW + IPCW                   |
| Stabilized IPCW             | 人工删失校正   | Person-period pooled logistic   |
| Bootstrap                   | 置信区间     | 500次（TTE）/ 200次（RCS）/ 300次（森林图） |
| E-value                     | 未测量混杂敏感性 | 主分析附加                           |


---

## 七、目录结构

```
HAIC_Immunotherapy_Decision_TTE/
├── README.md                          ← 本文件
├── scripts/                           ← 所有脚本（只读逻辑，不存放输出）
│   ├── rcs_interaction/               ← Step 1: RCS 交互分析
│   │   ├── build_all_pairs_cohorts.py
│   │   ├── rcs_all_pairs_dual_timescale.R
│   │   ├── build_cohort_psm.py
│   │   ├── RCS_PSM_dual_timescale.R
│   │   └── afp_pivka_composite/
│   │       ├── 00_build_composite_cohorts.py
│   │       └── 01_rcs_afp_pivka_composite.R
│   ├── publication_figures/           ← Step 2: 交互森林图
│   │   ├── 01_publication_figures.py
│   │   └── 02_publication_figures_ids06_IplusT.py
│   ├── tte_core/                      ← Step 3: 目标试验模拟
│   │   ├── tte_nlr_R_core_cohort_7group_psm02.R
│   │   ├── tte_piv_R_core_cohort_7group_psm02.R
│   │   └── tte_nlr_R_figures.py
│   ├── runners/                       ← 一键运行脚本
│   │   ├── run_tte_cohort_7group_psm02.sh
│   │   └── run_tte_cohort_7group_psm02_piv.sh
│   └── schematic_figures/             ← 策略示意图
│       ├── plot_tte_schematic_v3.py
│       └── plot_tte_nlr_strategy_flowchart.py
├── data/                              ← 输入数据（只读）
│   ├── tidy_data/                     ← 核心数据
│   │   ├── HAIC_NO_TACE_4_TIDY_baseline.csv
│   │   ├── HAIC_NO_TACE_4_TIDY_baseline_imputed.csv
│   │   ├── HAIC_NO_TACE_4_TIDY_longitudinal.csv
│   │   └── analysis_ready.csv
│   ├── psm_matched_ids/               ← PSM 匹配结果
│   │   ├── matched_ids_01_HAIC_alone_vs_HAIC+I_concurrent.csv
│   │   ├── matched_ids_02_HAIC_alone_vs_HAIC_then_I.csv
│   │   ├── matched_ids_03_HAIC_alone_vs_HAIC+T_concurrent.csv
│   │   ├── matched_ids_04_HAIC_alone_vs_HAIC_then_T.csv
│   │   ├── matched_ids_05_HAIC_alone_vs_HAIC+I+T_concurrent.csv
│   │   └── matched_ids_06_HAIC_alone_vs_HAIC_then_I+T.csv
│   ├── rcs_cohorts/                   ← RCS 分析用队列
│   │   ├── relaxed_dual/              ← IPTW 版
│   │   └── psm/                       ← PSM 版
│   └── publication_cohorts/           ← 发表图用数据
│       ├── 00_swimmer_plot_events.csv
│       └── matched_ids_06_HAIC_alone_vs_HAIC_then_I+T.csv
└── output/                            ← 所有脚本的运行输出（图表 + 结果 CSV）
    ├── step1_rcs_interaction/         ← Step 1 输出
    │   ├── iptw/                      ← IPTW 加权 RCS 交互图表
    │   │   └── <pair_label>/{landmark,total_os}/
    │   ├── psm/                       ← PSM 版 RCS 交互图表
    │   │   └── pair_XX_*/{landmark,total_os}/
    │   └── afp_pivka_composite/       ← AFP-PIVKA 复合指标 RCS 图表
    │       └── <cohort_key>/{landmark,total_os}/
    ├── step2_interaction_forest/      ← Step 2 输出
    │   ├── psm02_HAIC_then_I/        ← PSM02 森林图 + 表格
    │   └── psm06_HAIC_then_IplusT/   ← PSM06 森林图 + 表格
    ├── step3_tte/                     ← Step 3 输出
    │   ├── NLR_BASED_RULES_R/        ← NLR 版 TTE 结果 CSV + 发表级图
    │   │   └── cohort_7group_psm02/
    │   └── PIV_BASED_RULES_R/        ← PIV 版 TTE 结果 CSV + 发表级图
    │       └── cohort_7group_psm02/
    └── schematic_figures/             ← 策略示意图（PDF + PNG）
```

---

## 八、原始文件位置对照


| 本项目路径                                                    | 原始位置                                                                                               |
| -------------------------------------------------------- | -------------------------------------------------------------------------------------------------- |
| `scripts/rcs_interaction/rcs_all_pairs_dual_timescale.R` | `data/immunotherapy_decision_analysis_副本/true_pre_haic_3/nonlinear_analysis/rms_rcs_relaxed_dual/` |
| `scripts/rcs_interaction/RCS_PSM_dual_timescale.R`       | `data/immunotherapy_decision_analysis_副本/true_pre_haic_3/nonlinear_analysis/rms_rcs_psm/`          |
| `scripts/publication_figures/01_publication_figures.py`  | `data/immunotherapy_decision_analysis_副本/true_pre_haic_3/publication_relaxed/`                     |
| `scripts/tte_core/tte_nlr_R_core_cohort_7group_psm02.R`  | `data/immunotherapy_decision_analysis_副本/script/`                                                  |
| `scripts/tte_core/tte_piv_R_core_cohort_7group_psm02.R`  | `data/immunotherapy_decision_analysis_副本/script/`                                                  |
| `scripts/tte_core/tte_nlr_R_figures.py`                  | `data/immunotherapy_decision_analysis_副本/script/`                                                  |
| `data/tidy_data/HAIC_NO_TACE_4_TIDY_*.csv`               | `data/immunotherapy_decision_analysis_副本/true_pre_haic_3/nonlinear_analysis/rms_rcs_relaxed_dual/` |
| `data/psm_matched_ids/`                                  | `HAIC_NO_TACE_4_TIDY/update_group_7/results/psm_balance_tables_complete/`                          |
| `data/tidy_data/analysis_ready.csv`                      | `HAIC_NO_TACE_4_TIDY/update_group_7/data/`                                                         |


---

## 九、重要提示

1. **输出统一**：所有脚本的图表和结果 CSV 统一输出到 `output/` 目录，按 `step1_rcs_interaction/`、`step2_interaction_forest/`、`step3_tte/`、`schematic_figures/` 分子目录组织。`data/` 目录仅存放输入数据，`scripts/` 目录仅存放脚本。
2. **数据统一**：所有脚本已适配 TIDY CSV 格式，不再依赖 `HAIC_IMMUNE_*.xlsx`。TTE 脚本内部自动做列名兼容映射（`days_haic_to_immune` → `days_haic_to_immune_y` 等）。
3. **PSM 匹配**：PSM 匹配已在上游完成（`matched_ids_*.csv`），本项目直接使用匹配结果。
4. **双时间尺度**：所有生存分析均在两个时间尺度上进行——42天 Landmark 后残存 OS（主要）和总 OS（敏感性）。
5. **TTE 流程图中的硬编码**：`tte_nlr_R_figures.py` 中总数据库 N=1614 为硬编码，若数据更新需同步修改。
6. **AFP-PIVKA 复合分析**：为 RCS 交互分析的扩展，探索 AFP 与 PIVKA 的组合指标作为效应修饰因子。

