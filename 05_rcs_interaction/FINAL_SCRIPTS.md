# RCS 非线性交互分析 — 最终脚本索引

本目录下有多个版本的 RCS 分析脚本。**以下 4 个是研究最终采用的版本**，其他脚本为中间/归档版本，不用于正式结果。

## 最终版本（2 × 2 矩阵）

|  | **单指标组合图**<br/>（每个变量一张 combined 图） | **8×5 矩阵总览图**<br/>（一张大图覆盖 8 指标 × 5 时间列） |
|---|---|---|
| **路线 A：PSM 1:1 匹配，无权重 Cox**<br/>6 个 pair 对子 | `RCS_PSM_dual_timescale.R`<br/>（2026-04-17） | `RCS_PSM_matrix_panel.R`<br/>（2026-04-18） |
| **路线 B：composite 队列，IPTW 加权 Cox**<br/>5 个 arm vs HAIC_alone | `afp_pivka_composite/01_rcs_afp_pivka_composite.R`<br/>（2026-04-17） | `afp_pivka_composite/02_rcs_matrix_panel.R`<br/>（2026-04-18） |

- **路线 A vs B**：两种独立方法论（PSM 匹配 vs IPTW 加权），互为敏感性分析
- **single vs matrix**：同一路线内的两种出图形式
  - `single_indicator` → 分变量组合图（每个生物标志物一张图）
  - `matrix_panel` → 8×5 总览图（AFP / PIVKA / PIV / SII / NLR / PLR / MONOCYTE / ALBI × Baseline / Pre-HAIC-3 / Pre-IT / Pre-HAIC-3 Change Rate / Pre-IT Change Rate）

## 输出位置

- 路线 A → `output/step1_rcs_interaction/psm/pair_XX/{landmark,total_os}/`
- 路线 B → `output/step1_rcs_interaction/relaxed_dual/<cohort_key>/{landmark,total_os}/`

## 运行命令

```bash
# 前置：先构建队列
python build_cohort_psm.py                                # 路线 A 队列
python afp_pivka_composite/00_build_composite_cohorts.py  # 路线 B 队列

# 路线 A：PSM + 无权重
Rscript RCS_PSM_dual_timescale.R                          # 分变量图
Rscript RCS_PSM_matrix_panel.R                            # 矩阵大图

# 路线 B：IPTW 加权
Rscript afp_pivka_composite/01_rcs_afp_pivka_composite.R ALL
Rscript afp_pivka_composite/02_rcs_matrix_panel.R ALL
```

环境变量（可选）：`RMS_RCS_NK`（默认 3），`RMS_RCS_N_BOOT`（默认 200）

## 非最终版本（请勿用于正式结果）

- `rcs_all_pairs_dual_timescale.R` — 早期 all-pairs 版本，已被路线 A/B 拆分取代
- `build_all_pairs_cohorts.py` — 与上配套的旧队列构建脚本
- `extract_pre_it_for_psm.py` — pre-IT 实验室值抽取（辅助脚本，非主分析）
