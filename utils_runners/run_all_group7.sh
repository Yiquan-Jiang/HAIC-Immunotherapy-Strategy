#!/usr/bin/env bash
# ============================================================
# update_group_7 — 一键复现全部分析
# ============================================================
#
# 使用方法:
#   cd update_group_7/scripts
#   bash run_all.sh          # 运行全部步骤（含 Step7）
#   bash run_all.sh 3        # 从 Step3 开始运行
#   bash run_all.sh 4 5      # 只运行 Step4 和 Step5
#   bash run_all.sh 7        # 只运行 Step7 全部子步骤
#   bash run_all.sh 7a       # 只运行 Step7a（OW 分析，R）
#   bash run_all.sh 7b       # 只运行 Step7b（亚组可视化，Python）
#   bash run_all.sh 7c       # 只运行 Step7c（基线平衡表，R）
#   bash run_all.sh 7d       # 只运行 Step7d（平衡表可视化，Python）
#
# 前置依赖:
#   Python: pandas, numpy, lifelines, matplotlib, scikit-learn
#   R:      tidyverse, MatchIt, survival, survminer, cobalt,
#           gtsummary, flextable, officer, WeightIt, survey
# ============================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BASE_DIR="$(dirname "$SCRIPT_DIR")"

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
NC='\033[0m'

log()  { echo -e "${CYAN}[$(date '+%H:%M:%S')]${NC} $*"; }
ok()   { echo -e "${GREEN}[✓]${NC} $*"; }
warn() { echo -e "${YELLOW}[!]${NC} $*"; }
fail() { echo -e "${RED}[✗]${NC} $*"; exit 1; }

# 确保输出目录存在
mkdir -p "$BASE_DIR"/{data,figures/{km,loveplots,psm_final,psm_pub_quality,subgroup},results/{psm_balance_tables_complete,psm_vs_template,tables,subgroup_analysis},logs}

# 判断是否运行某个步骤
# 支持：精确匹配（如 "7a"）或前缀匹配（"7" 匹配 7a/7b/7c/7d）
should_run() {
    local step=$1
    if [ ${#STEPS[@]} -eq 0 ]; then
        return 0
    fi
    for s in "${STEPS[@]}"; do
        [ "$s" = "$step" ] && return 0
        # 前缀匹配: "3" 匹配 3b, "4" 匹配 4b, "7" 匹配 7a/7b/7c/7d
        [[ "$step" == "${s}"* ]] && [ ${#s} -lt ${#step} ] && return 0
    done
    return 1
}

STEPS=("$@")

echo ""
echo "============================================================"
echo "  update_group_7 — 7组 HAIC 治疗方案对比分析"
echo "  $(date '+%Y-%m-%d %H:%M:%S')"
echo "============================================================"
echo ""

# ── Step 0: 数据准备 ──────────────────────────────────────────
if should_run 0; then
    log "Step 0: 数据准备（合并 main_group，排除 grey_zone/before_haic）"
    python3 "$SCRIPT_DIR/step0_prepare_data.py"
    ok "Step 0 完成 → data/analysis_ready.csv"
    echo ""
fi

# ── Step 3: PSM 匹配 + 生存分析 ──────────────────────────────
if should_run 3; then
    log "Step 3: PSM 分析（21组两两对比）— 预计 5-10 分钟"
    Rscript "$SCRIPT_DIR/step3_psm_analysis.R"
    ok "Step 3 完成 → results/psm_balance_tables_complete/"
    echo ""
fi

# ── Step 3b: CBPS-IPTW 多组加权（WeightIt, ATE） ──
if should_run 3b; then
    log "Step 3b: CBPS-IPTW 多组加权（WeightIt ATE + 权重截断 + 全模型Cox + PH检验）"
    Rscript "$SCRIPT_DIR/step3b_psm_vs_template.R"
    ok "Step 3b 完成 → results/psm_vs_template/"
    echo ""
fi

# ── Step 4: KM 生存曲线（Python） ────────────────────────────
if should_run 4; then
    log "Step 4: KM 生存曲线（7组整体 + 21组两两 PSM前后）"
    python3 "$SCRIPT_DIR/step4_km_curves.py"
    ok "Step 4 完成 → figures/km/"
    echo ""
fi

# ── Step 4b: PSM 模板匹配后 7 组叠加 KM 曲线 ─────────────
if should_run 4b; then
    log "Step 4b: PSM 模板匹配后 7 组叠加 KM 曲线"
    python3 "$SCRIPT_DIR/step4b_km_template_matched.py"
    ok "Step 4b 完成 → figures/km/km_7groups_template_matched.*"
    echo ""
fi

# ── Step 5: 森林图 ───────────────────────────────────────────
if should_run 5; then
    log "Step 5: 森林图（21组 HR 对比，PSM前后）"
    python3 "$SCRIPT_DIR/step5_forest_plot.py"
    log "Step 5b: 森林图（各组 vs HAIC+I+T_concurrent，仅 PSM 后）"
    python3 "$SCRIPT_DIR/step5b_forest_vs_IT_concurrent.py"
    log "Step 5c: 森林图（各组 vs HAIC_alone，PSM 前后并排）"
    python3 "$SCRIPT_DIR/step5c_forest_vs_HAIC_alone.py"
    ok "Step 5 完成 → figures/psm_pub_quality/"
    echo ""
fi

# ── Step 6: Table 1 + PSM 平衡表 + Love Plot ────────────────
if should_run 6; then
    log "Step 6: 论文表格 + Love Plot — 预计 10-15 分钟"
    Rscript "$SCRIPT_DIR/step6_tables_and_loveplots.R"
    ok "Step 6 完成 → results/tables/ + figures/loveplots/"
    echo ""
fi

# ── Step 7: 高危亚组分析（Overlap Weighting）────────────────
# Step 7a: OW 分析（R）
if should_run 7a; then
    log "Step 7a: 高危亚组 OW 分析（R）— 预计 2 分钟"
    log "  方法: Overlap Weighting + 加权 Cox (robust SE)"
    log "  亚组: Composite / Tumor>3 / Diameter>10cm / PVTT Vp3/4 / Extrahepatic"
    Rscript "$SCRIPT_DIR/step7_subgroup_ow.R"
    ok "Step 7a 完成 → results/subgroup_analysis/ + figures/subgroup/ow_love_*"
    echo ""
fi

# Step 7b: 亚组森林图 + KM 曲线（Python）
if should_run 7b; then
    log "Step 7b: 亚组森林图 + KM 曲线（Python）— 预计 1 分钟"
    python3 "$SCRIPT_DIR/step7_subgroup_plots.py"
    ok "Step 7b 完成 → figures/subgroup/ow_subgroup_forest_plot.* + ow_subgroup_km_curves.*"
    echo ""
fi

# Step 7c: OW 加权后基线平衡表（R）
if should_run 7c; then
    log "Step 7c: OW 加权后基线平衡表（R）— 预计 2 分钟"
    Rscript "$SCRIPT_DIR/step7_ow_balance_table.R"
    ok "Step 7c 完成 → results/subgroup_analysis/ow_balance_table_*.csv"
    echo ""
fi

# Step 7d: 基线平衡表可视化（Python）
if should_run 7d; then
    log "Step 7d: 基线平衡表可视化（Python）— 预计 1 分钟"
    python3 "$SCRIPT_DIR/step7_ow_balance_figure.py"
    ok "Step 7d 完成 → figures/subgroup/ow_balance_table_*.pdf/png"
    echo ""
fi

# ── 清理 R 临时文件 ──────────────────────────────────────────
rm -f "$BASE_DIR/Rplots.pdf" "$BASE_DIR/.RData" "$BASE_DIR/.Rhistory"
rm -f "$SCRIPT_DIR/Rplots.pdf"

echo "============================================================"
echo -e "  ${GREEN}全部完成！${NC}"
echo "============================================================"
echo ""
echo "  输出目录结构:"
echo "  ├── data/analysis_ready.csv               分析数据（3,885例）"
echo "  ├── results/"
echo "  │   ├── psm_balance_tables_complete/       PSM匹配ID + 生存分析"
echo "  │   ├── tables/                            Table1 + 21组平衡表(.docx)"
echo "  │   └── subgroup_analysis/                 ★ Step7 亚组分析结果"
echo "  │       ├── ow_subgroup_results.csv         OW HR/CI/P/等价性"
echo "  │       ├── ow_smd_balance.csv              SMD（OW前后）"
echo "  │       ├── ow_interaction_tests.csv        交互作用 P 值"
echo "  │       └── ow_weighted_ids_*.csv           加权患者数据"
echo "  ├── figures/"
echo "  │   ├── km/                                KM曲线 (PDF+PNG)"
echo "  │   ├── psm_pub_quality/                   森林图 (PDF+PNG)"
echo "  │   ├── loveplots/                         Love Plot (PDF+PNG)"
echo "  │   └── subgroup/                          ★ Step7 亚组图表"
echo "  │       ├── ow_love_*.pdf/png               各亚组 Love Plot"
echo "  │       ├── ow_subgroup_forest_plot.*       亚组森林图"
echo "  │       ├── ow_subgroup_km_curves.*         亚组 KM 曲线"
echo "  │       └── ow_balance_table_*.pdf/png      OW 基线平衡表"
echo "  └── logs/                                  运行日志"
echo ""
