#!/bin/bash
# =============================================================================
# RCS Non-linear Interaction Analysis 运行脚本
# 分析 HAIC_alone 与 6 个治疗组的 RCS 非线性交互
# =============================================================================

set -e  # 遇到错误即退出

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "============================================"
echo "RCS Non-linear Interaction Analysis"
echo "HAIC_alone vs 6 Treatment Groups"
echo "============================================"

# 检查Python依赖
echo ""
echo "检查Python依赖..."
python3 -c "import pandas; import numpy" 2>/dev/null || {
    echo "错误: 缺少Python依赖 (pandas, numpy)"
    echo "请运行: pip install pandas numpy"
    exit 1
}
echo "  ✓ Python依赖已安装"

# 检查R依赖
echo ""
echo "检查R依赖..."
Rscript -e 'suppressPackageStartupMessages({library(survival); library(rms); library(ggplot2)})' 2>/dev/null || {
    echo "错误: 缺少R依赖"
    echo "请运行: install.packages(c('survival', 'rms', 'ggplot2', 'dplyr', 'gridExtra', 'glmnet'))"
    exit 1
}
echo "  ✓ R依赖已安装"

# 检查数据文件
echo ""
echo "检查数据文件..."
required_files=(
    "data/HAIC_NO_TACE_4_TIDY_baseline_imputed.csv"
    "data/HAIC_NO_TACE_4_TIDY_longitudinal.csv"
)

missing=0
for f in "${required_files[@]}"; do
    if [[ -f "$f" ]]; then
        echo "  ✓ $f"
    else
        echo "  ✗ $f (缺失)"
        missing=1
    fi
done

if [[ $missing -eq 1 ]]; then
    echo ""
    echo "警告: 部分数据文件缺失，将尝试在分析时自动生成。"
fi

# 定义分析队列
declare -A COHORTS=(
    ["THEN_I"]="HAIC alone vs HAIC then I"
    ["THEN_IT"]="HAIC alone vs HAIC then I+T"
    ["THEN_T"]="HAIC alone vs HAIC then T"
    ["I_CONC"]="HAIC alone vs HAIC+I concurrent"
    ["T_CONC"]="HAIC alone vs HAIC+T concurrent"
    ["I_T_CONC"]="HAIC alone vs HAIC+I+T concurrent"
)

# 函数：运行单个队列分析
run_cohort() {
    local key=$1
    local desc=$2
    local py_script="RCS_INT_HAIC_ALONE_AND_${key}_build_cohort.py"
    local r_script="RCS_INT_HAIC_ALONE_AND_${key}_dual_timescale.R"
    local csv_file="data/RCS_INT_HAIC_ALONE_AND_${key}_cohort.csv"
    
    echo ""
    echo "============================================"
    echo "[$key] $desc"
    echo "============================================"
    
    # 构建队列
    if [[ ! -f "$csv_file" ]]; then
        echo "构建队列..."
        python3 "$py_script"
    else
        echo "队列已存在: $csv_file"
    fi
    
    # 运行RCS分析
    echo "运行 RCS 分析..."
    Rscript "$r_script"
}

# 运行所有队列
idx=1
for key in "${!COHORTS[@]}"; do
    echo ""
    echo "============================================"
    echo "$idx. ${COHORTS[$key]}"
    echo "============================================"
    run_cohort "$key" "${COHORTS[$key]}"
    ((idx++))
done

# AFP-PIVKA 组合分析
echo ""
echo "============================================"
echo "7. AFP-PIVKA 组合指标分析 (可选)"
echo "============================================"
read -p "是否运行 AFP-PIVKA 组合分析? (y/n): " run_composite

if [[ "$run_composite" == "y" || "$run_composite" == "Y" ]]; then
    cd afp_pivka_composite
    
    # 构建组合队列
    echo "构建组合队列..."
    python3 00_build_composite_cohorts.py
    
    # 运行RCS分析
    echo "运行 AFP-PIVKA RCS 分析..."
    Rscript 01_rcs_afp_pivka_composite.R ALL
fi

echo ""
echo "============================================"
echo "分析完成！"
echo "输出目录: $SCRIPT_DIR/output/"
echo "============================================"
echo ""
echo "各队列输出位置:"
echo "  - output/RCS_INT_HAIC_ALONE_AND_THEN_I/"
echo "  - output/RCS_INT_HAIC_ALONE_AND_THEN_IT/"
echo "  - output/RCS_INT_HAIC_ALONE_AND_THEN_T/"
echo "  - output/RCS_INT_HAIC_ALONE_AND_I_CONC/"
echo "  - output/RCS_INT_HAIC_ALONE_AND_T_CONC/"
echo "  - output/RCS_INT_HAIC_ALONE_AND_I_T_CONC/"
echo ""