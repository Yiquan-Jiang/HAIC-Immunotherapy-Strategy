#!/bin/bash
# =============================================================================
# Interaction Analysis 运行脚本
# =============================================================================

set -e  # 遇到错误即退出

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "============================================"
echo "Interaction Analysis - PSM Forest Plots"
echo "============================================"

# 检查数据文件
echo ""
echo "检查数据文件..."
required_files=(
    "data/HAIC_NO_TACE_4_TIDY_baseline.csv"
    "data/HAIC_NO_TACE_4_TIDY_longitudinal.csv"
    "data/matched_ids_02_HAIC_alone_vs_HAIC_then_I.csv"
    "data/matched_ids_06_HAIC_alone_vs_HAIC_then_I+T.csv"
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
    echo "错误: 缺少必要的数据文件，请先准备数据。"
    exit 1
fi

# 运行PSM02分析
echo ""
echo "============================================"
echo "1. 运行 PSM02 分析 (HAIC alone vs HAIC then I)"
echo "============================================"
python 01_publication_figures.py

# 运行PSM06分析
echo ""
echo "============================================"
echo "2. 运行 PSM06 分析 (HAIC alone vs HAIC then I+T)"
echo "============================================"
python 02_publication_figures_ids06_IplusT.py

echo ""
echo "============================================"
echo "分析完成！"
echo "输出目录: $SCRIPT_DIR/output/"
echo "============================================"