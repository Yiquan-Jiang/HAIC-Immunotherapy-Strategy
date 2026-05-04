#!/usr/bin/env bash
# 一键运行 TTE cohort_7group_psm02（NLR 规则版）
# R 只写 CSV；本脚本在 R 成功后自动调用 tte_nlr_R_figures.py 生成 PDF/PNG。
# 用法:
#   bash run_tte_cohort_7group_psm02.sh
#   bash run_tte_cohort_7group_psm02.sh /path/to/custom_cohort_ids.csv
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
DATA_DIR="$PROJECT_ROOT/data/tidy_data"
R_OUT="$PROJECT_ROOT/output/step3_tte/NLR_BASED_RULES_R/cohort_7group_psm02"
if [[ $# -ge 1 ]]; then
  Rscript "$SCRIPT_DIR/../tte_core/tte_nlr_R_core_cohort_7group_psm02.R" "$DATA_DIR" "$1"
else
  Rscript "$SCRIPT_DIR/../tte_core/tte_nlr_R_core_cohort_7group_psm02.R" "$DATA_DIR"
fi
python3 "$SCRIPT_DIR/../tte_core/tte_nlr_R_figures.py" "$R_OUT"
