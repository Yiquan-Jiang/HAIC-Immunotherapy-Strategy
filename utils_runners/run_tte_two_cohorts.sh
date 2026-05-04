#!/usr/bin/env bash
# =============================================================================
# Target Trial Emulation — IT_RULES_v2 applied to two cohorts (Cohort A + B)
# =============================================================================
# Driver: tte_IT_R_two_cohorts.R   (R: CCW + IPCW + weighted Cox + RMST → CSV)
# Figures: tte_IT_R_figures_two_cohorts.py   (Python: KM, Δ-RMST forest, IPCW
#                                             diagnostics, sensitivity panels)
#
# Cohorts handled internally by tte_IT_R_two_cohorts.R:
#   (A) cohort_3matched     — matched_06 + all HAIC+I+T concurrent
#                             Add-on: ICI + antiangiogenic
#                             Triggers: PLR > 102.4, AFP drop < 32.5 %,
#                                       SII / PIVKA / lymph node
#   (B) cohort_7group_psm02 — matched_02 + all HAIC+I concurrent
#                             Add-on: ICI only
#                             Triggers: PLR > 98.7, AFP drop < 40 %,
#                                       NLR / PIVKA / lymph node
#
# Output base: <project>/output/step3_tte/IT_RULES_R_two_cohorts/{cohort_3matched,
#                                                                  cohort_7group_psm02}/
#
# Usage:
#   bash utils_runners/run_tte_two_cohorts.sh [DATA_DIR]
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

DATA_DIR="${1:-${PROJECT_ROOT}/data}"
TTE_CORE="${PROJECT_ROOT}/07_target_trial_emulation/tte_IT_R_two_cohorts.R"
TTE_FIG="${PROJECT_ROOT}/07_target_trial_emulation/tte_IT_R_figures_two_cohorts.py"

if [[ ! -d "$DATA_DIR" ]]; then
  echo "[error] DATA_DIR does not exist: $DATA_DIR"
  echo "        Pass an explicit data dir as first argument or place inputs in ./data/"
  exit 1
fi

echo "[1/2] R core: tte_IT_R_two_cohorts.R"
Rscript "$TTE_CORE" "$DATA_DIR"

echo "[2/2] Python figures: tte_IT_R_figures_two_cohorts.py"
python3 "$TTE_FIG"

echo
echo "Done. Outputs under: ${PROJECT_ROOT}/output/step3_tte/IT_RULES_R_two_cohorts/"
