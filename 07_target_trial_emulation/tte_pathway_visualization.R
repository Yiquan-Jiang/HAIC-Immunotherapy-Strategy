#!/usr/bin/env Rscript
# =============================================================================
# TTE companion: Real-world treatment pathways for the two TTE cohorts
# =============================================================================
# Purpose : For each cohort produced by tte_IT_R_two_cohorts.R, visualise the
#           ACTUAL treatment paths patients followed — assessing how often the
#           "Early Combination" vs "Adaptive On Demand" arms were realised in
#           practice, and whether downstream curative conversion occurred.
#
# Inputs (resolved from --data_dir, default = sibling of scripts/tte_core/):
#   Per cohort (output/step3_tte/IT_RULES_R_two_cohorts/<cohort>/):
#     tte_cohort_<name>_IT_ids.csv       — cohort membership
#     R_trigger_table.csv                — IT-Rules trigger day per patient
#     R_clone_dataset.csv                — CCW clones (for arm allegiance)
#   Shared:
#     HAIC_NO_TACE_4_TIDY/update_group_7/data/00_swimmer_plot_events.csv
#     HAIC_NO_TACE_4_TIDY/update_group_7/data/00_patient_treatment_summary.csv
#     HAIC_NO_TACE_4_TIDY_baseline.csv   — death/OS, first_haic_date
#
# Outputs (per cohort, written next to TTE results):
#   pathway_A_sankey.{pdf,png}           — Alluvial: arm → timing → curative
#   pathway_classification_table.csv     — Per-patient final labels (audit trail)
#   pathway_pathway_counts.csv           — Sankey flow counts
#
# Pathway taxonomy (3 axes):
#   Arm allegiance : Early Combination | Adaptive On Demand | Non-adherent
#   Timing         : Early <=14d / Bridging 15-42d / During-HAIC >42d <=lastHAIC /
#                    Maintenance >lastHAIC / Never
#   Curative event : Resection / Ablation / Liver Transplant / None
# =============================================================================

Sys.setenv(LANG = "en_US.UTF-8")
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(ggalluvial)
  library(ggrepel)
  library(scales)
  library(forcats)
  library(patchwork)
})

# ── CLI / path resolution ────────────────────────────────────────────────────
args_all <- commandArgs(trailingOnly = TRUE)
data_dir <- args_all[1]
get_script_dir <- function() {
  args_raw <- commandArgs(trailingOnly = FALSE)
  fa <- args_raw[grepl("^--file=", args_raw)]
  if (length(fa)) return(normalizePath(dirname(sub("^--file=", "", fa[1])),
                                       winslash = "/"))
  if (!is.null(sys.frames()[[1]]$ofile))
    return(normalizePath(dirname(sys.frames()[[1]]$ofile), winslash = "/"))
  NULL
}
if (is.na(data_dir) || data_dir == "") {
  sd <- get_script_dir()
  if (!is.null(sd)) {
    data_dir <- normalizePath(file.path(sd, "..", ".."), winslash = "/")
  } else {
    data_dir <- normalizePath(getwd(), winslash = "/")
  }
}
data_dir <- normalizePath(data_dir, winslash = "/", mustWork = TRUE)

project_root <- data_dir
base_out_dir <- file.path(project_root, "output", "step3_tte", "IT_RULES_R_two_cohorts")
if (!dir.exists(base_out_dir))
  stop("TTE output dir not found — run tte_IT_R_two_cohorts.R first: ", base_out_dir)

# Locate the events CSVs (HAIC_NO_TACE_4_TIDY data dir, two possible locations)
events_dir <- file.path(
  normalizePath(file.path(data_dir, ".."), winslash = "/", mustWork = TRUE),
  "HAIC_NO_TACE_4_TIDY", "update_group_7", "data"
)
if (!dir.exists(events_dir)) {
  events_dir_alt <- file.path(data_dir, "HAIC_NO_TACE_4_TIDY", "update_group_7", "data")
  if (dir.exists(events_dir_alt)) events_dir <- events_dir_alt
}
swimmer_path <- file.path(events_dir, "00_swimmer_plot_events.csv")
summary_path <- file.path(events_dir, "00_patient_treatment_summary.csv")
if (!file.exists(swimmer_path)) stop("Missing: ", swimmer_path)
if (!file.exists(summary_path)) stop("Missing: ", summary_path)

baseline_path <- file.path(data_dir, "HAIC_NO_TACE_4_TIDY_baseline.csv")
if (!file.exists(baseline_path)) stop("Missing baseline: ", baseline_path)

cat("======================================================================\n")
cat("TTE pathway visualization\n")
cat("======================================================================\n")
cat("data_dir   : ", data_dir, "\n")
cat("base_out   : ", base_out_dir, "\n")
cat("events_dir : ", events_dir, "\n")

# ── Cohort registry (mirrors tte_IT_R_two_cohorts.R) ─────────────────────────
COHORTS <- list(
  cohort_3matched = list(
    out_subdir = "cohort_3matched",
    label      = "HAIC alone vs HAIC+I+T (matched_06)",
    regimen_target = "I+T"
  ),
  cohort_7group_psm02 = list(
    out_subdir = "cohort_7group_psm02",
    label      = "HAIC alone vs HAIC+I (matched_02)",
    regimen_target = "I"
  )
)

# ── Constants matching the TTE protocol ──────────────────────────────────────
EARLY_GRACE_DAYS    <- 14
BRIDGING_END_DAYS   <- 42
DYNAMIC_GRACE_DAYS  <- 90
FLOW_LABEL_THRESHOLD <- 0.05   # show "n (%)" on flows >= 5% of cohort

# ── Colour palettes (Lancet-inspired) ────────────────────────────────────────
COL_ARM <- c(
  "Early Combination"   = "#00468B",
  "Adaptive On Demand"  = "#925E9F",
  "Non-adherent"        = "#9E9E9E"
)
# Sequential green palette for the time-of-add-on axis: darker = earlier
# clinical engagement. Distinct from the arm-blue and outcome-red/orange
# families so the three columns are visually separable.
COL_TIMING <- c(
  "Early concurrent (<=42d)"         = "#1B4332",
  "During-HAIC (>42d, <= last HAIC)" = "#52B788",
  "Maintenance (> last HAIC)"        = "#B7E4C7",
  "Never"                            = "#E5E5E5"
)
# Curative outcomes: bold NEJM/Lancet colours for conversion paths,
# very light grey for "None" so it visually fades into the background.
COL_CURATIVE <- c(
  "Resection"        = "#BC3C29",
  "Ablation"         = "#E18727",
  "Liver transplant" = "#0072B5",
  "None"             = "#DCDCDC"
)
# Muted-text colours used on right-margin endpoint cards
# (so "No curative event" doesn't visually compete with conversion endpoints)
COL_CURATIVE_TEXT <- COL_CURATIVE
COL_CURATIVE_TEXT["None"] <- "grey55"

# Inline-text colours for the middle (timing) column: white text on the
# dark-green strata, dark text on the lighter-green / grey strata, so
# labels are readable on every stratum without an obscuring box.
COL_TIMING_TEXT <- c(
  "Early concurrent (<=42d)"         = "white",
  "During-HAIC (>42d, <= last HAIC)" = "white",
  "Maintenance (> last HAIC)"        = "grey15",
  "Never"                            = "grey20"
)

# Short display names used on inline middle-column labels (full definitions
# go in the figure caption — keeps the inline text readable inside strata).
TIMING_SHORT <- c(
  "Early concurrent (<=42d)"         = "Early concurrent",
  "During-HAIC (>42d, <= last HAIC)" = "During HAIC",
  "Maintenance (> last HAIC)"        = "Maintenance",
  "Never"                            = "Never added"
)

# Combined lookup: every stratum across all three columns has a colour.
# Used to fill the stratum boxes via geom_rect (one colour per stratum).
STRAT_COLORS <- c(COL_ARM, COL_TIMING, COL_CURATIVE)

# ── Helpers ──────────────────────────────────────────────────────────────────
normalize_pid <- function(x) trimws(as.character(x))

classify_timing <- function(days_to_addon, last_haic_day) {
  if (is.na(days_to_addon) || days_to_addon >= 9999) return("Never")
  if (days_to_addon <= BRIDGING_END_DAYS) return("Early concurrent (<=42d)")
  if (!is.na(last_haic_day) && days_to_addon > last_haic_day)
    return("Maintenance (> last HAIC)")
  return("During-HAIC (>42d, <= last HAIC)")
}

label_curative <- function(has_resect, has_ablate, has_tx) {
  # Priority: transplant > resection > ablation > none
  if (isTRUE(has_tx)) return("Liver transplant")
  if (isTRUE(has_resect)) return("Resection")
  if (isTRUE(has_ablate)) return("Ablation")
  return("None")
}

# Arm-allegiance derivation from CCW logic in tte_IT_R_two_cohorts.R
# Inputs are per-patient: triggered (0/1), trigger_day (NA if not), actual_treatment_day
classify_arm <- function(triggered, trigger_day, actual_day) {
  triggered <- as.integer(triggered)
  if (is.na(actual_day)) actual_day <- 9999
  early <- actual_day <= EARLY_GRACE_DAYS
  if (early) return("Early Combination")
  # Adaptive arm compliant: never triggered & never added; OR triggered & added within grace after trigger
  if (triggered == 0 && actual_day >= 9999) return("Adaptive On Demand")
  if (triggered == 1 && !is.na(trigger_day) && actual_day >= trigger_day &&
      actual_day <= trigger_day + DYNAMIC_GRACE_DAYS)
    return("Adaptive On Demand")
  # Untriggered but added; or triggered but added too early/too late
  return("Non-adherent")
}

# ── Load shared data once ────────────────────────────────────────────────────
cat("\nLoading shared event/summary tables...\n")
swimmer_all  <- read.csv(swimmer_path,  stringsAsFactors = FALSE, check.names = FALSE,
                         fileEncoding = "UTF-8")
summary_all  <- read.csv(summary_path,  stringsAsFactors = FALSE, check.names = FALSE,
                         fileEncoding = "UTF-8")
baseline_all <- read.csv(baseline_path, stringsAsFactors = FALSE, check.names = FALSE)
cat(sprintf("  swimmer events: %d rows\n", nrow(swimmer_all)))
cat(sprintf("  patient summary: %d rows\n", nrow(summary_all)))
cat(sprintf("  baseline: %d rows\n", nrow(baseline_all)))

swimmer_all$patient_id <- normalize_pid(swimmer_all$patient_id)
summary_all$patient_id <- normalize_pid(summary_all$patient_id)
baseline_all$patient_id <- normalize_pid(baseline_all$patient_id)

# Pre-compute derived event categories: bundle HAIC and TACE
swimmer_all <- swimmer_all %>%
  mutate(
    treatment_category = case_when(
      treatment_category %in% c("HAIC", "TACE", "HAIC+TACE") ~ "HAIC/TACE",
      TRUE ~ treatment_category
    )
  )

# ── Per-cohort pipeline ──────────────────────────────────────────────────────
run_cohort_pathways <- function(cfg) {
  out_dir <- file.path(base_out_dir, cfg$out_subdir)
  if (!dir.exists(out_dir))
    stop(sprintf("Cohort output dir not found: %s", out_dir))

  cat("\n######################################################################\n")
  cat(sprintf("## Cohort: %s — %s\n", cfg$out_subdir, cfg$label))
  cat("######################################################################\n")

  ids_csv     <- file.path(out_dir, sprintf("tte_cohort_%s_IT_ids.csv", cfg$out_subdir))
  trig_csv    <- file.path(out_dir, "R_trigger_table.csv")
  clone_csv   <- file.path(out_dir, "R_clone_dataset.csv")
  for (p in c(ids_csv, trig_csv, clone_csv))
    if (!file.exists(p)) stop("Missing input: ", p)

  ids   <- normalize_pid(read.csv(ids_csv,  stringsAsFactors = FALSE)$patient_id)
  trig  <- read.csv(trig_csv,  stringsAsFactors = FALSE) %>%
    mutate(patient_id = normalize_pid(patient_id))
  clone <- read.csv(clone_csv, stringsAsFactors = FALSE) %>%
    mutate(patient_id = normalize_pid(patient_id))

  cat(sprintf("  Cohort patients: %d\n", length(ids)))

  # Subset baseline / events to cohort
  bl <- baseline_all %>% filter(patient_id %in% ids)
  ev <- swimmer_all  %>% filter(patient_id %in% ids)
  sm <- summary_all  %>% filter(patient_id %in% ids)

  # Per-patient timing facts
  bl <- bl %>% mutate(
    first_haic_date = as.Date(first_haic_date),
    first_immune_date = suppressWarnings(as.Date(first_immune_date)),
    days_to_immune  = suppressWarnings(as.numeric(days_haic_to_immune)),
    death_status_int = as.integer(
      death_status == 1 |
      tolower(trimws(as.character(death_status))) %in% c("yes", "1", "true")),
    os_months_num = suppressWarnings(as.numeric(os_months))
  )
  if ("days_haic_to_target" %in% names(bl)) {
    bl$days_to_target <- suppressWarnings(as.numeric(bl$days_haic_to_target))
  } else {
    bl$days_to_target <- NA_real_
  }
  bl$has_immune <- as.integer(
    bl$has_immunotherapy == 1 | tolower(as.character(bl$has_immunotherapy)) == "yes")
  bl$has_target <- as.integer(
    bl$has_target_therapy == 1 | tolower(as.character(bl$has_target_therapy)) == "yes")

  # Last HAIC day per patient (in days from first HAIC), from event table
  haic_events <- ev %>%
    filter(treatment_category == "HAIC/TACE") %>%
    group_by(patient_id) %>%
    summarise(
      last_haic_day = max(time_days, na.rm = TRUE),
      n_haic_events = n(),
      .groups = "drop"
    )

  # Curative events (date of first such)
  cur_events <- ev %>%
    filter(treatment_category %in% c("Resection", "Ablation", "Liver Transplant")) %>%
    group_by(patient_id, treatment_category) %>%
    summarise(first_day = min(time_days, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = treatment_category, values_from = first_day,
                values_fill = NA_real_)
  for (col in c("Resection", "Ablation", "Liver Transplant"))
    if (!col %in% names(cur_events)) cur_events[[col]] <- NA_real_

  # Build per-patient classification table
  pat <- bl %>%
    transmute(patient_id, first_haic_date, has_immune, has_target,
              days_to_immune, days_to_target,
              os_days, os_months = os_months_num,
              death_status = death_status_int,
              first_immune_date = if ("first_immune_date" %in% names(bl))
                first_immune_date else as.Date(NA)) %>%
    left_join(haic_events, by = "patient_id") %>%
    left_join(cur_events,  by = "patient_id") %>%
    left_join(trig %>% select(patient_id, trigger_day, trigger_stage, eff_triggered),
              by = "patient_id")

  # Add actual_treatment_day per cohort treatment_mode
  if (cfg$regimen_target == "I+T") {
    pat <- pat %>% mutate(
      treatment_added = as.integer(has_immune == 1 & has_target == 1),
      actual_day = ifelse(treatment_added == 1,
                          pmax(coalesce(days_to_immune, 9999),
                               coalesce(days_to_target, 9999)),
                          9999)
    )
  } else {
    pat <- pat %>% mutate(
      treatment_added = as.integer(has_immune == 1),
      actual_day = ifelse(treatment_added == 1,
                          coalesce(days_to_immune, 9999),
                          9999)
    )
  }

  pat <- pat %>% mutate(
    timing       = mapply(classify_timing, actual_day, last_haic_day),
    curative_outcome = mapply(label_curative,
                              !is.na(.data[["Resection"]]),
                              !is.na(.data[["Ablation"]]),
                              !is.na(.data[["Liver Transplant"]])),
    arm_allegiance   = mapply(classify_arm, eff_triggered, trigger_day, actual_day)
  ) %>%
    mutate(
      timing           = factor(timing,  levels = names(COL_TIMING)),
      # First level = top of stratum stack in ggalluvial; put curative
      # endpoints at the top so they are visually emphasised.
      curative_outcome = factor(curative_outcome,
                                 levels = c("Resection", "Ablation",
                                            "Liver transplant", "None")),
      arm_allegiance   = factor(arm_allegiance,   levels = names(COL_ARM))
    )

  # Persist classification table
  out_class <- pat %>%
    select(patient_id, arm_allegiance, timing, curative_outcome,
           actual_day, last_haic_day, eff_triggered, trigger_day,
           os_months, death_status,
           any_of(c("Resection", "Ablation", "Liver Transplant")))
  write.csv(out_class, file.path(out_dir, "pathway_classification_table.csv"),
            row.names = FALSE)
  cat("  Wrote pathway_classification_table.csv\n")

  # Print breakdown
  cat("\n  Arm allegiance x Timing:\n")
  print(table(pat$arm_allegiance, pat$timing))
  cat("\n  Arm allegiance x Curative outcome:\n")
  print(table(pat$arm_allegiance, pat$curative_outcome))

  # ── Sankey: full cohort (3 arms) and adherent-only (Early + Adaptive) ────
  build_sankey(pat, cfg, out_dir, suffix = "",
               subtitle_extra = "")
  pat_adh <- pat %>%
    filter(arm_allegiance %in% c("Early Combination", "Adaptive On Demand")) %>%
    mutate(arm_allegiance = droplevels(arm_allegiance))
  build_sankey(pat_adh, cfg, out_dir, suffix = "_adherent_only",
               subtitle_extra = " — Non-adherent excluded")

  invisible(pat)
}

# ── Sankey builder (curative-conversion focus) ──────────────────────────────
# Design rationale (publication-grade):
#   • Ribbons are coloured by *curative outcome* so resection/ablation/transplant
#     paths visually pop while "no curative event" fades into a neutral tone.
#   • Arm strata carry a thin colour tab (COL_ARM) on their left edge, giving
#     each TTE arm an instantly recognisable visual identity.
#   • Left/right margin labels use redistributed y-positions plus leader lines
#     so they never collide, even when stratum heights differ by >100×
#     (e.g. n=3 transplant beside n=1685 None).
#   • Subtitle reports overall conversion rate and a chi-square test of
#     arm × curative-vs-none, anchoring the figure with a statistical claim.
#   • Output: 12 × 7.2 in, cairo_pdf for proper Unicode (·, ×, χ²),
#     PNG at 600 DPI for ≥print-grade reproduction.
build_sankey <- function(pat, cfg, out_dir, suffix = "", subtitle_extra = "") {
  N_total <- nrow(pat)
  if (N_total == 0) {
    cat("  build_sankey: empty input, skipping ", suffix, "\n", sep = "")
    return(invisible(NULL))
  }

  # ── Per-arm summary statistics ──────────────────────────────────────────
  arm_summary <- pat %>%
    group_by(arm_allegiance, .drop = FALSE) %>%
    summarise(
      n_arm        = n(),
      n_resection  = sum(curative_outcome == "Resection"),
      n_ablation   = sum(curative_outcome == "Ablation"),
      n_transplant = sum(curative_outcome == "Liver transplant"),
      n_curative   = sum(curative_outcome != "None"),
      pct_arm      = ifelse(N_total > 0, n() / N_total * 100, 0),
      pct_curative = ifelse(n() > 0, n_curative / n() * 100, 0),
      .groups      = "drop"
    ) %>%
    filter(n_arm > 0)

  # Pearson chi-square: arm × (curative vs no curative event)
  pval_text <- tryCatch({
    tab <- with(pat, table(arm_allegiance, curative_outcome != "None"))
    tab <- tab[rowSums(tab) > 0, , drop = FALSE]
    if (nrow(tab) >= 2 && ncol(tab) >= 2) {
      pv <- suppressWarnings(stats::chisq.test(tab)$p.value)
      if (is.na(pv))      ""
      else if (pv < 0.001) "P < 0.001"
      else                 sprintf("P = %.3f", pv)
    } else ""
  }, error = function(e) "")

  total_curative     <- sum(arm_summary$n_curative)
  total_curative_pct <- if (N_total > 0) total_curative / N_total * 100 else 0

  alluv_df <- pat %>%
    count(arm_allegiance, timing, curative_outcome, name = "Freq") %>%
    filter(Freq > 0) %>%
    mutate(prop = Freq / N_total)

  # Stratum geometry per axis (ggalluvial stacks first factor level at top)
  axis_specs <- list(
    list(var = "arm_allegiance",   x = 1),
    list(var = "timing",           x = 2),
    list(var = "curative_outcome", x = 3)
  )
  # IMPORTANT: arrange(desc(stratum)) must operate on the FACTOR (not character)
  # to match ggalluvial's stacking order (first factor level → top of stack).
  # Convert to character only AFTER computing y-positions so bind_rows merges
  # cleanly across axes with different factor levels.
  strat_pos <- bind_rows(lapply(axis_specs, function(a) {
    pat %>%
      count(stratum = .data[[a$var]], name = "count") %>%
      filter(count > 0) %>%
      arrange(desc(stratum)) %>%
      mutate(
        cum_top = cumsum(count),
        y_top   = cum_top,
        y_bot   = cum_top - count,
        y       = cum_top - count / 2,
        x       = a$x,
        stratum = as.character(stratum)
      )
  }))

  # Display name helpers — wrap long names so they fit inside a 0.28-wide
  # stratum at the inline label font size.
  display_outcome <- function(s) ifelse(s == "None", "No curative\nevent", s)
  ARM_DISPLAY <- c(
    "Early Combination"   = "Early\nCombination",
    "Adaptive On Demand"  = "Adaptive\nOn Demand",
    "Non-adherent"        = "Non-adherent"
  )

  # Inline-text colours for arm and outcome strata (parallels COL_TIMING_TEXT
  # for the middle column): white on saturated fills, dark on pale fills.
  COL_ARM_TEXT <- c(
    "Early Combination"   = "white",
    "Adaptive On Demand"  = "white",
    "Non-adherent"        = "white"
  )
  COL_OUTCOME_TEXT_INLINE <- c(
    "Resection"        = "white",
    "Ablation"         = "white",
    "Liver transplant" = "white",
    "None"             = "grey20"
  )

  # ── Left column (arms): inline label on the coloured stratum ────────────
  # Format: <wrapped name> / n=N (P%) / Curative P%   — short enough to fit
  # inside the 0.28-wide stratum without overflow.
  arm_lbl <- strat_pos %>%
    filter(x == 1) %>%
    left_join(arm_summary %>%
                mutate(arm_allegiance = as.character(arm_allegiance)) %>%
                select(arm_allegiance, n_arm, n_curative, pct_curative),
              by = c("stratum" = "arm_allegiance")) %>%
    mutate(
      label = sprintf("%s\nn = %d (%.0f%%)\nCurative %.1f%%",
                      ARM_DISPLAY[stratum],
                      n_arm, n_arm / N_total * 100,
                      pct_curative),
      text_color = unname(COL_ARM_TEXT[stratum])
    )

  # ── Middle column (timing): inline label, defined earlier ──────────────
  timing_lbl <- strat_pos %>%
    filter(x == 2) %>%
    mutate(label = sprintf("%s\n%d  (%.0f%%)",
                           TIMING_SHORT[stratum], count,
                           count / N_total * 100))

  # Right column (outcomes): inline label only on strata wide enough to hold
  # readable text. Tiny strata (e.g. Liver transplant, n=3) are colour-coded
  # in the ribbons + bottom legend, so no extra callout is needed.
  TINY_FRAC <- 0.020
  outcome_inline_lbl <- strat_pos %>%
    filter(x == 3, count >= TINY_FRAC * N_total) %>%
    mutate(
      label      = sprintf("%s\n%d  (%.1f%%)",
                            display_outcome(stratum), count,
                            count / N_total * 100),
      text_color = unname(COL_OUTCOME_TEXT_INLINE[stratum])
    )

  axis_titles <- data.frame(
    x     = c(1, 2, 3),
    label = c("TTE arm allegiance", "Add-on initiation window",
              "Curative-intent endpoint")
  )

  # ── Title / subtitle ────────────────────────────────────────────────────
  display_label <- gsub("\\s*\\(matched_\\d+\\)\\s*", " (PSM-matched)", cfg$label)
  display_label <- trimws(display_label)
  is_adherent_only <- nzchar(suffix) && grepl("adherent", suffix)
  title_text <- if (is_adherent_only)
      "Real-world treatment pathways and curative-intent conversion (adherent arms)"
    else
      "Real-world treatment pathways and curative-intent conversion"
  # Use ASCII-safe punctuation only (· × χ² fail on systems without cairo_pdf
  # and would silently strip characters from the PDF subtitle)
  subtitle_text <- sprintf(
    "%s  |  n = %d patients  |  overall curative conversion %d / %d (%.1f%%)%s",
    display_label, N_total, total_curative, N_total, total_curative_pct,
    if (nzchar(pval_text))
      sprintf("  |  arm-by-outcome %s (Pearson chi-square)", pval_text)
    else "")

  # ── Combined fill scale ─────────────────────────────────────────────────
  # Every stratum (arm, timing, outcome) and every ribbon (curative outcome)
  # resolves through one scale_fill_manual call. `breaks` restricts the legend
  # to just the curative-outcome entries so the legend stays focused.
  fill_values <- c(COL_CURATIVE, COL_ARM, COL_TIMING)
  fill_breaks <- names(COL_CURATIVE)
  fill_labels <- vapply(fill_breaks, display_outcome, character(1))

  # ── Build the plot ──────────────────────────────────────────────────────
  p <- ggplot(alluv_df,
              aes(axis1 = arm_allegiance,
                  axis2 = timing,
                  axis3 = curative_outcome,
                  y     = Freq)) +
    # 1. Ribbons coloured by final curative outcome — variable alpha so
    #    curative paths visually pop while "No curative event" recedes
    geom_alluvium(aes(fill = curative_outcome, alpha = curative_outcome),
                  width = 0.28, knot.pos = 0.40, colour = NA) +
    # 2. Stratum nodes — each stratum filled with its own category colour
    #    so all three columns (arm / timing / outcome) are immediately
    #    distinguishable. geom_rect keyed off `stratum` resolves through
    #    the combined fill scale below (COL_ARM ∪ COL_TIMING ∪ COL_CURATIVE).
    geom_rect(data = strat_pos,
              aes(xmin = x - 0.140, xmax = x + 0.140,
                  ymin = y_bot,     ymax = y_top,
                  fill = stratum),
              inherit.aes = FALSE,
              colour = "grey25", linewidth = 0.45,
              show.legend = FALSE) +
    # 4. Left column (arms): inline 3-line label on the coloured stratum.
    #    Text colour passed directly (per row) so each arm gets the right
    #    contrast against its fill (white on dark blue/purple/grey).
    geom_text(data = arm_lbl,
              aes(x = x, y = y, label = label),
              inherit.aes = FALSE,
              colour = arm_lbl$text_color,
              size = 2.75, lineheight = 1.15, fontface = "bold") +
    # 5. Middle column (timing): inline label, text colour adapts per stratum
    #    (white on dark green, dark on mint/grey).
    geom_text(data = timing_lbl,
              aes(x = x, y = y, label = label, colour = stratum),
              inherit.aes = FALSE,
              size = 2.80, lineheight = 1.10, fontface = "bold",
              show.legend = FALSE) +
    # 6. Right column (outcomes): inline labels only — tiny strata are
    #    covered by the ribbon colour + bottom legend.
    geom_text(data = outcome_inline_lbl,
              aes(x = x, y = y, label = label),
              inherit.aes = FALSE,
              colour = outcome_inline_lbl$text_color,
              size = 2.85, lineheight = 1.10, fontface = "bold") +
    # 7. Column headers — bold text on white (no grey backdrop)
    geom_text(data = axis_titles,
              aes(x = x, y = N_total * 1.055, label = label),
              inherit.aes = FALSE,
              fontface = "bold", size = 3.7, colour = "grey15") +
    # ── Scales ────────────────────────────────────────────────────────────
    # All labels are inline within strata; tight x-limits keep the figure
    # compact and free of margin whitespace.
    scale_x_continuous(limits = c(0.55, 3.45), expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
    scale_fill_manual(values = fill_values,
                      breaks = fill_breaks,
                      labels = fill_labels,
                      name   = "Curative-intent endpoint",
                      drop   = FALSE,
                      guide  = guide_legend(
                        override.aes = list(alpha = 0.85, colour = NA),
                        keywidth  = grid::unit(0.95, "cm"),
                        keyheight = grid::unit(0.45, "cm"),
                        nrow      = 1)) +
    scale_alpha_manual(values = c("Resection"        = 0.92,
                                  "Ablation"         = 0.92,
                                  "Liver transplant" = 0.92,
                                  "None"             = 0.32),
                       guide = "none") +
    scale_colour_manual(values = c(COL_CURATIVE_TEXT, COL_TIMING_TEXT),
                        guide = "none") +
    labs(title = title_text, subtitle = subtitle_text) +
    coord_cartesian(clip = "off") +
    theme_void(base_size = 11) +
    theme(
      plot.background    = element_rect(fill = "white", colour = NA),
      panel.background   = element_rect(fill = "white", colour = NA),
      legend.position    = "bottom",
      legend.title       = element_text(face = "bold", size = 9.5,
                                         margin = margin(r = 6)),
      legend.text        = element_text(size = 9),
      legend.box.margin  = margin(t = 4),
      legend.margin      = margin(t = 0, b = 0),
      plot.title         = element_text(face = "bold", size = 14,
                                         colour = "grey10",
                                         margin = margin(b = 3)),
      plot.subtitle      = element_text(size = 10, colour = "grey25",
                                         margin = margin(b = 12),
                                         lineheight = 1.20),
      plot.margin        = margin(14, 14, 8, 14)
    )

  # ── Companion panel: per-arm curative-modality breakdown ────────────────
  # Horizontal stacked bar — one row per arm, segments = Resection / Ablation /
  # Liver transplant. Each segment annotated with N (% of arm). Bar end carries
  # the total curative N / N_arm (%). Per-modality between-arm chi-square in
  # the subtitle anchors the comparison statistically.
  p_bar <- build_curative_modality_bar(arm_summary, cfg)

  # ── Combined publication figure (Sankey + companion bar) ────────────────
  # patchwork stacks the two panels with a 2.6:1 height ratio; legends are
  # NOT collected so each panel keeps its own focused legend.
  p_combined <- (p / p_bar) +
    patchwork::plot_layout(heights = c(2.6, 1)) &
    theme(plot.background = element_rect(fill = "white", colour = NA))

  # ── Save outputs ────────────────────────────────────────────────────────
  pdf_device <- tryCatch({
    grDevices::cairo_pdf(tempfile(fileext = ".pdf"), width = 1, height = 1)
    grDevices::dev.off()
    cairo_pdf
  }, error = function(e) "pdf",
     warning = function(w) "pdf")

  # 1) Sankey alone
  ggsave(file.path(out_dir, sprintf("pathway_A_sankey%s.pdf", suffix)),
         p, width = 12, height = 7.2, device = pdf_device)
  ggsave(file.path(out_dir, sprintf("pathway_A_sankey%s.png", suffix)),
         p, width = 12, height = 7.2, dpi = 600)

  # 2) Companion bar alone
  ggsave(file.path(out_dir, sprintf("pathway_B_modality_bar%s.pdf", suffix)),
         p_bar, width = 11, height = 3.6, device = pdf_device)
  ggsave(file.path(out_dir, sprintf("pathway_B_modality_bar%s.png", suffix)),
         p_bar, width = 11, height = 3.6, dpi = 600)

  # 3) Combined publication figure (recommended for journal submission)
  ggsave(file.path(out_dir, sprintf("pathway_combined%s.pdf", suffix)),
         p_combined, width = 12, height = 9.4, device = pdf_device)
  ggsave(file.path(out_dir, sprintf("pathway_combined%s.png", suffix)),
         p_combined, width = 12, height = 9.4, dpi = 600)

  write.csv(alluv_df, file.path(out_dir,
                                 sprintf("pathway_pathway_counts%s.csv",
                                         suffix)),
            row.names = FALSE)
  write.csv(arm_summary,
            file.path(out_dir,
                      sprintf("pathway_curative_by_arm%s.csv", suffix)),
            row.names = FALSE)

  cat(sprintf("  Wrote pathway_{A_sankey,B_modality_bar,combined}%s.{pdf,png} (n=%d)\n",
              suffix, N_total))
  cat(sprintf("    curative by arm: %s\n",
              paste(sprintf("%s=%.1f%%",
                            substr(arm_summary$arm_allegiance, 1, 5),
                            arm_summary$pct_curative),
                    collapse = ", ")))
  invisible(p_combined)
}

# ── Companion bar: per-arm × per-modality conversion rate ───────────────────
build_curative_modality_bar <- function(arm_summary, cfg) {
  # Long format: one row per (arm, modality)
  bar_data <- arm_summary %>%
    select(arm_allegiance, n_arm, n_resection, n_ablation, n_transplant) %>%
    tidyr::pivot_longer(
      cols       = c(n_resection, n_ablation, n_transplant),
      names_to   = "modality_var",
      values_to  = "n_mod"
    ) %>%
    mutate(
      modality = factor(c(
        n_resection  = "Resection",
        n_ablation   = "Ablation",
        n_transplant = "Liver transplant"
      )[modality_var],
      levels = c("Resection", "Ablation", "Liver transplant")),
      pct_mod = ifelse(n_arm > 0, n_mod / n_arm * 100, 0),
      arm_allegiance = forcats::fct_rev(as.factor(arm_allegiance))
    )

  end_lbl <- arm_summary %>%
    mutate(arm_allegiance = forcats::fct_rev(as.factor(arm_allegiance)),
           end_label = sprintf("Total %d / %d  (%.1f%%)",
                                n_curative, n_arm, pct_curative))

  # Per-modality between-arm chi-square (Fisher fall-back if any expected < 5)
  test_one <- function(col) {
    n_with    <- arm_summary[[col]]
    n_without <- arm_summary$n_arm - n_with
    tab       <- cbind(n_with, n_without)
    rownames(tab) <- as.character(arm_summary$arm_allegiance)
    tab <- tab[rowSums(tab) > 0, , drop = FALSE]
    if (nrow(tab) < 2) return(NA_real_)
    expected_ok <- tryCatch(
      all(suppressWarnings(stats::chisq.test(tab)$expected) >= 5),
      error = function(e) FALSE)
    pv <- if (expected_ok) {
      suppressWarnings(stats::chisq.test(tab)$p.value)
    } else {
      tryCatch(suppressWarnings(stats::fisher.test(tab,
                                                   simulate.p.value = TRUE,
                                                   B = 5000)$p.value),
               error = function(e) NA_real_)
    }
    pv
  }
  pv_res <- test_one("n_resection")
  pv_abl <- test_one("n_ablation")
  pv_tx  <- test_one("n_transplant")

  fmt_p <- function(p, label) {
    if (is.na(p))     return(sprintf("%s: NA", label))
    if (p < 0.001)    return(sprintf("%s: P < 0.001", label))
    sprintf("%s: P = %.3f", label, p)
  }
  chi_text <- paste(c(fmt_p(pv_res, "Resection"),
                       fmt_p(pv_abl, "Ablation"),
                       fmt_p(pv_tx,  "Liver transplant")),
                     collapse = "    |    ")

  xmax <- max(arm_summary$pct_curative, na.rm = TRUE) * 1.42

  display_label <- gsub("\\s*\\(matched_\\d+\\)\\s*", " (PSM-matched)", cfg$label)
  display_label <- trimws(display_label)

  ggplot(bar_data,
         aes(y = arm_allegiance, x = pct_mod, fill = modality)) +
    # Stacked horizontal bars — one row per arm, segments by modality.
    # reverse = TRUE so Resection (first factor level) appears at the LEFT,
    # which matches natural reading order and the Sankey legend.
    geom_col(width = 0.62, colour = "white", linewidth = 0.6,
             position = position_stack(reverse = TRUE)) +
    # In-segment N (%) labels — only when segment is wide enough
    geom_text(aes(label = ifelse(pct_mod >= 1.5,
                                  sprintf("%d  (%.1f%%)", n_mod, pct_mod),
                                  "")),
              position = position_stack(vjust = 0.5, reverse = TRUE),
              size = 3.10, colour = "white", fontface = "bold",
              show.legend = FALSE) +
    # Tiny segments → just N inside (percentage redundant for n=1, n=2).
    # Must use the same x/fill aesthetics as geom_col so position_stack lines
    # the labels up with the right segment.
    geom_text(aes(label = ifelse(n_mod > 0 & pct_mod < 1.5,
                                  sprintf("%d", n_mod), "")),
              position = position_stack(vjust = 0.5, reverse = TRUE),
              size = 2.7, colour = "white", fontface = "bold",
              show.legend = FALSE) +
    # End-of-bar total
    geom_text(data = end_lbl,
              aes(y = arm_allegiance,
                  x = pct_curative + xmax * 0.012,
                  label = end_label),
              inherit.aes = FALSE,
              hjust = 0, size = 3.15,
              colour = "grey15", fontface = "bold") +
    scale_fill_manual(
      values = COL_CURATIVE[c("Resection", "Ablation", "Liver transplant")],
      name   = "Curative modality",
      guide  = guide_legend(nrow = 1,
                            keywidth  = grid::unit(0.85, "cm"),
                            keyheight = grid::unit(0.42, "cm"))) +
    scale_x_continuous(limits = c(0, xmax), expand = c(0, 0),
                       breaks = scales::pretty_breaks(6),
                       labels = function(x) paste0(x, "%")) +
    labs(
      x        = "Conversion rate (% of arm)",
      y        = NULL,
      title    = "Curative-intent conversion broken down by treatment modality",
      subtitle = sprintf(
        "%s    |    Between-arm test per modality:    %s",
        display_label, chi_text)) +
    theme_classic(base_size = 11) +
    theme(
      plot.background    = element_rect(fill = "white", colour = NA),
      plot.title         = element_text(face = "bold", size = 12,
                                         colour = "grey10",
                                         margin = margin(b = 2)),
      plot.subtitle      = element_text(size = 9, colour = "grey30",
                                         margin = margin(b = 10),
                                         lineheight = 1.20),
      axis.text.y        = element_text(face = "bold", size = 10.5,
                                         colour = "grey15"),
      axis.text.x        = element_text(size = 9, colour = "grey30"),
      axis.title.x       = element_text(size = 10, colour = "grey15",
                                         margin = margin(t = 6)),
      axis.line.y        = element_blank(),
      axis.ticks.y       = element_blank(),
      axis.line.x        = element_line(colour = "grey60", linewidth = 0.4),
      axis.ticks.x       = element_line(colour = "grey60", linewidth = 0.4),
      legend.position    = "bottom",
      legend.title       = element_text(face = "bold", size = 9.5),
      legend.text        = element_text(size = 9),
      legend.box.margin  = margin(t = 0),
      legend.margin      = margin(t = 0, b = 0),
      panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
      plot.margin        = margin(8, 18, 8, 18)
    )
}

# ── Driver ───────────────────────────────────────────────────────────────────
for (cn in names(COHORTS)) {
  res <- tryCatch(run_cohort_pathways(COHORTS[[cn]]),
                  error = function(e) {
                    cat(sprintf("\n[ERROR] cohort %s: %s\n", cn, conditionMessage(e)))
                    NULL
                  })
}

cat("\n======================================================================\n")
cat("Pathway visualization complete.\n")
cat(sprintf("Outputs under: %s\n", base_out_dir))
cat("======================================================================\n")
