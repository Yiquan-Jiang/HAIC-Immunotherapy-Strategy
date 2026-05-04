#!/usr/bin/env Rscript
# =============================================================================
# Companion to tte_pathway_visualization.R — alternative sample sets
# =============================================================================
# Purpose : Reuse only the timing × curative-endpoint half of the main pathway
#           pipeline, on four DIFFERENT samples requested by the PI:
#
#   Mode A — "Adaptive-On-Demand only" (filtered TTE cohorts)
#     • Restrict to patients whose TTE arm allegiance = "Adaptive On Demand".
#     • No sub-classification of the arm: column 1 = "Add-on timing",
#       column 2 = "Curative endpoint".  One figure per TTE cohort.
#
#   Mode B — "HAIC_then_I" and "HAIC_then_I+T" (7-group sequence labels)
#     • Two SEPARATE figures (one per group).  Same 2-column structure as
#       Mode A.  IDs are pulled from the 7-group sequence-label file then
#       restricted to the project inclusion cohort (HAIC_NO_TACE_4_TIDY
#       baseline, 4,234 patients) so per-group N matches the agreed
#       7-group spec: HAIC_then_I = 152, HAIC_then_I+T = 221.
#
# Output structure:
#   output/step3_tte/IT_RULES_R_two_cohorts/pathway_alt_samples/
#     ├── adaptive_only_cohort_3matched/
#     │     pathway_sankey.{pdf,png}
#     │     pathway_classification_table.csv
#     │     pathway_pathway_counts.csv
#     ├── adaptive_only_cohort_7group_psm02/  (same files)
#     ├── group7_haic_then_I/                  (same files)
#     └── group7_haic_then_IT/                 (same files)
# =============================================================================

Sys.setenv(LANG = "en_US.UTF-8")
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(ggalluvial)
  library(ggrepel)
})

# ── CLI / path resolution ───────────────────────────────────────────────────
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
base_tte_dir <- file.path(project_root, "output", "step3_tte",
                          "IT_RULES_R_two_cohorts")
if (!dir.exists(base_tte_dir))
  stop("TTE output dir not found — run tte_IT_R_two_cohorts.R first: ",
       base_tte_dir)

out_root <- file.path(base_tte_dir, "pathway_alt_samples")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

events_dir <- file.path(
  normalizePath(file.path(data_dir, ".."), winslash = "/", mustWork = TRUE),
  "HAIC_NO_TACE_4_TIDY", "update_group_7", "data"
)
if (!dir.exists(events_dir)) {
  events_dir_alt <- file.path(data_dir, "HAIC_NO_TACE_4_TIDY",
                              "update_group_7", "data")
  if (dir.exists(events_dir_alt)) events_dir <- events_dir_alt
}
swimmer_path     <- file.path(events_dir, "00_swimmer_plot_events.csv")
summary_path     <- file.path(events_dir, "00_patient_treatment_summary.csv")
seq_labels_path  <- file.path(events_dir, "patient_treatment_sequence_labels.csv")
for (p in c(swimmer_path, summary_path, seq_labels_path))
  if (!file.exists(p)) stop("Missing input: ", p)

baseline_path <- file.path(data_dir, "HAIC_NO_TACE_4_TIDY_baseline.csv")
if (!file.exists(baseline_path)) stop("Missing baseline: ", baseline_path)

cat("======================================================================\n")
cat("TTE pathway visualization — ALTERNATIVE samples (simplified 2-column)\n")
cat("======================================================================\n")
cat("data_dir   : ", data_dir,    "\n")
cat("base_tte   : ", base_tte_dir,"\n")
cat("events_dir : ", events_dir,  "\n")
cat("out_root   : ", out_root,    "\n")

# ── Constants ───────────────────────────────────────────────────────────────
EARLY_GRACE_DAYS    <- 14
DYNAMIC_GRACE_DAYS  <- 90
# NB: the original main script split "during HAIC" into two windows by a
# 42-day bridging cutoff. The PI requested a unified definition for these
# alt-sample figures because:
#   • the 7-group sequence-label cutoff is "≥0.1 month" (~3 days),
#     which is incompatible with 42 d
#   • Adaptive-On-Demand patients are by definition triggered LATE, so the
#     14–42 d window is essentially empty for them
# We therefore collapse "Early concurrent (≤42d)" into "During HAIC".

# ── Colour palettes ─────────────────────────────────────────────────────────
COL_TIMING <- c(
  "During HAIC (<= last HAIC)" = "#1B4332",
  "Maintenance (> last HAIC)"  = "#B7E4C7",
  "Never"                      = "#E5E5E5"
)
COL_CURATIVE <- c(
  "Resection"        = "#BC3C29",
  "Ablation"         = "#E18727",
  "Liver transplant" = "#0072B5",
  "None"             = "#DCDCDC"
)
COL_TIMING_TEXT <- c(
  "During HAIC (<= last HAIC)" = "white",
  "Maintenance (> last HAIC)"  = "grey15",
  "Never"                      = "grey20"
)
COL_OUTCOME_TEXT <- c(
  "Resection"        = "white",
  "Ablation"         = "white",
  "Liver transplant" = "white",
  "None"             = "grey20"
)
TIMING_SHORT <- c(
  "During HAIC (<= last HAIC)" = "During HAIC",
  "Maintenance (> last HAIC)"  = "Maintenance",
  "Never"                      = "Never added"
)

# Map for re-coding the old 4-category timing (still found in the
# pre-computed pathway_classification_table.csv from the main script) into
# the unified 3-category scheme used here.
TIMING_RECODE <- c(
  "Early concurrent (<=42d)"         = "During HAIC (<= last HAIC)",
  "During-HAIC (>42d, <= last HAIC)" = "During HAIC (<= last HAIC)",
  "Maintenance (> last HAIC)"        = "Maintenance (> last HAIC)",
  "Never"                            = "Never"
)

# ── Helpers (verbatim from main script) ─────────────────────────────────────
normalize_pid <- function(x) trimws(as.character(x))

classify_timing <- function(days_to_addon, last_haic_day) {
  if (is.na(days_to_addon) || days_to_addon >= 9999) return("Never")
  if (!is.na(last_haic_day) && days_to_addon > last_haic_day)
    return("Maintenance (> last HAIC)")
  return("During HAIC (<= last HAIC)")
}

label_curative <- function(has_resect, has_ablate, has_tx) {
  if (isTRUE(has_tx))     return("Liver transplant")
  if (isTRUE(has_resect)) return("Resection")
  if (isTRUE(has_ablate)) return("Ablation")
  return("None")
}

# ── Load shared data once ────────────────────────────────────────────────────
cat("\nLoading shared event/summary tables...\n")
swimmer_all  <- read.csv(swimmer_path,  stringsAsFactors = FALSE,
                         check.names = FALSE, fileEncoding = "UTF-8")
summary_all  <- read.csv(summary_path,  stringsAsFactors = FALSE,
                         check.names = FALSE, fileEncoding = "UTF-8")
seq_labels_all <- read.csv(seq_labels_path, stringsAsFactors = FALSE,
                           check.names = FALSE, fileEncoding = "UTF-8")
baseline_all <- read.csv(baseline_path, stringsAsFactors = FALSE,
                         check.names = FALSE)

cat(sprintf("  swimmer events : %d rows\n",  nrow(swimmer_all)))
cat(sprintf("  patient summary: %d rows\n",  nrow(summary_all)))
cat(sprintf("  seq labels     : %d rows\n",  nrow(seq_labels_all)))
cat(sprintf("  baseline       : %d rows\n",  nrow(baseline_all)))

swimmer_all$patient_id    <- normalize_pid(swimmer_all$patient_id)
summary_all$patient_id    <- normalize_pid(summary_all$patient_id)
seq_labels_all$patient_id <- normalize_pid(seq_labels_all$patient_id)
baseline_all$patient_id   <- normalize_pid(baseline_all$patient_id)

swimmer_all <- swimmer_all %>%
  mutate(
    treatment_category = case_when(
      treatment_category %in% c("HAIC", "TACE", "HAIC+TACE") ~ "HAIC/TACE",
      TRUE ~ treatment_category
    )
  )

# ── Per-patient feature builder ──────────────────────────────────────────────
# Source preference (covers ALL HAIC patients, including those outside the
# TTE-baseline 4234-row subset):
#   • days_haic_to_immune / days_haic_to_target / has_immune / has_target
#       <- 00_patient_treatment_summary.csv  (10,418 rows)
#   • last_haic_day, curative event days
#       <- 00_swimmer_plot_events.csv
build_patient_features <- function(ids, regimen_target = c("I", "I+T")) {
  regimen_target <- match.arg(regimen_target)

  ps <- summary_all %>% filter(patient_id %in% ids) %>%
    mutate(
      days_to_immune = suppressWarnings(as.numeric(days_haic_to_immune)),
      days_to_target = suppressWarnings(as.numeric(days_haic_to_target)),
      has_immune     = as.integer(immune_episodes > 0),
      has_target     = as.integer(target_episodes > 0)
    )

  ev <- swimmer_all %>% filter(patient_id %in% ids)

  haic_events <- ev %>%
    filter(treatment_category == "HAIC/TACE") %>%
    group_by(patient_id) %>%
    summarise(last_haic_day = max(time_days, na.rm = TRUE),
              n_haic_events = n(), .groups = "drop")

  cur_events <- ev %>%
    filter(treatment_category %in% c("Resection", "Ablation",
                                     "Liver Transplant")) %>%
    group_by(patient_id, treatment_category) %>%
    summarise(first_day = min(time_days, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = treatment_category, values_from = first_day,
                values_fill = NA_real_)
  for (col in c("Resection", "Ablation", "Liver Transplant"))
    if (!col %in% names(cur_events)) cur_events[[col]] <- NA_real_

  pat <- ps %>%
    transmute(patient_id, has_immune, has_target,
              days_to_immune, days_to_target) %>%
    left_join(haic_events, by = "patient_id") %>%
    left_join(cur_events,  by = "patient_id")

  # Optional OS/death info (for audit table; not used in the figure)
  bl <- baseline_all %>% filter(patient_id %in% ids) %>%
    transmute(
      patient_id,
      death_status = as.integer(
        death_status == 1 |
        tolower(trimws(as.character(death_status))) %in% c("yes", "1", "true")),
      os_months = suppressWarnings(as.numeric(os_months))
    )
  pat <- pat %>% left_join(bl, by = "patient_id")

  if (regimen_target == "I+T") {
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

  pat %>% mutate(
    timing = mapply(classify_timing, actual_day, last_haic_day),
    curative_outcome = mapply(label_curative,
                              !is.na(.data[["Resection"]]),
                              !is.na(.data[["Ablation"]]),
                              !is.na(.data[["Liver Transplant"]]))
  )
}

# ── Two-column Sankey (timing → curative endpoint) ──────────────────────────
# A single-cohort figure: minimal text, just the two columns + small caption.
build_sankey_2col <- function(pat, title_text, out_dir, suffix = "") {
  N_total <- nrow(pat)
  if (N_total == 0) {
    cat("  build_sankey_2col: empty input, skipping\n")
    return(invisible(NULL))
  }

  pat <- pat %>%
    mutate(
      timing           = factor(timing, levels = names(COL_TIMING)),
      curative_outcome = factor(curative_outcome,
                                 levels = c("Resection", "Ablation",
                                            "Liver transplant", "None"))
    )

  n_curative   <- sum(pat$curative_outcome != "None")
  pct_curative <- if (N_total > 0) n_curative / N_total * 100 else 0
  subtitle_text <- sprintf("n = %d   |   Curative conversion %d (%.1f%%)",
                           N_total, n_curative, pct_curative)

  alluv_df <- pat %>%
    count(timing, curative_outcome, name = "Freq") %>%
    filter(Freq > 0)

  axis_specs <- list(
    list(var = "timing",           x = 1),
    list(var = "curative_outcome", x = 2)
  )
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

  display_outcome    <- function(s) ifelse(s == "None", "No curative\nevent", s)
  display_outcome_1l <- function(s) ifelse(s == "None", "No curative event", s)

  # Inline labels: only on strata with enough vertical room (≥2.5%) so text
  # stays inside the box. Smaller strata get a callout label OUTSIDE the
  # column with a thin leader line, so every stratum (including 1-patient
  # slivers like "Early concurrent" in the alt-sample cohorts) is identified
  # — no mystery boxes.
  TINY_FRAC <- 0.025

  timing_inline <- strat_pos %>%
    filter(x == 1, count >= TINY_FRAC * N_total) %>%
    mutate(label = sprintf("%s\n%d (%.0f%%)",
                           TIMING_SHORT[stratum], count,
                           count / N_total * 100))

  timing_outside <- strat_pos %>%
    filter(x == 1, count < TINY_FRAC * N_total) %>%
    mutate(label = sprintf("%s — %d (%.1f%%)",
                           TIMING_SHORT[stratum], count,
                           count / N_total * 100))

  outcome_inline <- strat_pos %>%
    filter(x == 2, count >= TINY_FRAC * N_total) %>%
    mutate(
      label      = sprintf("%s\n%d (%.1f%%)",
                            display_outcome(stratum), count,
                            count / N_total * 100),
      text_color = unname(COL_OUTCOME_TEXT[stratum])
    )

  outcome_outside <- strat_pos %>%
    filter(x == 2, count < TINY_FRAC * N_total) %>%
    mutate(label = sprintf("%s — %d (%.1f%%)",
                           display_outcome_1l(stratum), count,
                           count / N_total * 100))

  axis_titles <- data.frame(
    x     = c(1, 2),
    label = c("Add-on timing", "Curative endpoint")
  )

  fill_values <- c(COL_CURATIVE, COL_TIMING)
  fill_breaks <- names(COL_CURATIVE)
  fill_labels <- vapply(fill_breaks, display_outcome, character(1))

  p <- ggplot(alluv_df,
              aes(axis1 = timing,
                  axis2 = curative_outcome,
                  y     = Freq)) +
    geom_alluvium(aes(fill = curative_outcome, alpha = curative_outcome),
                  width = 0.30, knot.pos = 0.40, colour = NA) +
    geom_rect(data = strat_pos,
              aes(xmin = x - 0.15, xmax = x + 0.15,
                  ymin = y_bot,    ymax = y_top,
                  fill = stratum),
              inherit.aes = FALSE,
              colour = "grey25", linewidth = 0.45,
              show.legend = FALSE) +
    geom_text(data = timing_inline,
              aes(x = x, y = y, label = label, colour = stratum),
              inherit.aes = FALSE,
              size = 3.0, lineheight = 1.10, fontface = "bold",
              show.legend = FALSE) +
    geom_text(data = outcome_inline,
              aes(x = x, y = y, label = label),
              inherit.aes = FALSE,
              colour = outcome_inline$text_color,
              size = 3.0, lineheight = 1.10, fontface = "bold") +
    # Outside callouts for tiny strata (column 1, LEFT). ggrepel pushes
    # labels apart vertically when several tiny strata cluster, drawing thin
    # leader lines back to the column.
    ggrepel::geom_text_repel(
      data = timing_outside,
      aes(x = 0.85, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 1, size = 2.6, colour = "grey25",
      fontface = "plain",
      direction      = "y",
      nudge_x        = -0.20,
      xlim           = c(NA, 0.78),
      segment.colour = "grey55",
      segment.size   = 0.3,
      box.padding    = 0.15,
      min.segment.length = 0,
      seed = 1) +
    # Outside callouts for tiny strata (column 2, RIGHT)
    ggrepel::geom_text_repel(
      data = outcome_outside,
      aes(x = 2.15, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 0, size = 2.6, colour = "grey25",
      fontface = "plain",
      direction      = "y",
      nudge_x        = 0.20,
      xlim           = c(2.22, NA),
      segment.colour = "grey55",
      segment.size   = 0.3,
      box.padding    = 0.15,
      min.segment.length = 0,
      seed = 1) +
    geom_text(data = axis_titles,
              aes(x = x, y = N_total * 1.06, label = label),
              inherit.aes = FALSE,
              fontface = "bold", size = 3.8, colour = "grey15") +
    # Slightly wider x-limits so callout labels fit without clipping
    scale_x_continuous(limits = c(0.30, 2.70), expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0.04, 0.10))) +
    scale_fill_manual(values = fill_values,
                      breaks = fill_breaks,
                      labels = fill_labels,
                      name   = NULL,
                      drop   = FALSE,
                      guide  = guide_legend(
                        override.aes = list(alpha = 0.85, colour = NA),
                        keywidth  = grid::unit(0.85, "cm"),
                        keyheight = grid::unit(0.40, "cm"),
                        nrow      = 1)) +
    scale_alpha_manual(values = c("Resection"        = 0.92,
                                  "Ablation"         = 0.92,
                                  "Liver transplant" = 0.92,
                                  "None"             = 0.32),
                       guide = "none") +
    scale_colour_manual(values = COL_TIMING_TEXT, guide = "none") +
    labs(title = title_text, subtitle = subtitle_text) +
    coord_cartesian(clip = "off") +
    theme_void(base_size = 11) +
    theme(
      plot.background    = element_rect(fill = "white", colour = NA),
      panel.background   = element_rect(fill = "white", colour = NA),
      legend.position    = "bottom",
      legend.text        = element_text(size = 9),
      legend.box.margin  = margin(t = 4),
      legend.margin      = margin(t = 0, b = 0),
      plot.title         = element_text(face = "bold", size = 13,
                                         colour = "grey10",
                                         margin = margin(b = 2)),
      plot.subtitle      = element_text(size = 10, colour = "grey30",
                                         margin = margin(b = 8)),
      plot.margin        = margin(12, 14, 8, 14)
    )

  pdf_device <- tryCatch({
    grDevices::cairo_pdf(tempfile(fileext = ".pdf"), width = 1, height = 1)
    grDevices::dev.off()
    cairo_pdf
  }, error = function(e) "pdf",
     warning = function(w) "pdf")

  ggsave(file.path(out_dir, sprintf("pathway_sankey%s.pdf", suffix)),
         p, width = 7.5, height = 5.6, device = pdf_device)
  ggsave(file.path(out_dir, sprintf("pathway_sankey%s.png", suffix)),
         p, width = 7.5, height = 5.6, dpi = 600)

  write.csv(alluv_df,
            file.path(out_dir,
                      sprintf("pathway_pathway_counts%s.csv", suffix)),
            row.names = FALSE)

  cat(sprintf("  Wrote pathway_sankey%s.{pdf,png} (n=%d, curative=%d %.1f%%)\n",
              suffix, N_total, n_curative, pct_curative))
  invisible(p)
}

# ── Driver helper: write classification CSV (audit) and call sankey ──────────
emit_figure <- function(pat, title_text, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  cat(sprintf("  n = %d\n", nrow(pat)))
  cat("  Timing × Curative outcome:\n")
  print(table(pat$timing, pat$curative_outcome))

  audit_cols <- intersect(
    c("patient_id", "timing", "curative_outcome",
      "actual_day", "last_haic_day", "treatment_added",
      "os_months", "death_status",
      "Resection", "Ablation", "Liver Transplant"),
    names(pat)
  )
  write.csv(pat[, audit_cols],
            file.path(out_dir, "pathway_classification_table.csv"),
            row.names = FALSE)

  build_sankey_2col(pat, title_text, out_dir)
}

# ============================================================================
# MODE A — Adaptive On Demand only (per TTE cohort)
# ============================================================================
COHORTS_TTE <- list(
  cohort_3matched = list(
    out_subdir     = "cohort_3matched",
    label          = "HAIC alone vs HAIC+I+T",
    regimen_target = "I+T"
  ),
  cohort_7group_psm02 = list(
    out_subdir     = "cohort_7group_psm02",
    label          = "HAIC alone vs HAIC+I",
    regimen_target = "I"
  )
)

run_adaptive_only <- function(cfg) {
  cat("\n## Mode A | Adaptive On Demand only — ", cfg$out_subdir, "\n", sep = "")

  src_dir <- file.path(base_tte_dir, cfg$out_subdir)
  pathway_csv <- file.path(src_dir, "pathway_classification_table.csv")
  if (!file.exists(pathway_csv))
    stop("Run tte_pathway_visualization.R first; missing: ", pathway_csv)

  pat_full <- read.csv(pathway_csv, stringsAsFactors = FALSE)
  pat_full$patient_id <- normalize_pid(pat_full$patient_id)

  # ggalluvial-friendly column names: rename "Liver Transplant" -> dot version
  if ("Liver.Transplant" %in% names(pat_full))
    names(pat_full)[names(pat_full) == "Liver.Transplant"] <- "Liver Transplant"

  pat_aod <- pat_full %>%
    filter(arm_allegiance == "Adaptive On Demand") %>%
    mutate(
      # Re-map the legacy 4-category timing into the unified 3-category scheme
      timing = unname(TIMING_RECODE[timing]),
      timing = factor(timing, levels = names(COL_TIMING)),
      curative_outcome = factor(curative_outcome,
                                 levels = c("Resection", "Ablation",
                                            "Liver transplant", "None"))
    )
  if (nrow(pat_aod) == 0) {
    cat("  No Adaptive On Demand patients — skipping\n")
    return(invisible(NULL))
  }

  out_dir <- file.path(out_root, sprintf("adaptive_only_%s", cfg$out_subdir))
  title_text <- sprintf("Adaptive On-Demand pathway (%s)", cfg$label)
  emit_figure(pat_aod, title_text, out_dir)
}

# ============================================================================
# MODE B — 7-group "HAIC_then_I" and "HAIC_then_I+T" — separate figures
# ============================================================================
run_group7_one <- function(group_name, regimen_target, title_text, out_subdir) {
  cat("\n## Mode B | 7-group ", group_name, "\n", sep = "")

  ids <- seq_labels_all %>%
    filter(main_group == group_name) %>%
    pull(patient_id) %>%
    normalize_pid()

  # Restrict to the project inclusion cohort ("first 4 HAIC without TACE",
  # 4,234 patients in HAIC_NO_TACE_4_TIDY_baseline.csv). Without this, ids
  # come from the full HAIC universe (10,391) and inflate the per-group N
  # vs. the agreed 7-group spec (e.g. HAIC_then_I 393 -> 152, I+T 704 -> 221).
  baseline_ids <- normalize_pid(baseline_all$patient_id)
  n_pre <- length(ids)
  ids <- intersect(ids, baseline_ids)
  cat(sprintf("  ids in seq_labels: %d  |  in HAIC_NO_TACE_4 baseline: %d\n",
              n_pre, length(ids)))

  pat <- build_patient_features(ids, regimen_target = regimen_target) %>%
    mutate(
      timing           = factor(timing, levels = names(COL_TIMING)),
      curative_outcome = factor(curative_outcome,
                                 levels = c("Resection", "Ablation",
                                            "Liver transplant", "None"))
    )

  out_dir <- file.path(out_root, out_subdir)
  emit_figure(pat, title_text, out_dir)
}

# ── Driver ───────────────────────────────────────────────────────────────────
for (cn in names(COHORTS_TTE)) {
  res <- tryCatch(run_adaptive_only(COHORTS_TTE[[cn]]),
                  error = function(e) {
                    cat(sprintf("\n[ERROR] mode-A cohort %s: %s\n",
                                cn, conditionMessage(e)))
                    NULL
                  })
}

res_b1 <- tryCatch(run_group7_one(
  group_name      = "HAIC_then_I",
  regimen_target  = "I",
  title_text      = "HAIC then Immunotherapy",
  out_subdir      = "group7_haic_then_I"),
  error = function(e) {
    cat(sprintf("\n[ERROR] mode-B HAIC_then_I: %s\n", conditionMessage(e)))
    NULL
  })

res_b2 <- tryCatch(run_group7_one(
  group_name      = "HAIC_then_I+T",
  regimen_target  = "I+T",
  title_text      = "HAIC then Immunotherapy + Targeted",
  out_subdir      = "group7_haic_then_IT"),
  error = function(e) {
    cat(sprintf("\n[ERROR] mode-B HAIC_then_I+T: %s\n", conditionMessage(e)))
    NULL
  })

cat("\n======================================================================\n")
cat("Alternative-sample pathway visualization complete.\n")
cat(sprintf("Outputs under: %s\n", out_root))
cat("======================================================================\n")
