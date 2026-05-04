#!/usr/bin/env Rscript
# Swimmer plot v3 — 8 treatment types, 60-month axis with break, NPG colors
# Each of 7 groups exported as separate PDF/PNG

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(tibble)
})
if (!requireNamespace("patchwork", quietly = TRUE))
  install.packages("patchwork", repos = "https://cloud.r-project.org")
library(patchwork)

BASE_DIR <- "/Users/yqj/Nutstore Files/我的坚果云/Liver_tumor_big_data/FIRST_LINE_HAIC_2025-12-30/HAIC_NO_TACE_4_TIDY/update_group_7"
DATA_DIR <- file.path(BASE_DIR, "data")
FIG_DIR  <- file.path(BASE_DIR, "figures", "swimmer_7groups")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

MAX_PER_GROUP <- 50L
RNG_SEED <- 42L
X_DISPLAY_MAX <- 60  # hard cap for display

GROUP_ORDER <- c(
  "HAIC_alone", "HAIC+I_concurrent", "HAIC_then_I", "HAIC+T_concurrent",
  "HAIC_then_T", "HAIC+I+T_concurrent", "HAIC_then_I+T"
)
GROUP_PRETTY <- c(
  "HAIC_alone"           = "HAIC alone",
  "HAIC+I_concurrent"    = "HAIC + I (concurrent)",
  "HAIC_then_I"          = "HAIC then I (sequential)",
  "HAIC+T_concurrent"    = "HAIC + T (concurrent)",
  "HAIC_then_T"          = "HAIC then T (sequential)",
  "HAIC+I+T_concurrent"  = "HAIC + I + T (concurrent)",
  "HAIC_then_I+T"        = "HAIC then I+T (sequential)"
)

# ── 8 treatment types — NPG palette ─────────────────────────────
TX_LEVELS <- c("HAIC", "Immunotherapy", "Targeted", "TACE", "HAIC+TACE",
               "Resection", "Ablation", "Transplant")
TX_COLORS <- c(
  "HAIC"          = "#E64B35",
  "Immunotherapy" = "#4DBBD5",
  "Targeted"      = "#3C5488",
  "TACE"          = "#E18727",
  "HAIC+TACE"     = "#DC0000",
  "Resection"     = "#00A087",
  "Ablation"      = "#8491B4",
  "Transplant"    = "#7E6148"
)
TX_SHAPES <- c(
  "HAIC"          = 15L,
  "Immunotherapy" = 16L,
  "Targeted"      = 17L,
  "TACE"          = 18L,
  "HAIC+TACE"     = 13L,
  "Resection"     = 8L,
  "Ablation"      = 10L,
  "Transplant"    = 12L
)

# Survival bar colors (muted, won't compete with event points)
LANE_ALIVE <- "#91D1C2"   # NPG light teal
LANE_DEAD  <- "#E64B35"   # NPG red-orange (drawn at alpha=0.45)

safe_fname <- function(x) gsub("^_|_$", "", gsub("[^A-Za-z0-9]+", "_", x))

# ── Read data ────────────────────────────────────────────────────
bl <- read_csv(
  file.path(DATA_DIR, "analysis_ready.csv"),
  col_types = cols(.default = col_character()),
  show_col_types = FALSE
) %>%
  type_convert(col_types = cols()) %>%
  filter(main_group %in% GROUP_ORDER) %>%
  mutate(
    first_haic_date = as.Date(first_haic_date),
    death_status = case_when(
      death_status %in% c("Yes", "1") ~ 1L,
      death_status %in% c("No", "0")  ~ 0L,
      TRUE ~ NA_integer_
    ),
    os_months     = as.numeric(os_months),
    haic_episodes = suppressWarnings(as.integer(haic_episodes))
  ) %>%
  filter(!is.na(first_haic_date), !is.na(os_months), os_months >= 0)

events_raw <- read_csv(
  file.path(DATA_DIR, "00_swimmer_plot_events.csv"),
  col_types = cols(.default = col_character()),
  show_col_types = FALSE
) %>% type_convert(col_types = cols())

# ── Right legend (annotate) ──────────────────────────────────────
build_legend_right <- function() {
  y_top <- 0.98
  y_step <- 0.038
  sz_pt <- 3.2
  sz_tx <- 3.2
  sz_hdr <- 3.5
  tx_y <- y_top - seq_len(length(TX_LEVELS)) * y_step
  tx_labels <- c("HAIC", "Immunotherapy", "Targeted therapy",
                 "TACE", "HAIC + TACE", "Resection", "Ablation", "Transplant")

  p <- ggplot() + xlim(0, 1) + ylim(0.15, 1)

  p <- p + annotate("text", x = 0.03, y = y_top, label = "Treatment event",
                    hjust = 0, size = sz_hdr, fontface = "bold")
  for (i in seq_along(TX_LEVELS)) {
    p <- p +
      annotate("point", x = 0.06, y = tx_y[i],
               shape = TX_SHAPES[TX_LEVELS[i]], size = sz_pt,
               color = TX_COLORS[TX_LEVELS[i]]) +
      annotate("text", x = 0.11, y = tx_y[i],
               label = tx_labels[i], hjust = 0, size = sz_tx)
  }

  y_surv_hdr <- tx_y[length(TX_LEVELS)] - y_step * 1.2
  p <- p +
    annotate("text", x = 0.03, y = y_surv_hdr, label = "Survival duration",
             hjust = 0, size = sz_hdr, fontface = "bold") +
    annotate("segment", x = 0.06, xend = 0.18, y = y_surv_hdr - y_step,
             yend = y_surv_hdr - y_step,
             color = LANE_ALIVE, linewidth = 2.2, alpha = 0.85) +
    annotate("text", x = 0.20, y = y_surv_hdr - y_step,
             label = "Alive / censored", hjust = 0, size = sz_tx) +
    annotate("segment", x = 0.06, xend = 0.18,
             y = y_surv_hdr - 2 * y_step, yend = y_surv_hdr - 2 * y_step,
             color = LANE_DEAD, linewidth = 2.2, alpha = 0.45) +
    annotate("text", x = 0.20, y = y_surv_hdr - 2 * y_step,
             label = "Deceased", hjust = 0, size = sz_tx)

  y_ep_hdr <- y_surv_hdr - 3.0 * y_step
  p <- p +
    annotate("text", x = 0.03, y = y_ep_hdr, label = "Endpoint",
             hjust = 0, size = sz_hdr, fontface = "bold") +
    annotate("point", x = 0.06, y = y_ep_hdr - y_step,
             shape = 18, size = 3.8, color = "#1A1A1A") +
    annotate("text", x = 0.11, y = y_ep_hdr - y_step,
             label = "Death", hjust = 0, size = sz_tx) +
    annotate("segment", x = 0.06, xend = 0.12,
             y = y_ep_hdr - 2 * y_step, yend = y_ep_hdr - 2 * y_step,
             arrow = arrow(length = unit(0.15, "cm")),
             color = "black", linewidth = 0.6) +
    annotate("text", x = 0.14, y = y_ep_hdr - 2 * y_step,
             label = "Ongoing (alive)", hjust = 0, size = sz_tx)

  y_brk_hdr <- y_ep_hdr - 3.0 * y_step
  p <- p +
    annotate("text", x = 0.03, y = y_brk_hdr, label = "Axis break",
             hjust = 0, size = sz_hdr, fontface = "bold") +
    annotate("text", x = 0.06, y = y_brk_hdr - y_step,
             label = "//  Truncated at 60 mo", hjust = 0, size = 3.0,
             color = "#666666")

  p + theme_void() +
    coord_fixed(ratio = 1.2) +
    theme(plot.background = element_rect(fill = "transparent", color = NA))
}

build_legend_bottom <- function() {
  sz <- 3.2
  pt <- 4.2
  cw <- 1 / 7  # 7 columns
  y1 <- 0.90; y2 <- 0.72; y3 <- 0.58; y4 <- 0.44
  dx_pt <- 0.015; dx_tx <- 0.04

  col_x <- function(i) (i - 1) * cw  # left edge of column i (1-based)

  ggplot() +
    xlim(0, 1) + ylim(0.35, 1) +
    # Col 1: Distant metastasis
    annotate("text", x = col_x(1), y = y1, label = "Distant metastasis",
             hjust = 0, size = sz, fontface = "bold") +
    annotate("point", x = col_x(1) + dx_pt, y = y2, shape = 22, size = pt, fill = "#00B26E", color = "black") +
    annotate("text", x = col_x(1) + dx_tx, y = y2, label = "Yes", hjust = 0, size = sz) +
    annotate("point", x = col_x(1) + dx_pt, y = y3, shape = 22, size = pt, fill = "#B2E1C4", color = "black") +
    annotate("text", x = col_x(1) + dx_tx, y = y3, label = "No", hjust = 0, size = sz) +
    # Col 2: Max tumor diameter
    annotate("text", x = col_x(2), y = y1, label = "Max tumor diameter",
             hjust = 0, size = sz, fontface = "bold") +
    annotate("point", x = col_x(2) + dx_pt, y = y2, shape = 22, size = pt, fill = "#BC59D3", color = "black") +
    annotate("text", x = col_x(2) + dx_tx, y = y2, label = ">10 cm", hjust = 0, size = sz) +
    annotate("point", x = col_x(2) + dx_pt, y = y3, shape = 22, size = pt, fill = "#E9CBF2", color = "black") +
    annotate("text", x = col_x(2) + dx_tx, y = y3, label = "<=10 cm", hjust = 0, size = sz) +
    # Col 3: Child-Pugh
    annotate("text", x = col_x(3), y = y1, label = "Child-Pugh",
             hjust = 0, size = sz, fontface = "bold") +
    annotate("point", x = col_x(3) + dx_pt, y = y2, shape = 22, size = pt, fill = "#EEDCD6", color = "black") +
    annotate("text", x = col_x(3) + dx_tx, y = y2, label = "A", hjust = 0, size = sz) +
    annotate("point", x = col_x(3) + dx_pt, y = y3, shape = 22, size = pt, fill = "#D38D7B", color = "black") +
    annotate("text", x = col_x(3) + dx_tx, y = y3, label = "B", hjust = 0, size = sz) +
    annotate("point", x = col_x(3) + dx_pt, y = y4, shape = 22, size = pt, fill = "#955442", color = "black") +
    annotate("text", x = col_x(3) + dx_tx, y = y4, label = "C / NA", hjust = 0, size = sz) +
    # Col 4: PVTT
    annotate("text", x = col_x(4), y = y1, label = "PVTT",
             hjust = 0, size = sz, fontface = "bold.italic") +
    annotate("point", x = col_x(4) + dx_pt, y = y2, shape = 22, size = pt, fill = "#DF3E9C", color = "black") +
    annotate("text", x = col_x(4) + dx_tx, y = y2, label = "Vp3/4", hjust = 0, size = sz) +
    annotate("point", x = col_x(4) + dx_pt, y = y3, shape = 22, size = pt, fill = "#FAC1E3", color = "black") +
    annotate("text", x = col_x(4) + dx_tx, y = y3, label = "Absent", hjust = 0, size = sz) +
    annotate("point", x = col_x(4) + dx_pt, y = y4, shape = 22, size = pt, fill = "#F5F5F5", color = "black") +
    annotate("text", x = col_x(4) + dx_tx, y = y4, label = "Vp1/2", hjust = 0, size = sz) +
    # Col 5: AFP
    annotate("text", x = col_x(5), y = y1, label = "AFP",
             hjust = 0, size = sz, fontface = "bold.italic") +
    annotate("point", x = col_x(5) + dx_pt, y = y2, shape = 22, size = pt, fill = "#E64B35", color = "black") +
    annotate("text", x = col_x(5) + dx_tx, y = y2, label = ">400", hjust = 0, size = sz) +
    annotate("point", x = col_x(5) + dx_pt, y = y3, shape = 22, size = pt, fill = "#F7C5BD", color = "black") +
    annotate("text", x = col_x(5) + dx_tx, y = y3, label = "<=400", hjust = 0, size = sz) +
    # Col 6: PIVKA-II
    annotate("text", x = col_x(6), y = y1, label = "PIVKA-II",
             hjust = 0, size = sz, fontface = "bold.italic") +
    annotate("point", x = col_x(6) + dx_pt, y = y2, shape = 22, size = pt, fill = "#3C5488", color = "black") +
    annotate("text", x = col_x(6) + dx_tx, y = y2, label = ">8000", hjust = 0, size = sz) +
    annotate("point", x = col_x(6) + dx_pt, y = y3, shape = 22, size = pt, fill = "#B8C5DB", color = "black") +
    annotate("text", x = col_x(6) + dx_tx, y = y3, label = "<=8000", hjust = 0, size = sz) +
    # Col 7: HAIC courses
    annotate("text", x = col_x(7), y = y1, label = "HAIC courses",
             hjust = 0, size = sz, fontface = "bold") +
    annotate("point", x = col_x(7) + dx_pt, y = y2, shape = 22, size = pt, fill = "#E2EEF5", color = "black") +
    annotate("text", x = col_x(7) + dx_tx, y = y2, label = "1", hjust = 0, size = sz) +
    annotate("point", x = col_x(7) + dx_pt, y = y3, shape = 22, size = pt, fill = "#8FB9D1", color = "black") +
    annotate("text", x = col_x(7) + dx_tx, y = y3, label = "2", hjust = 0, size = sz) +
    annotate("point", x = col_x(7) + dx_pt, y = y4, shape = 22, size = pt, fill = "#2389AF", color = "black") +
    annotate("text", x = col_x(7) + dx_tx, y = y4, label = "3+", hjust = 0, size = sz) +
    theme_void() +
    coord_fixed(ratio = 0.22) +
    theme(plot.background = element_rect(fill = "white", color = NA))
}

# ── Zigzag break helper ──────────────────────────────────────────
# Draws a small "//" zigzag at the end of truncated lanes
make_zigzag_df <- function(indices, x_start, amp = 0.25, n_teeth = 3) {
  if (length(indices) == 0) return(tibble(x = numeric(), y = numeric(), grp = integer()))
  tooth_w <- 1.0 / n_teeth
  purrr::map_dfr(indices, function(idx) {
    xs <- numeric()
    ys <- numeric()
    for (k in seq_len(n_teeth)) {
      xs <- c(xs,
              x_start + (k - 1) * tooth_w,
              x_start + (k - 0.5) * tooth_w,
              x_start + k * tooth_w)
      ys <- c(ys,
              idx - amp,
              idx + amp,
              idx - amp)
    }
    tibble(x = xs, y = ys, grp = idx)
  })
}

# ── Build one group ──────────────────────────────────────────────
build_one_group <- function(group_key, bl, events_raw, max_n, rng_seed) {
  pool <- bl %>% filter(main_group == group_key)
  n_take <- min(as.integer(max_n), nrow(pool))
  if (n_take == 0L) return(invisible(NULL))

  set.seed(rng_seed + which(GROUP_ORDER == group_key))
  sampled <- pool %>%
    slice_sample(n = n_take, replace = FALSE) %>%
    arrange(desc(os_months), patient_id) %>%
    mutate(
      index = row_number(),
      haic_line_cat = case_when(
        haic_episodes <= 1L ~ 1L,
        haic_episodes == 2L ~ 2L,
        TRUE ~ 3L
      ),
      truncated = os_months > X_DISPLAY_MAX,
      os_display = pmin(os_months, X_DISPLAY_MAX),
      lane_color = if_else(death_status == 1L, LANE_DEAD, LANE_ALIVE),
      lane_alpha = if_else(death_status == 1L, 0.45, 0.85)
    )

  # 1:1 mapping of treatment_category → event_type (no merging)
  ev <- events_raw %>%
    semi_join(sampled %>% distinct(patient_id), by = "patient_id") %>%
    left_join(sampled %>% select(patient_id, first_haic_date, os_months, os_display, index),
              by = "patient_id") %>%
    mutate(
      start_date = as.Date(start_date),
      time_m = as.numeric(start_date - first_haic_date, units = "days") / 30.44,
      event_type = case_when(
        treatment_category == "HAIC"              ~ "HAIC",
        treatment_category == "Immunotherapy"     ~ "Immunotherapy",
        treatment_category == "Targeted Therapy"  ~ "Targeted",
        treatment_category == "TACE"              ~ "TACE",
        treatment_category == "HAIC+TACE"         ~ "HAIC+TACE",
        treatment_category == "Resection"         ~ "Resection",
        treatment_category == "Ablation"          ~ "Ablation",
        treatment_category == "Liver Transplant"  ~ "Transplant",
        TRUE                                      ~ NA_character_
      ),
      event_type = factor(event_type, levels = TX_LEVELS)
    ) %>%
    filter(!is.na(time_m), !is.na(event_type),
           time_m >= 0, time_m <= X_DISPLAY_MAX)

  # Y-jitter for overlapping events on the same date
  # Group by patient + date, assign vertical offset so stacked types don't hide each other
  ev <- ev %>%
    group_by(patient_id, start_date) %>%
    mutate(
      n_same = n(),
      rank_in_day = row_number(),
      y_jitter = if_else(n_same > 1L,
                         index + (rank_in_day - (n_same + 1) / 2) * 0.18,
                         as.double(index))
    ) %>%
    ungroup()

  # Endpoint markers (clamped to display max)
  death_pts <- sampled %>%
    filter(death_status == 1L) %>%
    transmute(index, time_m = pmin(os_months, X_DISPLAY_MAX))

  alive_pts <- sampled %>%
    filter(death_status != 1L | is.na(death_status)) %>%
    transmute(index, time_m = pmin(os_months, X_DISPLAY_MAX))

  # Zigzag data for truncated lanes
  trunc_idx <- sampled$index[sampled$truncated]
  zz <- make_zigzag_df(trunc_idx, x_start = X_DISPLAY_MAX - 1.0, amp = 0.22, n_teeth = 3)

  # Left-side annotations (7 columns)
  anno <- sampled %>%
    transmute(index, distant_metastasis, tumor_size_category,
              child_pugh_grade, pvtt_classification,
              afp_high, pivka_high, haic_line_cat)

  idx_max <- max(sampled$index)
  pretty_lab <- GROUP_PRETTY[[group_key]]
  n_total <- nrow(pool)

  # ── Main plot ──────────────────────────────────────────────────
  # Draw alive and dead lanes separately to control alpha per group
  alive_lanes <- sampled %>% filter(death_status != 1L | is.na(death_status))
  dead_lanes  <- sampled %>% filter(death_status == 1L)

  p <- ggplot()

  if (nrow(alive_lanes) > 0) {
    p <- p + geom_segment(
      data = alive_lanes,
      aes(x = 0, xend = os_display, y = index, yend = index),
      color = LANE_ALIVE, linewidth = 2.5, alpha = 0.85, lineend = "round"
    )
  }
  if (nrow(dead_lanes) > 0) {
    p <- p + geom_segment(
      data = dead_lanes,
      aes(x = 0, xend = os_display, y = index, yend = index),
      color = LANE_DEAD, linewidth = 2.5, alpha = 0.45, lineend = "round"
    )
  }

  # Split events: common types first, then rare types (TACE, HAIC+TACE,
  # Resection, Ablation, Transplant) on top with slightly larger size
  rare_types <- c("TACE", "HAIC+TACE", "Resection", "Ablation", "Transplant")
  ev_common <- ev %>% filter(!event_type %in% rare_types)
  ev_rare   <- ev %>% filter(event_type %in% rare_types)

  p <- p +
    geom_point(
      data = ev_common,
      aes(x = time_m, y = y_jitter, color = event_type, shape = event_type),
      size = 2.0, alpha = 0.90, stroke = 0.3
    )
  if (nrow(ev_rare) > 0) {
    p <- p +
      geom_point(
        data = ev_rare,
        aes(x = time_m, y = y_jitter, color = event_type, shape = event_type),
        size = 2.8, alpha = 0.95, stroke = 0.4
      )
  }
  p <- p +
    scale_color_manual(values = TX_COLORS, drop = FALSE, guide = "none") +
    scale_shape_manual(values = TX_SHAPES, drop = FALSE, guide = "none")

  # Death endpoint
  if (nrow(death_pts) > 0) {
    p <- p + geom_point(
      data = death_pts, aes(x = time_m, y = index),
      shape = 18, size = 3, color = "#1A1A1A"
    )
  }
  # Alive endpoint (arrow)
  if (nrow(alive_pts) > 0) {
    p <- p + geom_segment(
      data = alive_pts,
      aes(x = time_m, xend = time_m + 1.2, y = index, yend = index),
      arrow = arrow(length = unit(0.15, "cm")),
      color = "black", linewidth = 0.4
    )
  }

  # Zigzag break marks for truncated lanes
  if (nrow(zz) > 0) {
    p <- p + geom_line(
      data = zz, aes(x = x, y = y, group = grp),
      color = "white", linewidth = 1.8
    ) +
    geom_line(
      data = zz, aes(x = x, y = y, group = grp),
      color = "#333333", linewidth = 0.6
    )
  }

  # ── Left-side 7-column annotation squares (spacing = 1.4) ─────
  xs <- -seq(1.4, by = 1.4, length.out = 7)  # x positions for 7 columns
  aps <- 3.2  # annotation point size

  p <- p +
    # Col 1: Distant metastasis
    geom_point(data = subset(anno, distant_metastasis == "Yes"),
               aes(x = xs[1], y = index), inherit.aes = FALSE,
               fill = "#00B26E", color = "black", shape = 22, size = aps) +
    geom_point(data = subset(anno, is.na(distant_metastasis) | distant_metastasis != "Yes"),
               aes(x = xs[1], y = index), inherit.aes = FALSE,
               fill = "#B2E1C4", color = "black", shape = 22, size = aps) +
    # Col 2: Tumor diameter
    geom_point(data = subset(anno, tumor_size_category == ">10cm"),
               aes(x = xs[2], y = index), inherit.aes = FALSE,
               fill = "#BC59D3", color = "black", shape = 22, size = aps) +
    geom_point(data = subset(anno, is.na(tumor_size_category) | tumor_size_category != ">10cm"),
               aes(x = xs[2], y = index), inherit.aes = FALSE,
               fill = "#E9CBF2", color = "black", shape = 22, size = aps) +
    # Col 3: Child-Pugh
    geom_point(data = subset(anno, child_pugh_grade == "A"),
               aes(x = xs[3], y = index), inherit.aes = FALSE,
               fill = "#EEDCD6", color = "black", shape = 22, size = aps) +
    geom_point(data = subset(anno, child_pugh_grade == "B"),
               aes(x = xs[3], y = index), inherit.aes = FALSE,
               fill = "#D38D7B", color = "black", shape = 22, size = aps) +
    geom_point(data = subset(anno, is.na(child_pugh_grade) | !child_pugh_grade %in% c("A", "B")),
               aes(x = xs[3], y = index), inherit.aes = FALSE,
               fill = "#955442", color = "black", shape = 22, size = aps) +
    # Col 4: PVTT
    geom_point(data = subset(anno, pvtt_classification == "Vp3/4"),
               aes(x = xs[4], y = index), inherit.aes = FALSE,
               fill = "#DF3E9C", color = "black", shape = 22, size = aps) +
    geom_point(data = subset(anno, is.na(pvtt_classification) | pvtt_classification == "Absent"),
               aes(x = xs[4], y = index), inherit.aes = FALSE,
               fill = "#FAC1E3", color = "black", shape = 22, size = aps) +
    geom_point(data = subset(anno, pvtt_classification == "Vp1/2"),
               aes(x = xs[4], y = index), inherit.aes = FALSE,
               fill = "#F5F5F5", color = "black", shape = 22, size = aps) +
    # Col 5: AFP (>400 / <=400)
    geom_point(data = subset(anno, afp_high == "Yes"),
               aes(x = xs[5], y = index), inherit.aes = FALSE,
               fill = "#E64B35", color = "black", shape = 22, size = aps) +
    geom_point(data = subset(anno, is.na(afp_high) | afp_high != "Yes"),
               aes(x = xs[5], y = index), inherit.aes = FALSE,
               fill = "#F7C5BD", color = "black", shape = 22, size = aps) +
    # Col 6: PIVKA-II (>8000 / <=8000)
    geom_point(data = subset(anno, pivka_high == "Yes"),
               aes(x = xs[6], y = index), inherit.aes = FALSE,
               fill = "#3C5488", color = "black", shape = 22, size = aps) +
    geom_point(data = subset(anno, is.na(pivka_high) | pivka_high != "Yes"),
               aes(x = xs[6], y = index), inherit.aes = FALSE,
               fill = "#B8C5DB", color = "black", shape = 22, size = aps) +
    # Col 7: HAIC courses
    geom_point(data = subset(anno, haic_line_cat == 1),
               aes(x = xs[7], y = index), inherit.aes = FALSE,
               fill = "#E2EEF5", color = "black", shape = 22, size = aps) +
    geom_point(data = subset(anno, haic_line_cat == 2),
               aes(x = xs[7], y = index), inherit.aes = FALSE,
               fill = "#8FB9D1", color = "black", shape = 22, size = aps) +
    geom_point(data = subset(anno, haic_line_cat >= 3),
               aes(x = xs[7], y = index), inherit.aes = FALSE,
               fill = "#2389AF", color = "black", shape = 22, size = aps)

  # Manual axes
  p <- p +
    geom_segment(data = tibble(x = 0, xend = 0, y = 0, yend = idx_max + 0.4),
                 aes(x = x, xend = xend, y = y, yend = yend),
                 color = "black", linewidth = 0.5) +
    geom_segment(data = tibble(x = 0, xend = X_DISPLAY_MAX, y = 0, yend = 0),
                 aes(x = x, xend = xend, y = y, yend = yend),
                 color = "black", linewidth = 0.5) +
    scale_x_continuous(breaks = seq(0, X_DISPLAY_MAX, 12)) +
    coord_cartesian(xlim = c(-10.5, X_DISPLAY_MAX + 1.5), clip = "off") +
    scale_y_continuous(expand = expansion(mult = c(0, 0), add = c(0, 0.5))) +
    xlab("Time from first HAIC (months)") +
    ylab("") +
    theme_classic(base_size = 10, base_family = "Helvetica") +
    theme(
      axis.text.y         = element_blank(),
      axis.line            = element_blank(),
      axis.ticks.x         = element_line(linewidth = 0.5, color = "black"),
      axis.ticks.length.x  = unit(0.15, "cm"),
      axis.ticks.y         = element_blank(),
      legend.position      = "none",
      axis.title.x         = element_text(size = 11, color = "black"),
      axis.text.x          = element_text(size = 9, color = "black"),
      plot.title           = element_text(face = "bold", size = 12, hjust = 0.5),
      plot.subtitle        = element_text(size = 8, color = "gray40", hjust = 0.5),
      plot.margin          = margin(8, 8, 4, 8)
    ) +
    labs(
      title = sprintf("Swimmer plot: %s", pretty_lab),
      subtitle = sprintf("Showing %d of %d patients (seed %d)", n_take, n_total, rng_seed)
    )

  top_part <- p +
    inset_element(build_legend_right(), left = 0.58, bottom = 0.18, right = 0.99, top = 0.99)

  top_part / build_legend_bottom() + plot_layout(heights = c(3.5, 1))
}

# ── Export all groups ────────────────────────────────────────────
saved <- character(0)
for (gk in GROUP_ORDER) {
  comb <- build_one_group(gk, bl, events_raw, MAX_PER_GROUP, RNG_SEED)
  if (is.null(comb)) next
  fn <- safe_fname(gk)
  out_pdf <- file.path(FIG_DIR, paste0("swimmer_", fn, ".pdf"))
  out_png <- file.path(FIG_DIR, paste0("swimmer_", fn, ".png"))
  ggsave(out_pdf, comb, width = 9, height = 8)
  ggsave(out_png, comb, width = 9, height = 8, dpi = 300)
  saved <- c(saved, out_pdf)
}

message("Saved ", length(saved), " figures under: ", FIG_DIR)
if (length(saved)) writeLines(basename(saved))
