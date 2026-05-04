# ---------------------------------------------------------------------------
# HAIC then I — interval from immunotherapy to first targeted therapy
#
# Reconstructed from Claude Code session 7e3ed3cb-... (foamy-munching-riddle),
# final iteration 2026-04-10 07:07 UTC. Original code was an inline heredoc
# piped to `Rscript --vanilla` and was never saved as a file.
#
# Output:
#   ../figures/haic_then_i_immuno_to_target_interval.pdf
#   ../figures/haic_then_i_immuno_to_target_interval.png
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
})

BASE_DIR <- "/Users/yqj/Nutstore Files/我的坚果云/Liver_tumor_big_data/FIRST_LINE_HAIC_2025-12-30/HAIC_NO_TACE_4_TIDY/update_group_7"
DATA_DIR <- file.path(BASE_DIR, "data")
FIG_DIR  <- file.path(BASE_DIR, "figures")

bl <- read_csv(file.path(DATA_DIR, "analysis_ready.csv"),
               col_types = cols(.default = col_character()),
               show_col_types = FALSE) %>%
  type_convert(col_types = cols())

ev <- read_csv(file.path(DATA_DIR, "00_swimmer_plot_events.csv"),
               col_types = cols(.default = col_character()),
               show_col_types = FALSE) %>%
  type_convert(col_types = cols())

haic_then_i <- bl %>% filter(main_group == "HAIC_then_I")
n_total <- nrow(haic_then_i)
ev_haic_i <- ev %>% filter(patient_id %in% haic_then_i$patient_id)

first_immuno <- ev_haic_i %>%
  filter(treatment_category == "Immunotherapy") %>%
  mutate(start_date = as.Date(start_date)) %>%
  group_by(patient_id) %>%
  summarize(first_immuno_date = min(start_date), .groups = "drop")

first_target <- ev_haic_i %>%
  filter(treatment_category == "Targeted Therapy") %>%
  mutate(start_date = as.Date(start_date)) %>%
  group_by(patient_id) %>%
  arrange(start_date) %>%
  slice(1) %>%
  ungroup() %>%
  select(patient_id, first_target_date = start_date, drug_ev = drug_name)

merged <- first_immuno %>%
  inner_join(first_target, by = "patient_id") %>%
  mutate(interval_days   = as.numeric(first_target_date - first_immuno_date),
         interval_months = interval_days / 30.44) %>%
  left_join(haic_then_i %>% select(patient_id, drug_bl = first_target_drug),
            by = "patient_id") %>%
  mutate(drug = coalesce(drug_bl, drug_ev))

n_target    <- nrow(merged)
n_no_target <- n_total - n_target

drug_counts <- merged %>%
  count(drug, sort = TRUE) %>%
  mutate(pct = n / sum(n) * 100)

drug_cols <- c(
  "Lenvatinib"  = "#E64B35",
  "Bevacizumab" = "#4DBBD5",
  "Apatinib"    = "#00A087",
  "Regorafenib" = "#3C5488",
  "Sorafenib"   = "#F39B7F",
  "Donafenib"   = "#8491B4"
)
used_drugs <- drug_counts$drug
cols_use   <- drug_cols[used_drugs]
cols_use[is.na(cols_use)] <- "#999999"
drug_counts$drug <- factor(drug_counts$drug, levels = used_drugs)

# ── Histogram (x-axis capped at 420 days; outliers squished) ──
p_hist <- ggplot(merged, aes(x = interval_days)) +
  geom_histogram(binwidth = 30, fill = "#4DBBD5",
                 color = "white", alpha = 0.85) +
  geom_vline(xintercept = median(merged$interval_days),
             linetype = "dashed", color = "#E64B35", linewidth = 0.8) +
  annotate("text",
           x = median(merged$interval_days) + 20, y = Inf,
           vjust = 2, hjust = 0,
           label = sprintf("Median = %d days (%.1f mo)",
                           median(merged$interval_days),
                           median(merged$interval_months)),
           color = "#E64B35", size = 3.5, fontface = "bold") +
  annotate("label",
           x = Inf, y = Inf, vjust = 1.3, hjust = 1.02,
           label = "Rule: Targeted therapy was permitted\n>= 1 month after immunotherapy initiation",
           size = 3.0, fontface = "italic", color = "#3C5488",
           fill = "#F0F4FA",
           label.size = 0.4, label.r = unit(0.2, "lines")) +
  scale_x_continuous(limits = c(0, 420),
                     breaks = seq(0, 420, 60),
                     oob = scales::squish) +
  labs(x = "Days from immunotherapy to targeted therapy",
       y = "Number of patients") +
  theme_classic(base_size = 11) +
  theme(plot.margin = margin(5, 8, 5, 5))

# ── Pie 1: targeted vs no targeted ──
pie1_df <- tibble(
  group = factor(c("Added targeted", "No targeted"),
                 levels = c("Added targeted", "No targeted")),
  n     = c(n_target, n_no_target),
  lab   = c(sprintf("%d (%.1f%%)", n_target,    n_target    / n_total * 100),
            sprintf("%d (%.1f%%)", n_no_target, n_no_target / n_total * 100))
)

p_pie1 <- ggplot(pie1_df, aes(x = "", y = n, fill = group)) +
  geom_col(width = 1, color = "white", linewidth = 1) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = c("Added targeted" = "#3C5488",
                               "No targeted"    = "#B8C5DB"),
                    labels = paste(pie1_df$group, pie1_df$lab),
                    name = NULL) +
  labs(title = sprintf("HAIC then I (n = %d)", n_total)) +
  theme_void(base_size = 11) +
  theme(legend.position = "bottom",
        legend.text     = element_text(size = 9),
        legend.key.size = unit(0.45, "cm"),
        plot.title      = element_text(face = "bold", size = 11, hjust = 0.5))

# ── Pie 2: drug types ──
p_pie2 <- ggplot(drug_counts, aes(x = "", y = n, fill = drug)) +
  geom_col(width = 1, color = "white", linewidth = 0.8) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = cols_use, name = NULL,
                    labels = sprintf("%s  %d (%.0f%%)",
                                     drug_counts$drug,
                                     drug_counts$n,
                                     drug_counts$pct)) +
  labs(title = sprintf("Targeted agents (n = %d)", n_target)) +
  theme_void(base_size = 11) +
  theme(legend.position = "bottom",
        legend.text     = element_text(size = 8.5),
        legend.key.size = unit(0.4, "cm"),
        plot.title      = element_text(face = "bold", size = 11, hjust = 0.5))

# ── Combine ──
p_right  <- p_pie1 / p_pie2
combined <- (p_hist | p_right) +
  plot_layout(widths = c(1.3, 1)) +
  plot_annotation(
    title    = "HAIC then I: Targeted Therapy After Immunotherapy",
    subtitle = sprintf(
      "%d of %d patients (%.1f%%) added targeted therapy; median interval %.1f months (IQR %.1f-%.1f)",
      n_target, n_total, n_target / n_total * 100,
      median(merged$interval_months),
      quantile(merged$interval_months, 0.25),
      quantile(merged$interval_months, 0.75)
    ),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 10, color = "gray30", hjust = 0.5)
    )
  )

ggsave(file.path(FIG_DIR, "haic_then_i_immuno_to_target_interval.pdf"),
       combined, width = 11, height = 6.5)
ggsave(file.path(FIG_DIR, "haic_then_i_immuno_to_target_interval.png"),
       combined, width = 11, height = 6.5, dpi = 300)
cat("Done\n")
