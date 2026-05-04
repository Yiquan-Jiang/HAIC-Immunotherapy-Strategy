#!/usr/bin/env Rscript
# AFP 阴性患者（AFP<20）的 AFP & PIVKA 随 HAIC 疗程变化折线图
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(gridExtra)
  library(scales)
  library(grid)
})

base_dir <- ".."

# ── 读取数据 ──
df_long <- read.csv(file.path(base_dir, "HAIC_NO_TACE_4_TIDY_longitudinal.csv"), stringsAsFactors = FALSE)
df_base <- read.csv(file.path(base_dir, "HAIC_NO_TACE_4_TIDY_baseline.csv"), stringsAsFactors = FALSE)

afp_neg_ids <- df_base$patient_id[df_base$afp < 20]
cat("AFP negative patients (baseline AFP<20):", length(afp_neg_ids), "\n")

df_sub <- df_long[df_long$patient_id %in% afp_neg_ids, ]
df_sub <- df_sub[!grepl("^post", df_sub$timepoint_type), ]

# 提取疗程序号
df_sub$haic_num <- NA_real_
df_sub$haic_num[df_sub$timepoint_type == "baseline"] <- 1
m <- grepl("^pre_haic_(\\d+)$", df_sub$timepoint_type)
df_sub$haic_num[m] <- as.numeric(sub("^pre_haic_(\\d+)$", "\\1", df_sub$timepoint_type[m])) + 1

# 限制到 haic_num <= 10
df_plot <- df_sub[df_sub$haic_num <= 10, ]

# 用 haic_num 作为 x 轴（数值型，更可靠）
# 同时保留可读标签
tp_labels <- c("1" = "Baseline", "3" = "Pre-HAIC 2", "4" = "Pre-HAIC 3",
               "5" = "Pre-HAIC 4", "6" = "Pre-HAIC 5", "7" = "Pre-HAIC 6",
               "8" = "Pre-HAIC 7", "9" = "Pre-HAIC 8", "10" = "Pre-HAIC 9",
               "11" = "Pre-HAIC 10")

# 基线值（只取有 baseline 记录的患者，去重）
baseline_vals <- df_plot[df_plot$timepoint_type == "baseline" & !is.na(df_plot$afp) & !is.na(df_plot$pivka),
                         c("patient_id", "afp", "pivka")]
baseline_vals <- baseline_vals[!duplicated(baseline_vals$patient_id), ]
colnames(baseline_vals) <- c("patient_id", "afp_bl", "pivka_bl")
cat("Baseline patients with valid afp+pivka:", nrow(baseline_vals), "\n")
df_plot <- merge(df_plot, baseline_vals, by = "patient_id", all.x = TRUE)
cat("After merge:", nrow(df_plot), "\n")
df_plot$pid_num <- as.numeric(factor(df_plot$patient_id))
df_plot$afp_ratio <- ifelse(!is.na(df_plot$afp) & !is.na(df_plot$afp_bl), df_plot$afp / df_plot$afp_bl, NA)
df_plot$pivka_ratio <- ifelse(!is.na(df_plot$pivka) & !is.na(df_plot$pivka_bl), df_plot$pivka / df_plot$pivka_bl, NA)

outdir <- file.path(base_dir, "afp_pivka_composite", "output")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# 统计
tp_counts <- df_plot %>%
  filter(!is.na(haic_num)) %>%
  group_by(haic_num) %>%
  summarise(n = n_distinct(patient_id), .groups = "drop")
cat("\nPatients by HAIC number:\n")
print(tp_counts)

# ═══════════════════════════════════════════════════════════════
# 图1: AFP 折线图
# ═══════════════════════════════════════════════════════════════
afp_sum <- df_plot %>%
  filter(!is.na(afp) & !is.na(haic_num)) %>%
  group_by(haic_num) %>%
  summarise(n = n(), med = median(afp, na.rm=TRUE),
            q25 = quantile(afp, 0.25, na.rm=TRUE),
            q75 = quantile(afp, 0.75, na.rm=TRUE), .groups="drop")

p1 <- ggplot(df_plot, aes(x = haic_num, y = afp, group = pid_num)) +
  geom_line(alpha = 0.04, color = "#4393C3", linewidth = 0.3) +
  geom_line(data = afp_sum, aes(x = haic_num, y = med, group = 1), inherit.aes = FALSE,
            color = "#DC0000", linewidth = 1.8) +
  geom_point(data = afp_sum, aes(x = haic_num, y = med), inherit.aes = FALSE,
             color = "#DC0000", size = 3) +
  geom_ribbon(data = afp_sum, aes(x = haic_num, ymin = q25, ymax = q75, group = 1),
              inherit.aes = FALSE, fill = "#DC0000", alpha = 0.12) +
  geom_hline(yintercept = 20, linetype = "dashed", color = "#333333", linewidth = 0.6) +
  scale_x_continuous(breaks = unique(afp_sum$haic_num), labels = tp_labels[as.character(unique(afp_sum$haic_num))]) +
  scale_y_log10(labels = comma, breaks = c(1, 5, 10, 20, 50, 100, 500, 1000, 5000)) +
  annotate("text", x = 11.5, y = 20, label = "AFP=20", color = "#333333", size = 3) +
  labs(title = "AFP trajectory in AFP-negative patients (baseline AFP<20)",
       subtitle = sprintf("Individual lines + median (red) + IQR ribbon | n_patients at BL: %d",
                          length(afp_neg_ids)),
       x = "", y = "AFP (ng/mL, log scale)") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        panel.grid.minor = element_blank())

# ═══════════════════════════════════════════════════════════════
# 图2: PIVKA 折线图
# ═══════════════════════════════════════════════════════════════
pivka_sum <- df_plot %>%
  filter(!is.na(pivka) & !is.na(haic_num)) %>%
  group_by(haic_num) %>%
  summarise(n = n(), med = median(pivka, na.rm=TRUE),
            q25 = quantile(pivka, 0.25, na.rm=TRUE),
            q75 = quantile(pivka, 0.75, na.rm=TRUE), .groups="drop")

p2 <- ggplot(df_plot, aes(x = haic_num, y = pivka, group = pid_num)) +
  geom_line(alpha = 0.04, color = "#D6604D", linewidth = 0.3) +
  geom_line(data = pivka_sum, aes(x = haic_num, y = med, group = 1), inherit.aes = FALSE,
            color = "#2166AC", linewidth = 1.8) +
  geom_point(data = pivka_sum, aes(x = haic_num, y = med), inherit.aes = FALSE,
             color = "#2166AC", size = 3) +
  geom_ribbon(data = pivka_sum, aes(x = haic_num, ymin = q25, ymax = q75, group = 1),
              inherit.aes = FALSE, fill = "#2166AC", alpha = 0.12) +
  geom_hline(yintercept = 20, linetype = "dashed", color = "#333333", linewidth = 0.6) +
  scale_x_continuous(breaks = unique(pivka_sum$haic_num), labels = tp_labels[as.character(unique(pivka_sum$haic_num))]) +
  scale_y_log10(labels = comma, breaks = c(1, 5, 10, 20, 50, 100, 500, 1000, 5000, 10000)) +
  annotate("text", x = 11.5, y = 20, label = "PIVKA=20", color = "#333333", size = 3) +
  labs(title = "PIVKA trajectory in AFP-negative patients (baseline AFP<20)",
       subtitle = sprintf("Individual lines + median (blue) + IQR ribbon | n_patients at BL: %d",
                          length(afp_neg_ids)),
       x = "", y = "PIVKA (mAU/mL, log scale)") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        panel.grid.minor = element_blank())

# ═══════════════════════════════════════════════════════════════
# 图3: 相对于 baseline 的变化倍数
# ═══════════════════════════════════════════════════════════════
afp_r_sum <- df_plot %>%
  filter(!is.na(afp_ratio) & !is.na(haic_num)) %>%
  group_by(haic_num) %>%
  summarise(n = n(), med = median(afp_ratio, na.rm=TRUE),
            q25 = quantile(afp_ratio, 0.25, na.rm=TRUE),
            q75 = quantile(afp_ratio, 0.75, na.rm=TRUE), .groups="drop")

pivka_r_sum <- df_plot %>%
  filter(!is.na(pivka_ratio) & !is.na(haic_num)) %>%
  group_by(haic_num) %>%
  summarise(n = n(), med = median(pivka_ratio, na.rm=TRUE),
            q25 = quantile(pivka_ratio, 0.25, na.rm=TRUE),
            q75 = quantile(pivka_ratio, 0.75, na.rm=TRUE), .groups="drop")

p3 <- ggplot() +
  geom_line(data = afp_r_sum, aes(x = haic_num, y = med, color = "AFP"), linewidth = 1.8) +
  geom_point(data = afp_r_sum, aes(x = haic_num, y = med, color = "AFP"), size = 3) +
  geom_ribbon(data = afp_r_sum, aes(x = haic_num, ymin = q25, ymax = q75, fill = "AFP"), alpha = 0.12) +
  geom_line(data = pivka_r_sum, aes(x = haic_num, y = med, color = "PIVKA"), linewidth = 1.8) +
  geom_point(data = pivka_r_sum, aes(x = haic_num, y = med, color = "PIVKA"), size = 3) +
  geom_ribbon(data = pivka_r_sum, aes(x = haic_num, ymin = q25, ymax = q75, fill = "PIVKA"), alpha = 0.12) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "#333333", linewidth = 0.8) +
  scale_x_continuous(breaks = unique(afp_r_sum$haic_num), labels = tp_labels[as.character(unique(afp_r_sum$haic_num))]) +
  scale_y_log10(labels = comma, breaks = c(0.1, 0.25, 0.5, 1, 2, 5, 10, 20, 50)) +
  scale_color_manual(name = "Marker", values = c("AFP" = "#4393C3", "PIVKA" = "#D6604D")) +
  scale_fill_manual(name = "Marker", values = c("AFP" = "#4393C3", "PIVKA" = "#D6604D")) +
  labs(title = "Relative change from baseline (ratio)",
       subtitle = "AFP-negative patients | baseline = 1.0",
       x = "", y = "Ratio to baseline (log scale)") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        panel.grid.minor = element_blank(), legend.position = "top")

# ═══════════════════════════════════════════════════════════════
# 图4: 各时间点阴性比例
# ═══════════════════════════════════════════════════════════════
neg_frac <- df_plot %>%
  filter(!is.na(haic_num) & !is.na(afp) & !is.na(pivka)) %>%
  group_by(haic_num) %>%
  summarise(
    n_total = n(),
    n_afp_neg = sum(afp < 20, na.rm = TRUE),
    n_pivka_neg = sum(pivka < 20, na.rm = TRUE),
    pct_afp = if(n_total > 0) 100 * n_afp_neg / n_total else NA_real_,
    pct_pivka = if(n_total > 0) 100 * n_pivka_neg / n_total else NA_real_,
    .groups = "drop"
  )
nf_long <- rbind(
  data.frame(haic_num = neg_frac$haic_num, Marker = "AFP<20",
             Pct = neg_frac$pct_afp, N = neg_frac$n_total),
  data.frame(haic_num = neg_frac$haic_num, Marker = "PIVKA<20",
             Pct = neg_frac$pct_pivka, N = neg_frac$n_total)
)

p4 <- ggplot(nf_long, aes(x = haic_num, y = Pct, color = Marker, group = Marker)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) +
  geom_text(aes(label = sprintf("%.0f%%", Pct)), vjust = -1.5, size = 2.8,
            show.legend = FALSE) +
  scale_x_continuous(breaks = unique(neg_frac$haic_num),
                     labels = tp_labels[as.character(unique(neg_frac$haic_num))]) +
  scale_color_manual(values = c("AFP<20" = "#4393C3", "PIVKA<20" = "#D6604D")) +
  scale_y_continuous(limits = c(0, 110), breaks = seq(0, 100, 20)) +
  labs(title = "Proportion remaining negative over HAIC courses",
       subtitle = sprintf("Among AFP-negative patients at baseline (n=%d)",
                          length(afp_neg_ids)),
       x = "", y = "% patients negative") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        legend.position = "top")

# ── 保存 ──
comb <- grid.arrange(p1, p2, p3, p4, ncol = 2,
  top = textGrob(sprintf("AFP-Negative Cohort (baseline AFP<20): AFP & PIVKA Dynamics (n=%d at baseline)",
                         length(afp_neg_ids)),
                 gp = gpar(fontsize = 15, fontface = "bold")))

ggsave(file.path(outdir, "afp_negative_dynamic_trajectory.png"),
       comb, width = 16, height = 14, dpi = 300)
ggsave(file.path(outdir, "afp_negative_dynamic_trajectory.pdf"),
       comb, width = 16, height = 14, device = "pdf")
cat("\nFigure saved.\n")
