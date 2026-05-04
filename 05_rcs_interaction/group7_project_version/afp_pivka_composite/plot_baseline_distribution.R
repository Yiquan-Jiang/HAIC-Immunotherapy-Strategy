#!/usr/bin/env Rscript
# AFP & PIVKA baseline distribution analysis
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(gridExtra)
  library(grid)
  library(scales)
})

SCRIPT_DIR <- "."
base_dir <- ".."
df <- read.csv(file.path(base_dir, "HAIC_NO_TACE_4_TIDY_baseline.csv"), stringsAsFactors = FALSE)
cat("Total n:", nrow(df), "\n\n")

# === 频率表 ===
cat("=== AFP baseline ===\n")
cat("Median:", median(df$afp, na.rm=TRUE), "| Mean:", round(mean(df$afp, na.rm=TRUE),1), "\n")
tbl_afp <- table(cut(df$afp, breaks=c(0,20,200,400,1000,Inf),
    labels=c("<20","20-200","200-400","400-1000",">1000"), right=FALSE))
print(tbl_afp)
cat("AFP<20:", sum(df$afp < 20, na.rm=TRUE), "(", round(100*mean(df$afp < 20, na.rm=TRUE),1), "%)\n\n")

cat("=== PIVKA baseline ===\n")
cat("Median:", median(df$pivka, na.rm=TRUE), "| Mean:", round(mean(df$pivka, na.rm=TRUE),1), "\n")
tbl_pivka <- table(cut(df$pivka, breaks=c(0,20,40,100,400,Inf),
    labels=c("<20","20-40","40-100","100-400",">400"), right=FALSE))
print(tbl_pivka)
cat("PIVKA<20:", sum(df$pivka < 20, na.rm=TRUE), "(", round(100*mean(df$pivka < 20, na.rm=TRUE),1), "%)\n\n")

# === 交叉表 ===
df$afp_neg <- df$afp < 20
df$pivka_neg <- df$pivka < 20
cross <- table(AFP=df$afp_neg, PIVKA=df$pivka_neg)
rownames(cross) <- c("AFP>=20 (pos)", "AFP<20 (neg)")
colnames(cross) <- c("PIVKA>=20 (pos)", "PIVKA<20 (neg)")
cat("=== AFP x PIVKA cross table ===\n")
print(cross)
cat("\nRow %:\n"); print(round(prop.table(cross, 1)*100, 1))
cat("\nCol %:\n"); print(round(prop.table(cross, 2)*100, 1))
cat("\nTotal %:\n"); print(round(prop.table(cross)*100, 1))

# === 绘图 ===
outdir <- file.path(SCRIPT_DIR, "output")
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

df$group4 <- with(df, case_when(
  afp_neg & pivka_neg ~ "Both neg\n(AFP<20, PIVKA<20)",
  afp_neg & !pivka_neg ~ "AFP neg only\n(AFP<20, PIVKA>=20)",
  !afp_neg & pivka_neg ~ "PIVKA neg only\n(AFP>=20, PIVKA<20)",
  !afp_neg & !pivka_neg ~ "Both pos\n(AFP>=20, PIVKA>=20)"
))
grp_ct <- df %>% count(group4) %>% mutate(pct=n/sum(n)*100)
grp_ct$label_str <- sprintf("%d\n(%.1f%%)", grp_ct$n, grp_ct$pct)

p1 <- ggplot(df, aes(x=afp)) +
  geom_histogram(bins=60, fill="#4393C3", alpha=0.7, color="white") +
  geom_vline(xintercept=20, linetype="dashed", color="#DC0000", linewidth=1) +
  annotate("text", x=30, y=Inf, label="AFP=20", color="#DC0000", vjust=1.8, size=3.5, fontface="bold") +
  scale_x_log10(labels=comma, breaks=c(1,5,10,20,50,100,500,1000,5000,10000)) +
  labs(title="AFP Baseline Distribution (log scale)", x="AFP (ng/mL)", y="Count") +
  theme_bw(base_size=13)

p2 <- ggplot(df, aes(x=pivka)) +
  geom_histogram(bins=60, fill="#D6604D", alpha=0.7, color="white") +
  geom_vline(xintercept=20, linetype="dashed", color="#DC0000", linewidth=1) +
  annotate("text", x=30, y=Inf, label="PIVKA=20", color="#DC0000", vjust=1.8, size=3.5, fontface="bold") +
  scale_x_log10(labels=comma, breaks=c(1,5,10,20,50,100,500,1000,5000)) +
  labs(title="PIVKA Baseline Distribution (log scale)", x="PIVKA (mAU/mL)", y="Count") +
  theme_bw(base_size=13)

dbl_neg <- sum(df$afp_neg & df$pivka_neg, na.rm=TRUE)
p3 <- ggplot(df, aes(x=afp, y=pivka)) +
  geom_point(alpha=0.12, size=0.7, color="#333333") +
  geom_hline(yintercept=20, linetype="dashed", color="#D6604D", linewidth=0.8) +
  geom_vline(xintercept=20, linetype="dashed", color="#4393C3", linewidth=0.8) +
  annotate("rect", xmin=1, xmax=20, ymin=1, ymax=20, fill="#2E7D32", alpha=0.15) +
  scale_x_log10(labels=comma, limits=c(1,20000), breaks=c(1,10,20,100,1000,10000)) +
  scale_y_log10(labels=comma, limits=c(1,10000), breaks=c(1,10,20,100,1000,5000)) +
  annotate("text", x=3, y=3, label=sprintf("Double Neg\nn=%d", dbl_neg),
           color="#2E7D32", size=3.5, fontface="bold") +
  labs(title="AFP vs PIVKA at Baseline", x="AFP (ng/mL)", y="PIVKA (mAU/mL)") +
  theme_bw(base_size=13)

bar_colors <- c("#2E7D32","#66BB6A","#FFA726","#EF5350")
p4 <- ggplot(grp_ct, aes(x=group4, y=n, fill=group4)) +
  geom_bar(stat="identity", color="white", linewidth=0.5) +
  geom_text(aes(label=label_str), vjust=-0.3, size=3.5, fontface="bold") +
  scale_fill_manual(values=bar_colors) +
  labs(title="Four-Quadrant Distribution", x="", y="Patient Count") +
  theme_bw(base_size=13) +
  theme(legend.position="none", axis.text.x=element_text(size=9))

comb <- grid.arrange(p1, p2, p3, p4, ncol=2,
  top=textGrob(sprintf("AFP & PIVKA Baseline Distribution (n=%d)", nrow(df)),
               gp=gpar(fontsize=15, fontface="bold")))
ggsave(file.path(outdir, "afp_pivka_baseline_distribution.png"), comb, width=14, height=12, dpi=300)
ggsave(file.path(outdir, "afp_pivka_baseline_distribution.pdf"), comb, width=14, height=12, device="pdf")
cat("\nFigure saved.\n")
