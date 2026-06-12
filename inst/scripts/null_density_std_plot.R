# =========================================================================
# null_density_std_plot.R
# Standardised DSGE null distribution — showing Z-score-like scale collapse
# across pathway sizes, for Supplementary Figure
# =========================================================================
#
# Demonstrates the effect of standardising each size group's null
# distribution as (null - mean(null)) / sd(null).  After standardisation
# all null distributions have mean ≈ 0, sd ≈ 1, making the observed
# DSGE_std values directly comparable across different pathway sizes.
#
# 用法（在项目根目录运行）：
#   source("inst/scripts/null_density_std_plot.R")
#
# 输出：
#   results/plots/null_density_std_curves.pdf   —— 代表性大小的标准化零分布密度
#   results/plots/null_std_thresholds.pdf       —— 标准化后的峰值 & 阈值 vs. 通路大小
#   results/plots/null_density_std_ridge.pdf    —— 标准化脊线图（可选）
# =========================================================================

# ---- 0. 加载依赖 & 读取真实数据 ----
cat("Loading dependencies...\n")
source("R/dsge.R")

cat("Loading limma results...\n")
res <- read.csv("inst/data_exp/limma_FLT3_IR_vs_FLT3.csv", stringsAsFactors = FALSE)
cat("  Genes:", nrow(res), "\n")

# ---- 从真实 limma 结果构建 z 分数池 ----
dsge_res <- calc_dsge(res$pvalue)
pool_z  <- dsge_res$z_scores
n_pool  <- length(pool_z)

cat(sprintf("Background gene pool: n=%d, mean(z)=%.4f, sd(z)=%.4f\n",
            n_pool, mean(pool_z), sd(pool_z)))

# ---- 1. 参数 ----
set.seed(42)
N_PERM     <- 100000L                          # 排列次数
SIZES      <- seq(5, 500, by = 5)             # 通路大小序列
SHOW_SIZES <- c(5, 10, 25, 50, 100, 200, 500) # 面板 A 展示的代表性大小

# ---- 2. 批量生成零分布 → 标准化 ----
peaks_std         <- numeric(length(SIZES))
medians_std       <- numeric(length(SIZES))
thresholds_95_std <- numeric(length(SIZES))
thresholds_99_std <- numeric(length(SIZES))
null_means        <- numeric(length(SIZES))
null_sds          <- numeric(length(SIZES))
null_store        <- vector("list", length(SIZES))   # raw null
null_std_store    <- vector("list", length(SIZES))   # standardised null
names(null_store)     <- as.character(SIZES)
names(null_std_store) <- as.character(SIZES)

cat("\nGenerating & standardising null distributions for", length(SIZES), "pathway sizes...\n")
pb <- utils::txtProgressBar(min = 0, max = length(SIZES), style = 3)

for (i in seq_along(SIZES)) {
  sz  <- SIZES[i]
  bat <- max(1L, floor(n_pool / sz))
  nul <- numeric(N_PERM)

  for (b in seq(1L, N_PERM, by = bat)) {
    nb <- min(bat, N_PERM - b + 1L)
    mat <- matrix(sample.int(n_pool, sz * nb, replace = FALSE), nrow = sz)
    nul[b:(b + nb - 1L)] <- compute_dsge_batch(mat, pool_z)
  }

  null_means[i]    <- mean(nul)
  null_sds[i]      <- sd(nul)
  null_store[[i]]  <- nul

  # Standardise: (null - mean) / sd
  nul_std <- (nul - null_means[i]) / null_sds[i]
  null_std_store[[i]] <- nul_std

  d <- density(nul_std, bw = "SJ")
  peaks_std[i]         <- d$x[which.max(d$y)]
  thresholds_95_std[i] <- as.numeric(quantile(nul_std, 0.95))
  thresholds_99_std[i] <- as.numeric(quantile(nul_std, 0.99))
  medians_std[i]       <- median(nul_std)

  utils::setTxtProgressBar(pb, i)
}
close(pb)
cat("Done.\n")

# ---- 3. 绘图 ----
dir.create("results/plots", showWarnings = FALSE)

# ========== 图 1：代表性大小的标准化零分布密度曲线 ==========
cols <- colorRampPalette(c("#56B4E9", "#0072B2", "#D55E00"))(length(SHOW_SIZES))

density_std <- lapply(as.character(SHOW_SIZES), function(key) {
  density(null_std_store[[key]], bw = "SJ")
})
names(density_std) <- as.character(SHOW_SIZES)
x_range <- range(sapply(density_std, function(d) range(d$x)))
y_max   <- max(sapply(density_std, function(d) max(d$y)))

pdf("results/plots/null_density_std_curves.pdf", width = 7, height = 6)
par(mar = c(4.5, 4.5, 3, 1), mgp = c(2.8, 0.8, 0))

first <- TRUE
for (k in seq_along(SHOW_SIZES)) {
  sz <- SHOW_SIZES[k]
  i  <- which(SIZES == sz)
  d  <- density_std[[as.character(sz)]]

  if (first) {
    plot(d, xlim = x_range, ylim = c(0, y_max * 1.15),
         main = "Standardised Null Distribution  (mean = 0, sd = 1)",
         xlab = "DSGE (standardised)", ylab = "Density",
         col = cols[k], lwd = 2.2, las = 1, yaxt = "n", yaxs = "i")
    first <- FALSE
  } else {
    lines(d, col = cols[k], lwd = 2.2)
  }

  pk   <- peaks_std[i]
  pk_y <- approx(d$x, d$y, xout = pk)$y
  points(pk, pk_y, col = cols[k], pch = 19, cex = 0.9)
}

# Reference line at x = 0
abline(v = 0, col = "#999999", lty = 2, lwd = 0.8)

legend("topright",
       legend = paste0("G = ", SHOW_SIZES),
       col = cols, lwd = 2.2, cex = 0.75, bty = "n", inset = 0.02)

dev.off()
cat("Saved: results/plots/null_density_std_curves.pdf\n")

# ========== 图 2：标准化后的峰值 & 显著性阈值 vs. 通路大小 ==========
pdf("results/plots/null_std_thresholds.pdf", width = 7, height = 6)
par(mar = c(4.5, 4.5, 3, 1), mgp = c(2.8, 0.8, 0))

ylim <- range(c(thresholds_95_std, thresholds_99_std, peaks_std))
plot(SIZES, peaks_std, type = "l", lwd = 2.5, col = "#0072B2",
     xlab = "Pathway size (number of matched genes, G)",
     ylab = "DSGE (standardised)",
     main = "Peak & One-tailed Significance Thresholds  (Standardised)",
     log = "x", las = 1, ylim = ylim)

lines(SIZES, thresholds_95_std, lwd = 2.5, col = "#D55E00", lty = 2)
lines(SIZES, thresholds_99_std, lwd = 1.8, col = "#CC79A7", lty = 3)

x_anno <- 300
text(x_anno, peaks_std[which(SIZES == x_anno)] + 0.03, "Peak (mode)",
     col = "#0072B2", cex = 0.85, pos = 3)
text(x_anno, thresholds_95_std[which(SIZES == x_anno)] + 0.03, "95% threshold",
     col = "#D55E00", cex = 0.85, pos = 3)
text(x_anno, thresholds_99_std[which(SIZES == x_anno)] + 0.03, "99% threshold",
     col = "#CC79A7", cex = 0.85, pos = 1)

# Expected reference line at 0 (standardised peak should be ~0)
abline(h = 0, col = "#999999", lty = 2, lwd = 0.8)
text(min(SIZES), 0, "0", col = "#666666", cex = 0.8, pos = 1)

legend("topright",
       legend = c("Peak (mode)", "95% threshold  (p < 0.05, one-tailed)",
                  "99% threshold  (p < 0.01)", "y = 0"),
       col = c("#0072B2", "#D55E00", "#CC79A7", "#999999"),
       lty = c(1, 2, 3, 2), lwd = c(2.5, 2.5, 1.8, 0.8),
       cex = 0.65, bty = "n", inset = 0.02)

dev.off()
cat("Saved: results/plots/null_std_thresholds.pdf\n")

# ---- 4. 统计摘要（用于论文文字引用） ----
cat("\n")
cat(paste(rep("=", 72), collapse = ""), "\n")
cat("  STANDARDISED NULL SUMMARY  (mean = 0, sd = 1)\n")
cat(paste(rep("=", 72), collapse = ""), "\n")
cat(sprintf("  %-6s %8s %8s %8s %8s\n", "Size", "Peak", "Median", "95%Thr", "99%Thr"))
cat(paste(rep("-", 72), collapse = ""), "\n")
for (sz in SHOW_SIZES) {
  i <- which(SIZES == sz)
  cat(sprintf("  %-6d %8.4f %8.4f %8.4f %8.4f\n",
              sz, peaks_std[i], medians_std[i],
              thresholds_95_std[i], thresholds_99_std[i]))
}
cat(paste(rep("-", 72), collapse = ""), "\n")
cat(sprintf("\n  Expected for standardised null:\n"))
cat(sprintf("    Peak (mode)     ≈ 0    (all sizes)\n"))
cat(sprintf("    95%% threshold   ≈ 1.645  (all sizes)\n"))
cat(sprintf("    99%% threshold   ≈ 2.326  (all sizes)\n"))
cat(sprintf("  Pool mean(z)       = %.4f\n", mean(pool_z)))
cat("\n")

# ---- 5. 脊线图（可选）----
DO_RIDGE <- TRUE
if (DO_RIDGE) {
  pdf("results/plots/null_density_std_ridge.pdf", width = 10, height = 7)

  ridge_sizes <- c(5, 10, 15, 20, 30, 40, 50, 75, 100, 125, 150,
                   175, 200, 250, 300, 350, 400, 450, 500)
  n_ridge <- length(ridge_sizes)
  ridge_col <- colorRampPalette(c("#56B4E9", "#0072B2", "#D55E00"))(n_ridge)

  x_dens <- range(sapply(as.character(ridge_sizes), function(key) {
    range(null_std_store[[key]])
  }))
  x_left_pad <- 0.12 * diff(x_dens)
  x_ridge <- c(x_dens[1] - x_left_pad, x_dens[2])

  par(mar = c(4.5, 5, 3, 1), mgp = c(2.8, 0.8, 0))
  plot(NULL, xlim = x_ridge, ylim = c(0, n_ridge + 1),
       xlab = "DSGE (standardised)", ylab = "", yaxt = "n",
       main = "Standardised Null Ridgeline  —  Pathway Size 5 to 500")

  for (k in seq_along(ridge_sizes)) {
    sz  <- ridge_sizes[k]
    nul <- null_std_store[[as.character(sz)]]
    d   <- density(nul, bw = "SJ")
    y_scaled <- d$y / max(d$y) * 1.2

    polygon(c(d$x, rev(d$x)),
            c(y_scaled + k, rep(k, length(y_scaled))),
            col = adjustcolor(ridge_col[k], alpha.f = 0.6),
            border = ridge_col[k])

    text(x_dens[1], k, paste0("G=", sz),
         cex = 0.55, col = "#333333", pos = 2)
  }

  abline(v = 0, col = "#999999", lty = 2, lwd = 0.8)

  dev.off()
  cat("Saved: results/plots/null_density_std_ridge.pdf\n")
}

cat("\nAll plots generated.\n")
