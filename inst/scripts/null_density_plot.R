# =========================================================================
# null_density_plot.R
# DSGE null distribution density plot -- illustrates how the null
# distribution evolves across pathway sizes, for Supplementary Figure
# =========================================================================
#
# Uses the z-score pool from real DESeq2 data to generate null
# distributions at various pathway sizes, records peak and one-tailed
# significance thresholds, and produces density curves and summary plots.
#
# Usage (run from project root):
#   source("inst/scripts/null_density_plot.R")
#
# Output:
#   plots/null_density_curves.pdf  —— density curves
#   plots/null_thresholds.pdf      —— threshold plots
#   plots/null_density_ridge.pdf   —— ridge plot (optional)
# =========================================================================

# ---- Output directory (use tempdir to avoid writing to user workspace) ----
out_dir <- file.path(tempdir(), "plots")
dir.create(out_dir, showWarnings = FALSE)

# ---- Save user graphical parameters (restore at script exit) ----
old_par <- par(no.readonly = TRUE)
on.exit(par(old_par))

cat("Loading dependencies...\n")
source("R/read_gaf.R")
source("R/read_obo.R")
source("R/get_pathway_genes.R")
source("R/dsge.R")
#library(DSGEr)
cat("Loading DESeq2 results...\n")
res <- read.csv("inst/data_exp/limma_FLT3_IR_vs_FLT3.csv", stringsAsFactors = FALSE)
#res <- subset(res, geneType == "protein_coding" & geneName != ".")
cat("  Protein-coding genes:", nrow(res), "\n")

# ---- Build z-score pool from real DESeq2 results ----
dsge_res <- calc_dsge(res$pvalue)
pool_z  <- dsge_res$z_scores
n_pool  <- length(pool_z)

cat(sprintf("Background gene pool: n=%d, mean(z)=%.4f, sd(z)=%.4f\n",
            n_pool, mean(pool_z), sd(pool_z)))

# ---- 1. Parameters ----
set.seed(666)
N_PERM     <- 100000L                          # Number of permutations
SIZES      <- seq(5, 500, by = 5)             # Pathway size sequence
SHOW_SIZES <- c(5, 10, 25, 50, 100, 200, 500) # Representative sizes for Panel A

# ---- 2. Batch generate null distributions (compute once, reuse) ----
peaks         <- numeric(length(SIZES))
thresholds_95 <- numeric(length(SIZES))
thresholds_99 <- numeric(length(SIZES))
medians       <- numeric(length(SIZES))
null_means    <- numeric(length(SIZES))
null_sds      <- numeric(length(SIZES))
null_store    <- vector("list", length(SIZES))  # Store null vectors for re-use in plotting
names(null_store) <- as.character(SIZES)

cat("\nGenerating null distributions for", length(SIZES), "pathway sizes...\n")
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

  d <- density(nul, bw = "SJ")
  peaks[i]         <- d$x[which.max(d$y)]
  thresholds_95[i] <- as.numeric(quantile(nul, 0.95))
  thresholds_99[i] <- as.numeric(quantile(nul, 0.99))
  medians[i]       <- median(nul)
  null_means[i]    <- mean(nul)
  null_sds[i]      <- sd(nul)
  null_store[[i]]  <- nul

  utils::setTxtProgressBar(pb, i)
}
close(pb)
cat("Done.\n")

# ---- 3. Plotting ----

# ========== Fig 1: Null density curves for representative sizes ==========
cols <- colorRampPalette(c("#56B4E9", "#0072B2", "#D55E00"))(length(SHOW_SIZES))

density_info <- lapply(as.character(SHOW_SIZES), function(key) {
  density(null_store[[key]], bw = "SJ")
})
names(density_info) <- as.character(SHOW_SIZES)
x_range <- range(sapply(density_info, function(d) range(d$x)))
y_max   <- max(sapply(density_info, function(d) max(d$y)))

pdf(file.path(out_dir, "null_density_curves.pdf"), width = 7, height = 6)
par(mar = c(4.5, 4.5, 3, 1), mgp = c(2.8, 0.8, 0))

first <- TRUE
for (k in seq_along(SHOW_SIZES)) {
  sz <- SHOW_SIZES[k]
  i  <- which(SIZES == sz)
  d  <- density_info[[as.character(sz)]]

  if (first) {
    plot(d, xlim = x_range, ylim = c(0, y_max * 1.15),
         main = "Null Distribution by Pathway Size",
         xlab = "DSGE (mean z-score)", ylab = "Density",
         col = cols[k], lwd = 2.2, las = 1, yaxt = "n", yaxs = "i")
    first <- FALSE
  } else {
    lines(d, col = cols[k], lwd = 2.2)
  }

  pk   <- peaks[i]
  pk_y <- approx(d$x, d$y, xout = pk)$y
  points(pk, pk_y, col = cols[k], pch = 19, cex = 0.9)
}

legend("topright",
       legend = paste0("G = ", SHOW_SIZES),
       col = cols, lwd = 2.2, cex = 0.75, bty = "n", inset = 0.02)

dev.off()
cat(sprintf("Saved: %s/null_density_curves.pdf\n", out_dir))

# ========== Fig 2: Peak & significance threshold vs. pathway size ==========
pdf(file.path(out_dir, "null_thresholds.pdf"), width = 7, height = 6)
par(mar = c(4.5, 4.5, 3, 1), mgp = c(2.8, 0.8, 0))

ylim <- range(c(thresholds_95, thresholds_99, peaks))
plot(SIZES, peaks, type = "l", lwd = 2.5, col = "#0072B2",
     xlab = "Pathway size (number of matched genes, G)",
     ylab = "DSGE  (mean z-score)",
     main = "Peak & One-tailed Significance Thresholds",
     log = "x", las = 1, ylim = ylim)

lines(SIZES, thresholds_95, lwd = 2.5, col = "#D55E00", lty = 2)
lines(SIZES, thresholds_99, lwd = 1.8, col = "#CC79A7", lty = 3)

x_anno <- 300
text(x_anno, peaks[which(SIZES == x_anno)] + 0.02, "Peak (mode)",
     col = "#0072B2", cex = 0.85, pos = 3)
text(x_anno, thresholds_95[which(SIZES == x_anno)] + 0.02, "95% threshold",
     col = "#D55E00", cex = 0.85, pos = 3)
text(x_anno, thresholds_99[which(SIZES == x_anno)] + 0.02, "99% threshold",
     col = "#CC79A7", cex = 0.85, pos = 1)


abline(h = mean(pool_z), col = "#999999", lty = 2, lwd = 0.8)
text(min(SIZES), mean(pool_z), expression(bar(z)[pool]), col = "#666666",
     cex = 0.8, pos = 1)

legend("topright",
       legend = c("Peak (mode)", "95% threshold (p < 0.05, one-tailed)",
                  "99% threshold (p < 0.01)", "Pool mean"),
       col = c("#0072B2", "#D55E00", "#CC79A7", "#999999"),
       lty = c(1, 2, 3, 2), lwd = c(2.5, 2.5, 1.8, 0.8),
       cex = 0.65, bty = "n", inset = 0.02)

dev.off()
cat(sprintf("Saved: %s/null_thresholds.pdf\n", out_dir))

# ---- 4. Summary statistics table (for manuscript text references) ----
cat("\n")
cat(paste(rep("=", 72), collapse = ""), "\n")
cat("  NULL DISTRIBUTION SUMMARY  (based on real DESeq2 data)\n")
cat(paste(rep("=", 72), collapse = ""), "\n")
cat(sprintf("  %-6s %8s %8s %8s %8s\n", "Size", "Peak", "Median", "95%Thr", "99%Thr"))
cat(paste(rep("-", 72), collapse = ""), "\n")
for (sz in SHOW_SIZES) {
  i <- which(SIZES == sz)
  cat(sprintf("  %-6d %8.4f %8.4f %8.4f %8.4f\n",
              sz, peaks[i], medians[i], thresholds_95[i], thresholds_99[i]))
}
cat(paste(rep("-", 72), collapse = ""), "\n")

cat(sprintf("\n  Pool mean(z)            = %.4f\n", mean(pool_z)))
cat(sprintf("  Peak at G = 500         = %.4f  (deviation: %.6f)\n",
            peaks[length(peaks)], peaks[length(peaks)] - mean(pool_z)))
cat(sprintf("  Peak at G = 5           = %.4f\n", peaks[1]))
cat(sprintf("  95%% threshold at G=5     = %.4f\n", thresholds_95[1]))
cat(sprintf("  95%% threshold at G=500   = %.4f\n", thresholds_95[length(thresholds_95)]))
cat(sprintf("  Threshold gap at G=5    = %.4f  (95%% - peak)\n",
            thresholds_95[1] - peaks[1]))
cat(sprintf("  Threshold gap at G=500  = %.4f  (95%% - peak)\n",
            thresholds_95[length(thresholds_95)] - peaks[length(thresholds_95)]))
cat("\n")

# ---- 5. Ridge plot (optional) ----
DO_RIDGE <- TRUE
if (DO_RIDGE) {
  pdf(file.path(out_dir, "null_density_ridge.pdf"), width = 10, height = 7)

  ridge_sizes <- c(5, 10, 15, 20, 30, 40, 50, 75, 100, 125, 150,
                   175, 200, 250, 300, 350, 400, 450, 500)
  n_ridge <- length(ridge_sizes)
  ridge_col <- colorRampPalette(c("#56B4E9", "#0072B2", "#D55E00"))(n_ridge)

  # Use stored null vectors to determine x-axis range; leave space on left for labels
  x_dens <- range(sapply(as.character(ridge_sizes), function(key) {
    range(null_store[[key]])
  }))
  x_left_pad <- 0.12 * diff(x_dens)
  x_ridge <- c(x_dens[1] - x_left_pad, x_dens[2])

  par(mar = c(4.5, 5, 3, 1), mgp = c(2.8, 0.8, 0))
  plot(NULL, xlim = x_ridge, ylim = c(0, n_ridge + 1),
       xlab = "DSGE (mean z-score)", ylab = "", yaxt = "n",
       main = "Null Distribution Ridgeline — Pathway Size 5 to 500")

  for (k in seq_along(ridge_sizes)) {
    sz  <- ridge_sizes[k]
    nul <- null_store[[as.character(sz)]]
    d   <- density(nul, bw = "SJ")
    y_scaled <- d$y / max(d$y) * 1.2

    polygon(c(d$x, rev(d$x)),
            c(y_scaled + k, rep(k, length(y_scaled))),
            col = adjustcolor(ridge_col[k], alpha.f = 0.6),
            border = ridge_col[k])

    text(x_dens[1], k, paste0("G=", sz),
         cex = 0.55, col = "#333333", pos = 2)
  }

  dev.off()
  cat(sprintf("Saved: %s/null_density_ridge.pdf\n", out_dir))
}

cat("\nAll plots generated.\n")
