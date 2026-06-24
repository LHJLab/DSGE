# =========================================================================
# gpd_diagnostic_plot.R
# GPD fit diagnostic plot -- visualises tail GPD fit quality of DSGE null
# distribution, for validating the reliability of extreme value extrapolation
# =========================================================================
#
# Usage (run from project root):
#   source("inst/scripts/gpd_diagnostic_plot.R")
#
# Output:
#   plots/gpd_diagnostic_Gm{GM}.pdf  —— PDF (for publication)
#   plots/gpd_diagnostic_Gm{GM}.jpg  —— JPG (quick preview)
# =========================================================================

# ---- Output directory (use tempdir to avoid writing to user workspace) ----
out_dir <- file.path(tempdir(), "plots")
dir.create(out_dir, showWarnings = FALSE)

# ---- Save user graphical parameters (restore at script exit) ----
old_par <- par(no.readonly = TRUE)
on.exit(par(old_par))

# ---- 0. Parameters ----
GM      <- 50          # Pathway size
N_PERM  <- 100000L     # Number of permutations
TAU     <- 0.99        # GPD fit threshold quantile
SEED    <- 666

# ---- 1. Load dependencies & real data ----
cat("Loading dependencies...\n")
source("R/dsge.R")

cat("Loading DESeq2 results...\n")
res <- read.csv("inst/data_exp/limma_FLT3_IR_vs_FLT3.csv",
                stringsAsFactors = FALSE)
cat("  Genes loaded:", nrow(res), "\n")

# ---- 2. Build z-score pool ----
dsge_res <- calc_dsge(res$pvalue)
pool_z   <- dsge_res$z_scores
n_pool   <- length(pool_z)

cat(sprintf("Background gene pool: n=%d, mean(z)=%.3f, sd(z)=%.3f\n",
            n_pool, mean(pool_z), sd(pool_z)))
cat(sprintf("max(z)=%.2f, p(z>2)=%.1f%%\n",
            max(pool_z), mean(pool_z > 2) * 100))

# ---- 3. Generate null distribution ----
set.seed(SEED)
cat(sprintf("Generating null (Gm=%d, n_perm=%d)...\n", GM, N_PERM))

bat <- max(1L, floor(n_pool / GM))
nul <- numeric(N_PERM)

for (b in seq(1L, N_PERM, by = bat)) {
  nb <- min(bat, N_PERM - b + 1L)
  mat <- matrix(sample.int(n_pool, GM * nb, replace = FALSE), nrow = GM)
  nul[b:(b + nb - 1L)] <- compute_dsge_batch(mat, pool_z)
}

cat(sprintf("Done.\n  mean(null)=%.3f, sd(null)=%.3f\n", mean(nul), sd(nul)))

# ---- 4. GPD fitting ----
gpd_fit <- fit_gpd_tail(nul, tail = TAU)

if (is.null(gpd_fit)) {
  stop("GPD fit returned NULL — cannot generate diagnostic plot.")
}

n_excess <- sum(nul > gpd_fit$u)
cat(sprintf("GPD fit at tau=%.2f: scale=%.4f, shape=%.4f\n",
            TAU, gpd_fit$scale, gpd_fit$shape))
cat(sprintf("  u=%.3f (%dth pctl), pat=%.4f\n",
            gpd_fit$u, as.integer(TAU * 100), gpd_fit$pat))

# ---- 5. Prepare Q-Q plot data ----
# Tail excess: values exceeding the threshold
excess <- nul[nul > gpd_fit$u] - gpd_fit$u
excess_sort <- sort(excess)
n_excess <- length(excess_sort)

# GPD theoretical quantiles
gpd_quant <- evd::qgpd(ppoints(n_excess),
                        scale = gpd_fit$scale,
                        shape = gpd_fit$shape)
# Remove NaN (qgpd returns NaN when xi < 0 and quantile > theoretical upper bound)
valid <- is.finite(gpd_quant) & is.finite(excess_sort)
gpd_quant   <- gpd_quant[valid]
excess_sort <- excess_sort[valid]

# ---- 6. Plotting function (shared between PDF and JPG) ----
draw_diagnostic <- function() {
  layout(matrix(1:2, nrow = 1), widths = c(1.7, 1.3))

  # ========== Panel A: Tail Zoom — GPD Fit Highlight ==========
  par(mar = c(4.5, 4.5, 3, 1), mgp = c(2.8, 0.8, 0))

  x_min <- gpd_fit$u - 0.5
  x_max <- max(nul) * 1.08

  # Plot full data with xlim clipping to match density scale with GPD curve
  h <- hist(nul, breaks = 200, plot = FALSE)
  hist(nul, breaks = 200, col = "#F0F0F0", border = "#D0D0D0",
       xlim = c(x_min, x_max),
       ylim = c(0, max(h$density[h$mids >= x_min], na.rm = TRUE) * 1.8),
       main = sprintf("A  GPD Tail Fit (%dth percentile threshold)",
                      as.integer(TAU * 100)),
       xlab = expression(DSGE[std]),
       ylab = "Density",
       las = 1, yaxs = "i", xaxs = "i",
       freq = FALSE)

  # Threshold vertical line
  abline(v = gpd_fit$u, col = "#2166AC", lwd = 1.8, lty = 2)

  # Threshold label
  text(gpd_fit$u, par("usr")[4] * 0.95,
       sprintf("u = %.2f", gpd_fit$u),
       col = "#2166AC", cex = 0.9, pos = 2)

  # GPD fitted curve (from u to right end)
  x_gpd <- seq(gpd_fit$u, x_max, length.out = 500)
  y_gpd <- gpd_fit$pat * evd::dgpd(x_gpd - gpd_fit$u,
                                    scale = gpd_fit$scale,
                                    shape = gpd_fit$shape)

  # Fill GPD fit area
  polygon(c(x_gpd, rev(x_gpd)),
          c(y_gpd, rep(0, length(y_gpd))),
          col = adjustcolor("#D6604D", alpha.f = 0.25),
          border = NA)

  # GPD fit curve outline
  lines(x_gpd, y_gpd, col = "#D6604D", lwd = 2.5)

  # GPD parameter annotation
  if (gpd_fit$shape < 0) {
    upper_bound <- gpd_fit$u - gpd_fit$scale / gpd_fit$shape
    bound_line <- sprintf("upper bound = %.2f", upper_bound)
  } else {
    bound_line <- "heavy-tailed (no finite upper bound)"
  }

  legend("topright", bty = "n", inset = c(0.02, 0.02),
         legend = c(
           sprintf("GPD:  scale = %.3f,  shape = %.3f",
                   gpd_fit$scale, gpd_fit$shape),
           bound_line,
           sprintf("tail obs. = %d  (%.1f%%)", n_excess, gpd_fit$pat * 100)
         ),
         text.col = c("#D6604D", "#666666", "#666666"),
         cex = 0.8)

  # ========== Panel B: Q-Q Plot ==========
  par(mar = c(4.5, 4.5, 3, 1), mgp = c(2.8, 0.8, 0))

  qq_lim <- range(c(gpd_quant, excess_sort))
  plot(gpd_quant, excess_sort, pch = 19, cex = 0.4, col = "#BB0000",
       xlim = qq_lim, ylim = qq_lim,
       main = "B  GPD Q-Q Plot",
       xlab = "Theoretical quantiles (GPD)",
       ylab = "Empirical quantiles (tail excess)",
       las = 1, xaxs = "i", yaxs = "i")

  abline(0, 1, col = "#333333", lwd = 1.5, lty = 2)

  # Annotate fit quality
  r_sq <- round(cor(gpd_quant, excess_sort)^2, 4)
  legend("bottomright",
         legend = c("y = x", sprintf("R² = %.4f", r_sq)),
         col = c("#333333", NA),
         lty = c(2, NA),
         lwd = c(1.5, NA),
         cex = 0.8, bty = "n", inset = c(0.02, 0.02),
         text.col = c("#333333", "#444444"))
}

# ---- 7. Output PDF + JPG ----
# PDF (write to temp file first to avoid Windows file lock)
pdf_tmp <- tempfile(fileext = ".pdf")
pdf(pdf_tmp, width = 10, height = 5.5)
draw_diagnostic()
dev.off()
file.copy(pdf_tmp, file.path(out_dir, sprintf("gpd_diagnostic_Gm%d.pdf", GM)),
          overwrite = TRUE)
unlink(pdf_tmp)
cat(sprintf("Saved: %s/gpd_diagnostic_Gm%d.pdf\n", out_dir, GM))

jpeg(file.path(out_dir, sprintf("gpd_diagnostic_Gm%d.jpg", GM)),
     width = 10, height = 5.5, units = "in", res = 150, quality = 95)
draw_diagnostic()
dev.off()
cat(sprintf("Saved: %s/gpd_diagnostic_Gm%d.jpg\n", out_dir, GM))

cat("Done.\n")
