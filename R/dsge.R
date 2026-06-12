# =========================================================================
# DSGE: Disruption Score of Gene Expression
# Pathway-level transcriptional perturbation analysis
# =========================================================================
#
# Core idea:
#   1. Convert DESeq2 differential expression p-values to absolute z-scores
#      (|z-score|); larger z means stronger transcriptional perturbation.
#   2. For each GO pathway (gene set), compute the mean of member gene
#      z-scores as the pathway DSGE.
#' @useDynLib DSGE, .registration = TRUE
#' @importFrom Rcpp evalCpp
#   3. Generate null distributions via permutation test: randomly draw
#      equally-sized gene sets from the background pool, compute random
#      DSGE, repeat n times.
#   4. For extreme observations (above the 90th percentile), fit a
#      Generalized Pareto Distribution (GPD) to the tail of the null
#      distribution for extrapolation, yielding precise extreme p-values.
#   5. Apply Benjamini-Hochberg FDR correction across all pathways.
#
# Key formulas (computed on each drawn gene subset):
#   Per-gene z-score:    z_i = |Φ⁻¹(1 - p_i/2)|
#   Pathway DSGE:        DSGE = mean(z_i)
#
# Dependencies:
#   evd  - pgpd(): survival function of the Generalized Pareto Distribution
#   POT  - fitgpd(): fit GPD to the tail of the null distribution
#
# References:
#   Coles, S. (2001). An Introduction to Statistical Modeling of
#     Extreme Values. Springer.
# =========================================================================

# ---- Internal helper: p-values to absolute z-scores ----
#
# Convert two-sided test p-values to absolute standard normal z-scores.
# p=0 (numerical underflow) is replaced with machine epsilon to avoid
# infinite z-scores.
#
# Note: This function only computes raw z_i = |Φ⁻¹(1 - p_i/2)|.
#
# Args:
#   pvalue - numeric vector, DESeq2 p-values
#
# Returns: vector of absolute z-scores
compute_zscore <- function(pvalue) {
  # Convert two-sided p-value to absolute z-score
  # qnorm(p/2, lower.tail = FALSE) = |Φ⁻¹(1 - p/2)|
  z <- abs(qnorm(pvalue / 2, lower.tail = FALSE))

  # Handle p=0 (DESeq2 may produce p=0 due to numerical precision);
  # qnorm returns Inf in this case; replace with machine epsilon
  # .Machine$double.xmin ≈ 2.2e-308
  z[!is.finite(z)] <- abs(qnorm(.Machine$double.xmin / 2, lower.tail = FALSE))
  z
}


# ---- Internal helper: compute DSGE on a gene subset ----
#
# Given a gene index subset, compute DSGE = mean(z_i).
#
# Args:
#   idx    - indices of genes in pool_z
#   pool_z - raw z-scores of the full gene pool
#   pool_z - raw z-scores of the full gene pool
#
# Returns: scalar DSGE
compute_dsge <- function(idx, pool_z) {
  mean(pool_z[idx])
}


# ---- Internal helper: batch-compute DSGE on a subset matrix ----
#
# Used in pathway_dsge for vectorized null distribution generation.
# Computes DSGE for each column of a sz × nb matrix (one permutation
# sample per column).
#
# Args:
#   mat    - sz × nb integer matrix, each column is one permutation sample
#   pool_z - raw z-scores of the full gene pool
#
# Returns: numeric vector of length nb (one DSGE per column)
compute_dsge_batch <- function(mat, pool_z) {
  val_z <- pool_z[mat]; dim(val_z) <- dim(mat)
  colMeans(val_z)
}


# ---- Internal helper: fit GPD to the upper tail of a null distribution ----
#
# Fits a Generalized Pareto Distribution to the upper tail of a
# permutation-generated null distribution. GPD is the core distribution
# of extreme value theory, used to extrapolate p-values beyond the
# permutation range, avoiding p=0 issues from direct counting.
#
# GPD cumulative distribution function (excess above threshold u):
#   H(y) = 1 - [1 + ξ*y/σ]^(-1/ξ),  y = x - u > 0
#   where ξ (shape) is the shape parameter, σ (scale) the scale parameter
#
#   ξ < 0: distribution has finite upper bound u - σ/ξ; p=0 when obs exceeds it
#   ξ = 0: reduces to exponential distribution
#   ξ > 0: heavy-tailed distribution, suitable for biological extremes
#
# Args:
#   null - numeric vector, null distribution (permutation DSGE values)
#   tail - quantile threshold, default 0.90 (fit GPD above 90th percentile)
#
# Returns:
#   On success: list(u, scale, shape, pat); on failure: NULL
#   u     - threshold (90th percentile)
#   scale - GPD scale parameter σ
#   shape - GPD shape parameter ξ
#   pat   - empirical proportion of null exceeding threshold, P(X > u)
#
# Note: support-constrained adjustment (arXiv:2602.22975) is applied
# downstream in eval_gpd_p(), not at fitting time, because GPD is
# fitted per size group while the constraint is pathway-specific.
fit_gpd_tail <- function(null, tail = 0.99, gpd_method = "mle") {
  u <- unname(quantile(null, tail))

  # Too few samples above threshold (< 10): tail fit unreliable, skip GPD
  if (sum(null > u) < 10) return(NULL)

  # Fit GPD using POT::fitgpd (warnings suppressed: benign S3 method
  # overwrite messages from evd/POT namespace loading in parallel workers)
  # Estimation methods: "mle" (default), "mple", "moments", "pwmu", "pwmb",
  # "mdpd", "med", "pickands", "lme", "mgf"
  fit <- tryCatch(
    suppressWarnings(POT::fitgpd(null, u, est = gpd_method)),
    error = function(e) {
      warning("GPD fit failed: ", e$message, call. = FALSE)
      NULL
    }
  )
  if (is.null(fit) || is.null(fit$param)) return(NULL)

  list(u = u, scale = fit$param[1], shape = fit$param[2], pat = fit$pat)
}


# ---- Internal helper: compute p-value from fitted GPD ----
#
# Given an observed DSGE and fitted GPD parameters, compute the
# right-tail p-value:
#   p = P(X > obs) = P(X > u) * P(X > obs | X > u)
#                  = pat * pgpd(obs - u, lower.tail = FALSE)
#
# pgpd is the GPD survival function from the evd package, using the
# same Coles (2001) parametrization as POT::fitgpd, guaranteeing
# compatibility.
#
# When ξ < 0 and obs exceeds the theoretical GPD upper bound,
# a support-constrained adjustment (arXiv:2602.22975) is applied,
# ensuring a non-zero p-value.
#
# Args:
#   fit - GPD parameter list from fit_gpd_tail()
#   obs - observed DSGE value
#
# Returns: right-tail p-value (range [0, 1])
#' Evaluate GPD-based p-value with support-constrained adjustment
#'
#' Computes the right-tail p-value for an observed DSGE value using
#' a fitted Generalized Pareto Distribution. When the unconstrained GPD
#' would give p = 0 due to a finite upper bound (\eqn{\xi < 0}), a
#' support-constrained adjustment is applied (arXiv:2602.22975) to yield
#' a small but non-zero p-value.
#'
#' @param fit GPD parameter list from \code{\link{fit_gpd_tail}()}.
#'   Must contain elements \code{u}, \code{scale}, \code{shape}, \code{pat}.
#' @param obs Observed DSGE value.
#' @param safety_margin Safety margin for support-constrained adjustment.
#'   Default \code{1.2}. Larger values produce more conservative p-values.
#'
#' @return Right-tail p-value in \eqn{[0, 1]}.
#' @noRd
eval_gpd_p <- function(fit, obs, safety_margin = 1.2) {
  # pat = P(X > u), pgpd(*, lower.tail = FALSE) = P(X > obs | X > u)
  p <- fit$pat * evd::pgpd(obs - fit$u, scale  = fit$scale,
                            shape = fit$shape, lower.tail = FALSE)

  # ---- Support-constrained adjustment ----
  # When pgpd returns 0 (finite upper bound with ξ < 0, or numerical
  # underflow for extreme observations under any ξ), we adjust the shape
  # parameter ξ upward so the new upper bound sits at obs + (safety_margin
  # - 1) × (obs - u), guaranteeing a non-zero p-value.
  # Larger safety_margin → larger adjusted p-value (more conservative).
  # safety_margin = 1.0 disables the margin (tight bound at obs).
  if (isTRUE(p == 0)) {
    shape_adj <- -fit$scale / ((obs - fit$u) * safety_margin)
    if (isTRUE(shape_adj < 0)) {
      p <- fit$pat * evd::pgpd(obs - fit$u, scale = fit$scale,
                                shape = shape_adj, lower.tail = FALSE)
    }
  }

  pmin(pmax(p, 0), 1)
}


# ---- Internal helper: compute Gini coefficient of a gene set ----
#
# Measures the unevenness of gene z-score distribution within a pathway.
# Gini ∈ [0, 1]: 0 = perfectly uniform perturbation across all genes,
# 1 = extreme inequality (single gene dominates).
#
# Formula (after sorting): G = (2·Σ(i·xᵢ)) / (n·Σxᵢ) - (n+1)/n
#
# Args:
#   x - non-negative numeric vector (gene z-scores)
#
# Returns: scalar Gini coefficient
compute_gini <- function(x) {
  n <- length(x)
  if (n < 2L) return(0)
  s <- sum(x)
  if (s == 0) return(0)
  x_sorted <- sort(x)
  (2 * sum(seq_len(n) * x_sorted)) / (n * s) - (n + 1) / n
}


# ---- Internal helper: compute coefficient of variation of a gene set ----
#
# CV = sd(x) / mean(x), measures relative dispersion.
# For non-negative z-scores, larger CV indicates more uneven perturbation.
#
# Args:
#   x - non-negative numeric vector
#
# Returns: scalar CV
compute_cv <- function(x) {
  m <- mean(x)
  if (m == 0) return(NA_real_)
  stats::sd(x) / m
}


# ---- Internal helper: compute Normalised Direction Score (NDS) ----
#
# Computes NDS = (U - D) / max(U, D) for a single pathway, where U is
# the sum of z-scores of up-regulated genes and D is the sum of
# z-scores of down-regulated genes.
#
# By default, only the top 25% most-perturbed genes (ranked by absolute
# direction value, e.g. |log2FC|) are retained; the bottom 75% is
# treated as noise.  Within this subset, U is the sum of z-scores of
# all up-regulated (dir > 0) genes and D the sum of all down-regulated
# (dir < 0) genes.  If the top-25% subset contains fewer than 10 genes
# the function falls back to using all genes — the legacy behaviour.
#
# Args:
#   z_sub - numeric vector of absolute z-scores for pathway genes
#   d_sub - numeric vector of direction values (e.g. log2FC) of the
#           same length; NAs are silently removed
#
# Returns:
#   Numeric scalar in [-1, 1].  Positive → net up-regulation;
#   negative → net down-regulation; NA when no valid direction info.
compute_nds <- function(z_sub, d_sub) {
  # Drop genes with missing direction information
  valid <- !is.na(d_sub)
  z_sub <- z_sub[valid]
  d_sub <- d_sub[valid]

  n <- length(z_sub)
  if (n == 0L) return(NA_real_)

  n_quarter <- floor(n * 0.25)

  if (n_quarter >= 10L) {
    # Keep only the top 25% most-perturbed genes (largest |log2FC|)
    ord    <- order(abs(d_sub), decreasing = TRUE)
    keep   <- ord[seq_len(n_quarter)]
    z_keep <- z_sub[keep]
    d_keep <- d_sub[keep]
    up     <- z_keep[d_keep > 0]
    down   <- z_keep[d_keep < 0]
    if (length(up) == 0L || length(down) == 0L) {
      if (length(up) > 0L)  return(1)
      if (length(down) > 0L) return(-1)
      return(NA_real_)
    }
    u <- sum(up)
    d <- sum(down)
  } else {
    # Fallback: use all directional genes (legacy behaviour)
    up   <- z_sub[d_sub > 0]
    down <- z_sub[d_sub < 0]
    if (length(up) == 0L || length(down) == 0L) {
      if (length(up) > 0L)  return(1)
      if (length(down) > 0L) return(-1)
      return(NA_real_)
    }
    u <- sum(up)
    d <- sum(down)
  }

  if (u == 0 && d == 0) return(NA_real_)
  max_ud <- max(u, d)
  if (max_ud == 0) return(NA_real_)
  (u - d) / max_ud
}


# ---- Internal helper: batch-compute Gini coefficient ----
#
# Computes Gini for each column of a sz × nb matrix (one permutation
# sample per column).
#
# Args:
#   mat    - sz × nb integer matrix
#   pool_z - raw z-scores of the full gene pool
#
# Returns: numeric vector of length nb
compute_gini_batch <- function(mat, pool_z) {
  val_z <- pool_z[mat]
  dim(val_z) <- dim(mat)
  apply(val_z, 2, compute_gini)
}


# ---- Internal helper: batch-compute coefficient of variation ----
#
# Computes CV for each column of a sz × nb matrix.
#
# Args:
#   mat    - sz × nb integer matrix
#   pool_z - raw z-scores of the full gene pool
#
# Returns: numeric vector of length nb
compute_cv_batch <- function(mat, pool_z) {
  val_z <- pool_z[mat]
  dim(val_z) <- dim(mat)
  apply(val_z, 2, compute_cv)
}


# =========================================================================
# Exported function 1: calc_dsge -- genome-wide DSGE
# =========================================================================

#' Compute DSGE (Disruption Score of Gene Expression)
#'
#' Computes the mean absolute z-score for all genes passing filters
#' from differential expression analysis results. Higher DSGE indicates
#' stronger global transcriptional perturbation.
#'
#' @param pvalue Numeric vector of p-values from differential expression
#'   analysis (e.g., DESeq2, edgeR, Seurat FindMarkers).
#' @param base_mean Numeric vector of mean expression values
#'   (same length as pvalue), e.g. DESeq2 baseMean or Seurat avg_log2FC
#'   corresponding expression level. If NULL, no expression filtering
#'   is applied.
#' @param base_mean_cutoff Expression filter threshold, default 0.1.
#'   Genes with mean expression at or below this value are excluded as
#'   lowly expressed. Ignored when base_mean = NULL.
#'
#' @return A list with elements:
#'   \item{dsge}{scalar, genome-wide DSGE}
#'   \item{n_genes}{integer, number of genes passing filters}
#'   \item{z_scores}{named numeric vector of per-gene raw z-scores}
#' @importFrom stats qnorm quantile sd setNames
#' @export
#'
#' @examples
#' \dontrun{
#' res <- DESeq2::results(dds)
#' calc_dsge(res$pvalue, res$baseMean)
#' }
calc_dsge <- function(pvalue, base_mean = NULL, base_mean_cutoff = 0.1) {
  # ---- Input validation ----
  if (length(pvalue) == 0)
    stop("'pvalue' is empty", call. = FALSE)
  if (!is.null(base_mean) && length(base_mean) != length(pvalue))
    stop("'base_mean' must have the same length as 'pvalue'", call. = FALSE)

  # ---- Gene filtering ----
  keep <- !is.na(pvalue)
  if (!is.null(base_mean))
    keep <- keep & !is.na(base_mean) & base_mean > base_mean_cutoff
  if (sum(keep) == 0)
    stop("No genes pass the baseMean > ", base_mean_cutoff, " filter", call. = FALSE)

  # ---- Compute raw z-scores ----
  pool_z <- compute_zscore(pvalue[keep])
  if (!is.null(names(pvalue)))
    names(pool_z) <- names(pvalue)[keep]

  # ---- Compute DSGE ----
  dsge_val <- compute_dsge(seq_along(pool_z), pool_z)

  list(dsge = dsge_val, n_genes = sum(keep), z_scores = pool_z)
}


# =========================================================================
# Exported function 2: dsge_perm_test -- permutation test for a single gene set
# =========================================================================

#' Permutation test for DSGE enrichment of a single gene set
#'
#' Tests whether a target gene set has a significantly higher DSGE
#' than expected by chance. Generates a null distribution via sampling
#' without replacement and computes an empirical right-tail p-value.
#'
#' @details
#' Algorithm steps:
#' \enumerate{
#'   \item Filter gene pool: baseMean > cutoff, non-missing p-value.
#'   \item Convert p-values to per-gene raw z-scores.
#'   \item Match \code{gene_list} to the filtered gene pool; compute
#'         observed DSGE (mean z-score).
#'   \item Generate null distribution: randomly sample (without
#'         replacement) equally-sized gene sets from the pool,
#'         computing random DSGE via the same mean z-score formula, repeated
#'         n_perm times.
#'   \item Empirical right-tail p-value: count(DSGE_null > DSGE_obs) / n_perm.
#' }
#'
#' @param gene_list Character vector of target gene identifiers
#'   (must match values in gene_names).
#' @param pvalue Numeric vector of p-values from differential expression
#'   analysis (e.g., DESeq2, edgeR, Seurat FindMarkers) for all genes
#'   in the background pool.
#' @param base_mean Numeric vector of mean expression values,
#'   same length as pvalue. Pass NULL to skip expression-level filtering.
#' @param gene_names Character vector of gene names, same length as
#'   pvalue, must be unique.
#' @param base_mean_cutoff baseMean filter threshold, default 0.1.
#' @param n_perm Number of permutations, default 10000.
#' @param seed Optional integer random seed for reproducibility.
#' @param progress Whether to show a progress bar, default TRUE.
#' @param heterogeneity Whether to compute perturbation heterogeneity
#'   indices (Gini, CV) and two-sided p-values. Default \code{FALSE}.
#'   When enabled, also computes Gini and CV null distributions within
#'   the permutation loop, increasing runtime by approximately
#'   30\%-50\%.
#' @param directional Whether to compute Normalized Direction Score
#'   (NDS) for each pathway. Default \code{FALSE}. When \code{TRUE},
#'   \code{direction_vec} must be provided. NDS ranges from -1 to 1:
#'   positive values indicate net up-regulation, negative values net
#'   down-regulation. Requires \code{use_gpd = TRUE} (default) for
#'   reliable extreme p-values alongside directional information.
#' @param direction_vec Numeric vector of direction indicators
#'   (e.g., log fold changes), same length as \code{pvalue}. Used to
#'   classify genes as up-regulated (positive direction) or
#'   down-regulated (negative direction). Values of exactly 0 are
#'   treated as up-regulated. Ignored when \code{directional = FALSE}.
#' @param use_std Whether to compute and use standardised DSGE.
#'   Default \code{TRUE}. When enabled:
#'   \itemize{
#'     \item The returned list includes \code{dsge_std = (observed -
#'           mean(null)) / sd(null)}.
#'     \item p-value is computed from the standardised null
#'           distribution (\code{mean(null_std >= dsge_std)}).
#'   }
#'   When \code{FALSE}, p-value is computed from the raw null
#'   distribution.
#' @param use_gpd Whether to use GPD tail extrapolation for extreme
#'   p-values. Default \code{TRUE}. When \code{TRUE} and the observed
#'   DSGE exceeds the \code{gpd_threshold} quantile of the null, the
#'   p-value is extrapolated via a fitted Generalized Pareto
#'   Distribution with support-constrained adjustment (Peschel et al.
#'   2025, arXiv:2602.22975) ensuring a non-zero p-value. When
#'   \code{FALSE}, always uses empirical ECDF
#'   (p-value always >= \code{1/n_perm}).
#' @param gpd_threshold Tail quantile threshold for GPD fitting,
#'   between 0 and 1. Default 0.99. Lower = more tail samples (lower
#'   variance, higher bias); higher = fewer samples (higher variance,
#'   lower bias).
#' @param gpd_method GPD estimation method passed to
#'   \code{POT::fitgpd}. Default \code{"mle"} (maximum likelihood).
#'   Other options: \code{"mple"}, \code{"moments"}, \code{"pwmu"},
#'   \code{"pwmb"}, \code{"mdpd"}, \code{"med"}, \code{"pickands"},
#'   \code{"lme"}, \code{"mgf"}.
#' @param safety_margin Safety margin for GPD support-constrained adjustment.
#'   Default \code{1.05} (5 % margin). Larger values (e.g., \code{2}) produce
#'   larger (more conservative) p-values for extremely extreme observations,
#'   avoiding double-precision underflow to zero at the cost of increased bias.
#'
#' @return A list with elements:
#'   \item{observed}{observed DSGE value}
#'   \item{n_genes}{number of target genes matched}
#'   \item{null}{permutation null distribution vector}
#'   \item{p_value}{right-tail p-value (GPD tail extrapolation when
#'         observed falls above 90th percentile; empirical ECDF otherwise)}
#'   \item{ecdf}{empirical cumulative distribution function of the null}
#'   \item{dsge_std}{(only when \code{use_std = TRUE}) standardised DSGE}
#'   \item{nds}{(only when \code{directional = TRUE}) Normalized
#'         Direction Score, ranging from -1 (pure down) to +1 (pure up)}
#'   \item{gini_observed}{(only when \code{heterogeneity = TRUE}) observed Gini coefficient}
#'   \item{cv_observed}{(only when \code{heterogeneity = TRUE}) observed CV}
#'   \item{null_gini}{(only when \code{heterogeneity = TRUE}) Gini null distribution vector}
#'   \item{null_cv}{(only when \code{heterogeneity = TRUE}) CV null distribution vector}
#'   \item{het_p_value}{(only when \code{heterogeneity = TRUE}) two-sided permutation p-value based on Gini}
#' @export
#'
#' @examples
#' \dontrun{
#' res <- DESeq2::results(dds)
#' dsge_perm_test(gene_list = forebrain_genes,
#'                pvalue = res$pvalue, base_mean = res$baseMean,
#'                gene_names = rownames(res), seed = 42)
#' # Without standardised DSGE
#' dsge_perm_test(gene_list = forebrain_genes,
#'                pvalue = res$pvalue, base_mean = res$baseMean,
#'                gene_names = rownames(res), seed = 42, use_std = FALSE)
#' }
dsge_perm_test <- function(gene_list, pvalue, base_mean, gene_names,
                           base_mean_cutoff = 0.1,
                           n_perm           = 10000L,
                           seed             = NULL,
                           progress         = TRUE,
                           heterogeneity    = FALSE,
                           directional      = FALSE,
                           direction_vec    = NULL,
                           use_std           = TRUE,
                           use_gpd           = TRUE,
                           gpd_threshold     = 0.99,
                           gpd_method        = "mle",
                           safety_margin     = 1.05) {
  # ---- Input validation ----
  stopifnot(is.character(gene_list), length(gene_list) > 0)
  if (!is.null(base_mean))
    stopifnot(length(pvalue) == length(base_mean))
  stopifnot(length(pvalue) == length(gene_names))
  if (anyDuplicated(gene_names))
    stop("'gene_names' must be unique", call. = FALSE)
  if (isTRUE(directional)) {
    stopifnot("'direction_vec' must be provided when directional = TRUE" =
                !is.null(direction_vec))
    stopifnot(length(direction_vec) == length(pvalue))
  }

  # ---- Build filtered gene pool ----
  keep <- !is.na(pvalue)
  if (!is.null(base_mean))
    keep <- keep & !is.na(base_mean) & base_mean > base_mean_cutoff
  pool_z <- compute_zscore(pvalue[keep])
  names(pool_z) <- gene_names[keep]
  n_pool <- length(pool_z)
  if (n_pool == 0) {
    if (!is.null(base_mean))
      stop("No genes pass the baseMean > ", base_mean_cutoff, " filter", call. = FALSE)
    else
      stop("No genes pass the non-missing pvalue filter", call. = FALSE)
  }

  # ---- Direction pool (optional) ----
  if (isTRUE(directional)) {
    pool_dir_raw <- direction_vec[keep]                # raw log2FC for top-25%
    pool_dir     <- sign(pool_dir_raw)
    pool_dir[pool_dir == 0] <- 1
    keep_dir <- !is.na(pool_dir)
    if (sum(keep_dir) == 0)
      stop("No genes with non-NA direction after filtering", call. = FALSE)
    pool_z       <- pool_z[keep_dir]
    pool_dir     <- pool_dir[keep_dir]
    pool_dir_raw <- pool_dir_raw[keep_dir]
    n_pool       <- length(pool_z)
  }

  # ---- Match target genes to gene pool ----
  idx <- match(gene_list, names(pool_z))
  idx <- idx[!is.na(idx)]
  if (length(idx) == 0)
    stop("No genes from 'gene_list' found in the filtered gene pool", call. = FALSE)

  n_target <- length(idx)
  observed <- compute_dsge(idx, pool_z)

  if (heterogeneity) {
    observed_gini <- compute_gini(pool_z[idx])
    observed_cv   <- compute_cv(pool_z[idx])
  }
  if (isTRUE(directional)) {
    z_sub <- pool_z[idx]
    d_sub <- pool_dir_raw[idx]                        # raw log2FC for abs-ranking
    nds_observed <- compute_nds(z_sub, d_sub)
  }

  # ---- Permutation: generate null distribution ----
  if (!is.null(seed)) set.seed(seed)
  null <- numeric(n_perm)
  if (heterogeneity) {
    null_gini <- numeric(n_perm)
    null_cv   <- numeric(n_perm)
  }

  if (progress) pb <- utils::txtProgressBar(min = 0, max = n_perm, style = 3)
  for (i in seq_len(n_perm)) {
    samp <- sample.int(n_pool, n_target, replace = FALSE)
    null[i] <- compute_dsge(samp, pool_z)
    if (heterogeneity) {
      z_sub <- pool_z[samp]
      null_gini[i] <- compute_gini(z_sub)
      null_cv[i]   <- compute_cv(z_sub)
    }
    if (progress && i %% max(1L, floor(n_perm / 100)) == 0)
      utils::setTxtProgressBar(pb, i)
  }
  if (progress) { utils::setTxtProgressBar(pb, n_perm); close(pb) }

  # ---- Fit GPD (optional), compute standardised DSGE & p-value ----
  if (isTRUE(use_std)) {
    out <- list(observed  = observed,
                n_genes   = n_target,
                null      = null,
                dsge_std  = (observed - mean(null)) / sd(null),
                nds       = if (isTRUE(directional)) nds_observed else NULL,
                ecdf      = stats::ecdf(null))

    null_std <- (null - mean(null)) / sd(null)
    if (isTRUE(use_gpd)) {
      gpd_std  <- fit_gpd_tail(null_std, tail = gpd_threshold, gpd_method = gpd_method)
      if (!is.null(gpd_std) && out$dsge_std > gpd_std$u) {
        out$p_value <- eval_gpd_p(gpd_std, out$dsge_std, safety_margin = safety_margin)
      } else {
        out$p_value <- sum(null_std >= out$dsge_std) / n_perm
      }
    } else {
      out$p_value <- sum(null_std >= out$dsge_std) / n_perm
    }
  } else {
    out <- list(observed = observed,
                n_genes  = n_target,
                null     = null,
                nds      = if (isTRUE(directional)) nds_observed else NULL,
                ecdf     = stats::ecdf(null))

    if (isTRUE(use_gpd)) {
      gpd <- fit_gpd_tail(null, tail = gpd_threshold, gpd_method = gpd_method)
      if (!is.null(gpd) && observed > gpd$u) {
        out$p_value <- eval_gpd_p(gpd, observed, safety_margin = safety_margin)
      } else {
        out$p_value <- sum(null > observed) / n_perm
      }
    } else {
      out$p_value <- sum(null > observed) / n_perm
    }
  }

  if (heterogeneity) {
    p_left  <- sum(null_gini <= observed_gini) / n_perm
    p_right <- sum(null_gini >= observed_gini) / n_perm
    out$gini_observed <- observed_gini
    out$cv_observed   <- observed_cv
    out$null_gini     <- null_gini
    out$null_cv       <- null_cv
    out$het_p_value   <- 2 * min(p_left, p_right)
  }

  out
}


# =========================================================================
# Exported function 3: pathway_dsge -- batch pathway DSGE with
#   size-grouped permutation + FDR
# =========================================================================

#' Pathway-level DSGE permutation test with FDR correction
#'
#' Computes DSGE for each GO pathway, generates permutation null
#' distributions grouped by number of matched genes, computes p-values
#' using GPD tail extrapolation, and applies Benjamini-Hochberg FDR
#' correction.
#'
#' Pathways sharing the same number of matched genes reuse the same
#' null distribution, greatly reducing computation.
#'
#' @details
#' **Algorithm steps (9 steps)**
#' \enumerate{
#'   \item \strong{Build gene pool}: filter DESeq2 results
#'         (baseMean > cutoff), compute per-gene raw z-scores.
#'   \item \strong{Match pathway genes}: match each pathway's genes to
#'         the gene pool via the \code{gene_id_col} column.
#'   \item \strong{Filter pathways}: retain pathways with
#'         \code{min_size <= n_matched <= max_size}.
#'   \item \strong{Compute observed DSGE}: per pathway as the mean of
#'         member gene z-scores: DSGE = mean(z_i).
#'   \item \strong{Generate size-grouped null distributions}: each unique
#'         matched-gene count gets its own permutation run (vectorized
#'         batch sampling), with
#'         simultaneous GPD tail fitting.
#'   \item \strong{Compute standardised DSGE} (step 6): when
#'         \code{use_std = TRUE}, compute \code{dsge_std = (observed -
#'         mean(null)) / sd(null)} for each pathway, using the
#'         size-grouped null distribution.
#'   \item \strong{Compute p-values}: when \code{use_std = TRUE},
#'         empirical ECDF of the standardised null vs. dsge_std; when
#'         \code{use_std = FALSE}, GPD tail extrapolation (above 90th
#'         percentile) + empirical ECDF on the raw null distribution.
#'   \item \strong{BH FDR correction}: Benjamini-Hochberg multiple
#'         testing correction on all pathway p-values.
#'   \item \strong{Sort and return}: ordered by p_adj ascending.
#' }
#'
#' @param pathway_genes Named list returned by
#'   \code{\link{get_pathway_genes}()}. Each element is a data.frame
#'   containing gene information for one pathway.
#' @param pvalue Numeric vector of p-values from differential expression
#'   analysis (e.g., DESeq2, edgeR, Seurat FindMarkers) for all genes
#'   in the background pool.
#' @param base_mean Numeric vector of mean expression values,
#'   same length as pvalue. Pass NULL to skip expression-level filtering.
#' @param gene_names Character vector of gene names (same length as
#'   pvalue), must be unique.
#' @param gene_id_col Column name in pathway_genes data.frames used to
#'   match gene names, default \code{"db_object_symbol"}.
#' @param base_mean_cutoff baseMean filter threshold, default 0.1.
#' @param min_size Minimum number of matched genes; pathways below this
#'   are not tested. Default 5.
#' @param max_size Maximum number of matched genes; pathways above this
#'   are excluded. Default 500. Set to Inf to retain all.
#' @param n_perm Number of permutations per size group, default 10000.
#' @param seed Optional integer random seed.
#' @param return_null Whether to include null distribution data in the
#'   result (for plotting). Default \code{FALSE}. When \code{TRUE},
#'   returns a list instead of a data.frame.
#' @param progress Whether to show a progress bar, default TRUE.
#' @param heterogeneity Whether to compute perturbation heterogeneity
#'   indices (Gini, CV) and two-sided p-values. Default \code{FALSE}.
#'   When enabled, additionally computes Gini and CV during null
#'   distribution generation, increasing runtime by approximately
#'   30\%-50\%. The result data.frame will include extra columns
#'   \code{gini}, \code{cv}, \code{het_p_value}, \code{het_p_adj}.
#' @param directional Whether to compute Normalized Direction Score
#'   (NDS) for each pathway. Default \code{FALSE}. When \code{TRUE},
#'   \code{direction_vec} must be provided. NDS ranges from -1 to 1:
#'   positive values indicate net up-regulation, negative values net
#'   down-regulation. The result data.frame will include an
#'   \code{nds} column.
#' @param direction_vec Numeric vector of direction indicators
#'   (e.g., log fold changes), same length as \code{pvalue}. Used to
#'   classify genes as up-regulated (positive direction) or
#'   down-regulated (negative direction). Values of exactly 0 are
#'   treated as up-regulated. Ignored when \code{directional = FALSE}.
#' @param use_std Whether to compute and use standardised DSGE.
#'   Default \code{TRUE}. When enabled:
#'   \itemize{
#'     \item \code{dsge_std = (observed - mean(null)) / sd(null)}
#'           using the size-grouped null distribution is added as a
#'           result column.
#'     \item p-values are computed from the standardised null
#'           distribution (\code{mean(null_std >= dsge_std)}),
#'           without GPD tail extrapolation.
#'   }
#'   When \code{FALSE}, p-values are computed from the raw null
#'   distribution (using GPD when \code{use_gpd = TRUE}).
#' @param use_gpd Whether to use GPD tail extrapolation for extreme
#'   p-values. Default \code{TRUE}. When \code{TRUE} and the observed
#'   DSGE exceeds the \code{gpd_threshold} quantile of the null, the
#'   p-value is extrapolated via a fitted Generalized Pareto
#'   Distribution with support-constrained adjustment (Peschel et al.
#'   2025, arXiv:2602.22975) ensuring a non-zero p-value. When
#'   \code{FALSE}, always uses empirical ECDF
#'   (p-value always >= \code{1/n_perm}).
#' @param gpd_threshold Tail quantile threshold for GPD fitting,
#'   between 0 and 1. Default 0.99. Lower = more tail samples (lower
#'   variance, higher bias); higher = fewer samples (higher variance,
#'   lower bias).
#' @param gpd_method GPD estimation method passed to
#'   \code{POT::fitgpd}. Default \code{"mle"} (maximum likelihood).
#'   Other options: \code{"mple"}, \code{"moments"}, \code{"pwmu"},
#'   \code{"pwmb"}, \code{"mdpd"}, \code{"med"}, \code{"pickands"},
#'   \code{"lme"}, \code{"mgf"}.
#' @param safety_margin Safety margin for GPD support-constrained adjustment.
#'   Default \code{1.05} (5 % margin). Larger values (e.g., \code{2}) produce
#'   larger (more conservative) p-values for extremely extreme observations,
#'   avoiding double-precision underflow to zero at the cost of increased bias.
#' @param n_cores Number of CPU cores for parallel null distribution
#'   generation. Default \code{1} (sequential). Set to
#'   \code{parallel::detectCores()} to use all available cores.
#'   Uses \code{parallel::mclapply} on Unix/macOS (only effective
#'   when multiple unique pathway sizes exist). On Windows, falls
#'   back to sequential with a message.
#' @param dsge_std \strong{[Deprecated]}. Use \code{use_std} instead.
#' @param p_adjust_method Multiple testing correction method. Passed to
#'   \code{stats::p.adjust()}. Default \code{"BY"} (Benjamini-Yekutieli)
#'   which controls FDR under arbitrary dependence. Use \code{"BH"} for
#'   Benjamini-Hochberg (controls FDR under positive dependence).
#'   All methods supported by \code{p.adjust} are valid:
#'   \code{"holm"}, \code{"hochberg"}, \code{"hommel"},
#'   \code{"bonferroni"}, \code{"BH"}, \code{"BY"}, \code{"fdr"},
#'   \code{"none"}.
#'
#' @return By default, a \code{data.frame} sorted by \code{p_adj}
#'   ascending, with columns:
#'   \code{go_id}, \code{go_name}, \code{aspect}, \code{n_pathway},
#'   \code{n_matched}, \code{dsge}, \code{p_value}, \code{p_adj}.
#'   \code{aspect} is the GO ontology classification:
#'   \code{"BP"} (Biological Process), \code{"MF"} (Molecular Function),
#'   \code{"CC"} (Cellular Component). If \code{get_pathway_genes()} was
#'   called without \code{go_names}, the \code{aspect} column is empty.
#'
#'   When \code{heterogeneity = TRUE}, additionally includes:
#'   \code{gini}, \code{cv}, \code{het_p_value}, \code{het_p_adj}.
#'
#'   When \code{directional = TRUE}, additionally includes:
#'   \code{nds} (Normalized Direction Score, ranging from -1 to +1).
#'
#'   When \code{use_std = TRUE} (default), additionally includes:
#'   \code{dsge_std} = (observed DSGE - mean(null)) / sd(null),
#'   standardised using the size-grouped null distribution.
#'
#'   When \code{return_null = TRUE}, returns a \code{list} with:
#'   \item{table}{the data.frame described above}
#'   \item{null_raw}{named list, keys are pathway sizes (as character),
#'         values are DSGE null distribution vectors}
#'   \item{null_gpd}{named list, keys are pathway sizes,
#'         values are GPD fit parameters (or NULL)}
#'
#'   When both \code{heterogeneity = TRUE} and \code{return_null = TRUE},
#'   the list also includes:
#'   \item{null_gini_raw}{named list, keys are pathway sizes,
#'         values are Gini null distribution vectors}
#'   \item{null_cv_raw}{named list, keys are pathway sizes,
#'         values are CV null distribution vectors}
#' @export
#'
#' @examples
#' \dontrun{
#' res <- DESeq2::results(dds)
#' gaf <- read_gaf("goa_human.gaf")
#' go  <- read_obo("go.obo")
#'
#' pw <- get_pathway_genes(gaf, go_names = go, min_size = 5)
#' result <- pathway_dsge(pw, res$pvalue, res$baseMean, rownames(res),
#'                        seed = 42)
#' head(result)
#' }
pathway_dsge <- function(pathway_genes, pvalue, base_mean = NULL, gene_names,
                         gene_id_col       = "db_object_symbol",
                         base_mean_cutoff  = 0.1,
                         min_size          = 5L,
                         max_size          = 500L,
                         n_perm            = 10000L,
                         seed              = NULL,
                         return_null       = FALSE,
                         progress          = TRUE,
                         heterogeneity     = FALSE,
                         directional       = FALSE,
                         direction_vec     = NULL,
                         use_std           = TRUE,
                         use_gpd           = TRUE,
                         gpd_threshold     = 0.99,
                         gpd_method        = "mle",
                         safety_margin     = 1.05,
                         n_cores           = 1L,
                         dsge_std          = NULL,
                         p_adjust_method   = "BY") {
  # ---- Backward compatibility: dsge_std → use_std ----
  if (!is.null(dsge_std)) {
    warning("'dsge_std' is deprecated; use 'use_std' instead", call. = FALSE)
    use_std <- dsge_std
  }

  # ---- Input validation ----
  stopifnot(is.list(pathway_genes), length(pathway_genes) > 0)
  if (!is.null(base_mean))
    stopifnot(length(pvalue) == length(base_mean))
  stopifnot(length(pvalue) == length(gene_names))
  stopifnot(is.character(gene_id_col), length(gene_id_col) == 1L)
  stopifnot(min_size >= 1L, n_perm >= 1L)
  p_adjust_method <- match.arg(p_adjust_method,
                               c("holm", "hochberg", "hommel",
                                 "bonferroni", "BH", "BY", "fdr", "none"))
  if (isTRUE(directional)) {
    stopifnot("'direction_vec' must be provided when directional = TRUE" =
                !is.null(direction_vec))
    stopifnot(length(direction_vec) == length(pvalue))
  }

  # =========================================================================
  # Step 1: Build DESeq2 gene pool
  # =========================================================================
  # Handle duplicate gene names: keep first occurrence, drop rest
  if (anyDuplicated(gene_names)) {
    warning(sum(duplicated(gene_names)),
            " duplicate 'gene_names' found; keeping first occurrence",
            call. = FALSE)
    uniq <- !duplicated(gene_names)
    pvalue    <- pvalue[uniq]
    if (!is.null(base_mean))
      base_mean <- base_mean[uniq]
    if (isTRUE(directional))
      direction_vec <- direction_vec[uniq]
    gene_names <- gene_names[uniq]
  }

  # Filter: p-value non-missing, baseMean non-missing, baseMean > cutoff
  keep <- !is.na(pvalue)
  if (!is.null(base_mean))
    keep <- keep & !is.na(base_mean) & base_mean > base_mean_cutoff
  pool_z <- compute_zscore(pvalue[keep])        # raw z-scores
  names(pool_z) <- gene_names[keep]
  n_pool <- length(pool_z)
  if (n_pool == 0) {
    if (!is.null(base_mean))
      stop("No genes pass the baseMean > ", base_mean_cutoff, " filter", call. = FALSE)
    else
      stop("No genes pass the non-missing pvalue filter", call. = FALSE)
  }

  # ---- Direction pool (optional) ----
  if (isTRUE(directional)) {
    pool_dir_raw <- direction_vec[keep]                # raw log2FC for top-25%
    pool_dir     <- sign(pool_dir_raw)
    # sign(0) = 0; treat as up (NA/0 directions removed from pool)
    pool_dir[pool_dir == 0] <- 1
    keep_dir <- !is.na(pool_dir)
    if (sum(keep_dir) == 0)
      stop("No genes with non-NA direction after filtering", call. = FALSE)
    pool_z       <- pool_z[keep_dir]
    pool_dir     <- pool_dir[keep_dir]
    pool_dir_raw <- pool_dir_raw[keep_dir]
    n_pool       <- length(pool_z)
    if (isTRUE(progress))
      cat("  directional: after filtering", sum(!keep_dir),
          "genes with NA/zero direction\n")
  }

  # =========================================================================
  # Step 2: Match each pathway's genes to the gene pool
  # =========================================================================
  go_ids    <- names(pathway_genes)
  n_pathway <- vapply(pathway_genes, nrow, integer(1L))
  go_name   <- character(length(go_ids))
  go_aspect <- character(length(go_ids))
  matched   <- vector("list", length(go_ids))

  for (i in seq_along(pathway_genes)) {
    df <- pathway_genes[[i]]
    g  <- if (gene_id_col %in% names(df)) df[[gene_id_col]] else character(0L)
    if ("go_name" %in% names(df) && nrow(df) > 0)
      go_name[i] <- df$go_name[1]
    if ("go_namespace" %in% names(df) && nrow(df) > 0)
      go_aspect[i] <- df$go_namespace[1]
    matched[[i]] <- match(g, names(pool_z), nomatch = 0L)
    matched[[i]] <- matched[[i]][matched[[i]] > 0L]
  }

  # =========================================================================
  # Step 3: Filter pathways by min_size and max_size
  # =========================================================================
  n_matched <- lengths(matched)
  keep_pw   <- n_matched >= min_size & n_matched <= max_size
  if (sum(keep_pw) == 0)
    stop("No pathways with ", min_size, " <= n_matched <= ", max_size, call. = FALSE)

  go_ids    <- go_ids[keep_pw]
  go_name   <- go_name[keep_pw]
  go_aspect <- go_aspect[keep_pw]
  n_pathway <- n_pathway[keep_pw]
  n_matched <- n_matched[keep_pw]
  matched   <- matched[keep_pw]

  # =========================================================================
  # Step 4: Compute observed DSGE per pathway (and optional Gini, CV, NDS)
  # =========================================================================
  # DSGE = mean(z_i)
  observed <- vapply(matched, function(idx) compute_dsge(idx, pool_z),
                     numeric(1L))
  if (heterogeneity) {
    gini_obs <- vapply(matched, function(idx) compute_gini(pool_z[idx]),
                       numeric(1L))
    cv_obs   <- vapply(matched, function(idx) compute_cv(pool_z[idx]),
                       numeric(1L))
  }
  if (isTRUE(directional)) {
    nds_obs <- vapply(matched, function(idx) {
      compute_nds(pool_z[idx], pool_dir_raw[idx])   # raw log2FC for abs-ranking
    }, numeric(1L))
  }

  # =========================================================================
  # Step 5: Generate size-grouped null distributions + optional GPD tail fitting
  # =========================================================================
  # Pathways sharing the same matched-gene count reuse a single null.
  # Batch vectorized sampling: draw bat sets per iteration, compute
  # all at once with compute_dsge_batch.
  #
  # When use_gpd = TRUE, a Generalized Pareto Distribution is fitted to
  # the upper tail of each null distribution for extreme-value p-value
  # extrapolation. When use_gpd = FALSE, p-values are computed purely
  # from the empirical ECDF (with 1/n_perm floor).
  #
  sizes    <- unique(n_matched)
  null_raw <- vector("list", length(sizes))
  names(null_raw) <- as.character(sizes)
  if (isTRUE(use_gpd)) {
    null_gpd <- vector("list", length(sizes))
    names(null_gpd) <- as.character(sizes)
    if (isTRUE(use_std)) {
      null_gpd_std <- vector("list", length(sizes))
      names(null_gpd_std) <- as.character(sizes)
    }
  }
  if (heterogeneity) {
    null_gini_raw <- vector("list", length(sizes))
    null_cv_raw   <- vector("list", length(sizes))
    names(null_gini_raw) <- names(null_cv_raw) <- as.character(sizes)
  }

  n_cores_effective <- if (n_cores > 1L) {
    if (.Platform$OS.type == "windows") {
      if (progress) cat("Note: n_cores > 1 is not supported on Windows;",
                        "running sequentially.\n")
      1L
    } else {
      min(n_cores, length(sizes))
    }
  } else 1L

  if (!is.null(seed)) set.seed(seed)

  # ---- Debug: print reproducible fingerprint of the null-generation state ----
  if (isTRUE(progress)) {
    fingerprint <- sum(n_pool * 31 + as.integer(sort(sizes)) * 7 + n_perm)
    cat("  perm fingerprint:", fingerprint,
        "(n_pool =", n_pool, ", n_sizes =", length(sizes), ")\n")
  }

  if (n_cores_effective > 1L) {
    # ---- Parallel null generation ----
    if (progress)
      cat("Generating null distributions for", length(sizes),
          "unique sizes using", n_cores_effective, "cores\n")

    seed_base <- if (is.null(seed)) -1L else seed
    results <- suppressWarnings(
      parallel::mclapply(seq_along(sizes), function(s) {
      suppressMessages(library(POT))   # forestall S3 registration noise
      sz  <- sizes[s]
      res <- permute_null_cpp(pool_z, sz, n_perm,
                              seed_base + s,
                              compute_gini = heterogeneity,
                              compute_cv   = heterogeneity)
      gpd_obj <- gpd_std_obj <- NULL
      if (isTRUE(use_gpd)) {
        gpd_obj <- fit_gpd_tail(res$null, tail = gpd_threshold, gpd_method = gpd_method)
        if (isTRUE(use_std)) {
          nul_std <- (res$null - mean(res$null)) / sd(res$null)
          gpd_std_obj <- fit_gpd_tail(nul_std, tail = gpd_threshold, gpd_method = gpd_method)
        }
      }
      list(null = res$null, null_gpd = gpd_obj, null_gpd_std = gpd_std_obj,
           null_gini = if (heterogeneity) res$null_gini else NULL,
           null_cv   = if (heterogeneity) res$null_cv   else NULL)
    }, mc.cores = n_cores_effective, mc.preschedule = TRUE)
    )  # suppressWarnings

    for (s in seq_along(sizes)) {
      null_raw[[s]] <- results[[s]]$null
      if (isTRUE(use_gpd)) {
        null_gpd[[s]] <- results[[s]]$null_gpd
        if (isTRUE(use_std)) null_gpd_std[[s]] <- results[[s]]$null_gpd_std
      }
      if (heterogeneity) {
        null_gini_raw[[s]] <- results[[s]]$null_gini
        null_cv_raw[[s]]   <- results[[s]]$null_cv
      }
    }

  } else {
    # ---- Sequential null generation ----
    if (progress) {
      cat("Generating null distributions for", length(sizes),
          "unique pathway sizes\n")
      pb <- utils::txtProgressBar(min = 0, max = length(sizes), style = 3)
    }

    seed_base <- if (is.null(seed)) -1L else seed

    for (s in seq_along(sizes)) {
      sz  <- sizes[s]
      res <- permute_null_cpp(pool_z, sz, n_perm,
                              seed_base + s,
                              compute_gini = heterogeneity,
                              compute_cv   = heterogeneity)

      null_raw[[s]] <- res$null
      if (isTRUE(use_gpd)) {
        null_gpd[[s]] <- fit_gpd_tail(res$null, tail = gpd_threshold, gpd_method = gpd_method)
        if (isTRUE(use_std)) {
          nul_std <- (res$null - mean(res$null)) / sd(res$null)
          null_gpd_std[[s]] <- fit_gpd_tail(nul_std, tail = gpd_threshold, gpd_method = gpd_method)
        }
      }
      if (heterogeneity) {
        null_gini_raw[[s]] <- res$null_gini
        null_cv_raw[[s]]   <- res$null_cv
      }
      if (progress) utils::setTxtProgressBar(pb, s)
    }
    if (progress) { utils::setTxtProgressBar(pb, length(sizes)); close(pb) }
  }

  # =========================================================================
  # Step 6: Standardised DSGE (must precede p-value; use_std = TRUE reuses it)
  # =========================================================================
  if (isTRUE(use_std)) {
    dsge_std_vals <- numeric(length(observed))
    for (i in seq_along(observed)) {
      key  <- as.character(n_matched[i])
      nul  <- null_raw[[key]]
      dsge_std_vals[i] <- (observed[i] - mean(nul)) / sd(nul)
    }
  }

  # =========================================================================
  # Step 7: Compute p-value per pathway
  # =========================================================================
  # Two p-value methods:
  #   use_gpd = TRUE  — GPD tail extrapolation when observed > 90th %ile u,
  #                     with support-constrained adjustment
  #                     (arXiv:2602.22975) to avoid p = 0; empirical ECDF
  #                     otherwise (with 1/n_perm floor)
  #   use_gpd = FALSE — pure empirical ECDF always (with 1/n_perm floor)
  # In both cases, the standardised null is used when use_std = TRUE.
  #
  p_val <- numeric(length(observed))
  for (i in seq_along(observed)) {
    key <- as.character(n_matched[i])

    if (isTRUE(use_std)) {
      nul <- null_raw[[key]]
      nul_std <- (nul - mean(nul)) / sd(nul)
      if (isTRUE(use_gpd)) {
        gpd <- null_gpd_std[[key]]
        if (!is.null(gpd) && dsge_std_vals[i] > gpd$u) {
          p_val[i] <- eval_gpd_p(gpd, dsge_std_vals[i], safety_margin = safety_margin)
        } else {
          p_val[i] <- sum(nul_std >= dsge_std_vals[i]) / n_perm
          if (p_val[i] == 0) p_val[i] <- 1 / n_perm
        }
      } else {
        p_val[i] <- sum(nul_std >= dsge_std_vals[i]) / n_perm
        if (p_val[i] == 0) p_val[i] <- 1 / n_perm
      }
    } else {
      if (isTRUE(use_gpd)) {
        gpd <- null_gpd[[key]]
        if (!is.null(gpd) && observed[i] > gpd$u) {
          p_val[i] <- eval_gpd_p(gpd, observed[i], safety_margin = safety_margin)
        } else {
          p_val[i] <- sum(null_raw[[key]] > observed[i]) / n_perm
          if (p_val[i] == 0) p_val[i] <- 1 / n_perm
        }
      } else {
        p_val[i] <- sum(null_raw[[key]] > observed[i]) / n_perm
        if (p_val[i] == 0) p_val[i] <- 1 / n_perm
      }
    }
  }

  # ---- Heterogeneity two-sided p-values (Gini-based empirical ECDF, no GPD) ----
  # Two-sided test: left tail = uniform perturbation, right tail = heterogeneous
  if (heterogeneity) {
    het_p <- numeric(length(gini_obs))
    for (i in seq_along(gini_obs)) {
      key  <- as.character(n_matched[i])
      ng   <- null_gini_raw[[key]]
      p_left  <- sum(ng <= gini_obs[i]) / n_perm
      p_right <- sum(ng >= gini_obs[i]) / n_perm
      het_p[i] <- 2 * min(p_left, p_right)
    }
  }

  # =========================================================================
  # Step 8: Multiple testing correction
  # =========================================================================
  p_adj <- stats::p.adjust(p_val, method = p_adjust_method)
  if (heterogeneity)
    het_p_adj <- stats::p.adjust(het_p, method = p_adjust_method)

  # =========================================================================
  # Step 9: Assemble results, sorted by p_adj ascending
  # =========================================================================
  result <- data.frame(
    go_id      = go_ids,
    go_name    = go_name,
    aspect     = go_aspect,
    n_pathway  = as.integer(n_pathway),
    n_matched  = as.integer(n_matched),
    dsge       = observed,
    p_value    = p_val,
    p_adj      = p_adj,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  if (isTRUE(use_std)) {
    # Insert dsge_std immediately after dsge column
    dsge_pos <- which(names(result) == "dsge")
    result <- data.frame(
      result[, seq_len(dsge_pos), drop = FALSE],
      dsge_std = dsge_std_vals,
      result[, seq(dsge_pos + 1, ncol(result)), drop = FALSE],
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }
  if (isTRUE(directional)) {
    # Insert nds after dsge (or dsge_std), before p_value
    nds_pos <- if (isTRUE(use_std)) which(names(result) == "dsge_std") + 1L
               else which(names(result) == "dsge") + 1L
    result <- data.frame(
      result[, seq_len(nds_pos - 1L), drop = FALSE],
      nds = nds_obs,
      result[, seq(nds_pos, ncol(result)), drop = FALSE],
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }
  if (heterogeneity) {
    result$gini        <- gini_obs
    result$cv          <- cv_obs
    result$het_p_value <- het_p
    result$het_p_adj   <- het_p_adj
  }
  result <- result[order(result$p_adj), ]
  rownames(result) <- NULL

  # If null distribution data is not needed, return data.frame directly
  if (!isTRUE(return_null)) return(result)

  # Return list with table and null distributions (for plot_dsge())
  out <- list(table = result, null_raw = null_raw,
              safety_margin = safety_margin)
  if (isTRUE(use_gpd)) {
    out$null_gpd <- null_gpd
    if (isTRUE(use_std)) out$null_gpd_std <- null_gpd_std
  }
  if (heterogeneity) {
    out$null_gini_raw <- null_gini_raw
    out$null_cv_raw   <- null_cv_raw
  }
  out
}


# =========================================================================
# Exported function 4: plot_dsge -- null distribution density plots
# =========================================================================

#' Plot null distribution vs. observed DSGE for selected pathways
#'
#' For pathways in a \code{\link{pathway_dsge}()} result, plots the
#' density of the permutation null distribution with a dashed red line
#' marking the observed DSGE. If GPD tail fitting is available for the
#' size group, the tail region is highlighted in semi-transparent orange.
#'
#' @param result List returned by \code{\link{pathway_dsge}()} with
#'   \code{return_null = TRUE} (contains \code{table}, \code{null_raw},
#'   \code{null_gpd}).
#' @param n When \code{go_ids = NULL}, the top \code{n} most significant
#'   pathways (by p_adj ascending) are plotted. Default 9. Ignored when
#'   \code{go_ids} is specified.
#' @param go_ids Optional character vector of GO term IDs to plot
#'   (e.g. \code{c("GO:0007156", "GO:0007268")}). When provided,
#'   directly plots these pathways, overriding \code{n}. Unmatched IDs
#'   are skipped with a warning.
#' @param col_null Color of the null distribution density curve.
#'   Default \code{"steelblue"}.
#' @param col_obs Color of the observed DSGE vertical line.
#'   Default \code{"red"}.
#' @param col_tail Color of the GPD tail region highlight.
#'   Default \code{"orange"} (semi-transparent).
#' @param cex_main Title font scaling. Default \code{0.85}.
#' @param use_std Whether to plot standardised DSGE. Defaults to
#'   \code{TRUE} when the result table contains a \code{dsge_std}
#'   column (i.e. \code{pathway_dsge()} was run with
#'   \code{use_std = TRUE}), and \code{FALSE} otherwise.
#'
#' @return No return value; called for its side effect (plotting).
#' @importFrom graphics abline legend lines par rect title
#' @importFrom stats density na.omit sd
#' @export
#'
#' @examples
#' \dontrun{
#' result <- pathway_dsge(pw, res$pvalue, res$baseMean, res$geneName,
#'                        seed = 42, return_null = TRUE, heterogeneity = TRUE)
#' # Top 9
#' plot_dsge(result, n = 9)
#' # Selected GO terms
#' plot_dsge(result, go_ids = c("GO:0007156", "GO:0007268"))
#' # Standardised scale
#' plot_dsge(result, n = 9, use_std = TRUE)
#' }
plot_dsge <- function(result, n = 9L,
                       go_ids   = NULL,
                       col_null = "#2166AC",
                       col_tail = "#E41A1C",
                       col_obs  = "#D73027",
                       col_thr  = "#999999",
                       safety_margin = NULL,
                       cex_main = 0.85,
                       use_std  = NULL) {
  # ---- Input check ----
  if (!is.list(result) || !all(c("table", "null_raw") %in% names(result)))
    stop("'result' must be from pathway_dsge(..., return_null = TRUE)", call. = FALSE)

  # Auto-read safety_margin from result; user can still override
  if (is.null(safety_margin))
    safety_margin <- if (!is.null(result$safety_margin)) result$safety_margin else 1.05

  # ---- Auto-detect: if dsge_std column exists, default to standardised ----
  if (is.null(use_std))
    use_std <- "dsge_std" %in% names(result$table)
  if (isTRUE(use_std) && !"dsge_std" %in% names(result$table))
    stop("'use_std = TRUE' requires 'dsge_std' column; ",
         "run pathway_dsge() with use_std = TRUE", call. = FALSE)

  tbl <- result$table

  # ---- Select pathways to plot ----
  if (!is.null(go_ids)) {
    idx <- match(go_ids, tbl$go_id)
    if (all(is.na(idx)))
      stop("None of the specified 'go_ids' found in result", call. = FALSE)
    miss <- go_ids[is.na(idx)]
    if (length(miss) > 0)
      warning("GO id(s) not found: ", paste(miss, collapse = ", "), call. = FALSE)
    top <- tbl[na.omit(idx), ]
  } else {
    n   <- min(n, nrow(tbl))
    top <- tbl[seq_len(n), ]
  }
  n <- nrow(top)

  # ---- Layout: near-square grid ----
  cols <- ceiling(sqrt(n))
  rows <- ceiling(n / cols)
  old_par <- par(mfrow = c(rows, cols),
                 mar = c(3.2, 3.2, 2.8, 1),
                 mgp = c(1.8, 0.5, 0),
                 oma = c(0, 0, 2, 0),
                 bg  = "#FAFAFA")
  on.exit(par(old_par))

  for (i in seq_len(n)) {
    key     <- as.character(top$n_matched[i])
    null    <- result$null_raw[[key]]
    nm      <- top$go_name[i]
    go_id   <- top$go_id[i]
    p_val   <- top$p_value[i]
    p_adj   <- top$p_adj[i]
    asp     <- if ("aspect" %in% names(top)) top$aspect[i] else ""

    # ---- Standardise if requested ----
    if (isTRUE(use_std)) {
      null <- (null - mean(null)) / sd(null)
      obs  <- top$dsge_std[i]
    } else {
      obs  <- top$dsge[i]
    }

    # ---- Kernel density of null ----
    d_emp <- density(null)

    # ---- GPD fit for this size group ----
    gpd <- if (isTRUE(use_std)) {
      if ("null_gpd_std" %in% names(result)) result$null_gpd_std[[key]]
    } else {
      result$null_gpd[[key]]
    }

    # ---- Build combined density (empirical bulk + GPD tail) ----
    x_max_display <- max(d_emp$x, obs)
    x_max_display <- x_max_display + 0.08 * diff(range(c(d_emp$x, obs, null)))

    if (!is.null(gpd)) {
      # Safety margin: when shape < 0, adjust shape so GPD tail extends
      # to x_max_display, matching the logic in eval_gpd_p()
      shape_plot <- gpd$shape
      if (isTRUE(shape_plot < 0)) {
        theoretical_max <- gpd$u - gpd$scale / shape_plot
        if (isTRUE(x_max_display > theoretical_max)) {
          shape_plot <- -gpd$scale / ((x_max_display - gpd$u) * safety_margin)
        }
      }

      x_tail <- seq(gpd$u, x_max_display, length.out = 300)
      y_tail <- gpd$pat * evd::dgpd(x_tail - gpd$u,
                                     scale = gpd$scale, shape = shape_plot)
      keep   <- d_emp$x <= gpd$u
      d <- list(x = c(d_emp$x[keep], x_tail),
                y = c(d_emp$y[keep], y_tail),
                tail_idx = (length(which(keep)) + 1):(length(which(keep)) + length(x_tail)))
    } else {
      d <- list(x = d_emp$x, y = d_emp$y, tail_idx = integer(0))
    }

    # ---- Plot ----
    xlim <- range(c(d$x, obs, null))
    xlim[2] <- xlim[2] + 0.05 * diff(xlim)
    xlab <- if (isTRUE(use_std)) "DSGE (std)" else "DSGE"

    # Empty plot with light background
    plot(NA, xlim = xlim, ylim = c(0, max(d$y) * 1.08),
         main = "", xlab = xlab, ylab = "Density",
         las = 1, yaxt = "n",
         col.axis = "#333333", col.lab = "#333333",
         bty = "n")

    # Light panel grid
    grid_col <- "#E8E8E8"
    usr <- par("usr")
    abline(h = pretty(par("usr")[3:4]), col = grid_col, lwd = 0.4)
    abline(v = pretty(par("usr")[1:2]), col = grid_col, lwd = 0.4)

    # ---- Density fill + line ----
    # Fill under curve (semi-transparent null color)
    polygon(c(d$x, rev(d$x)), c(d$y, rep(0, length(d$y))),
            col = paste0(col_null, "20"), border = NA)
    # GPD tail portion filled with tail color
    if (length(d$tail_idx) > 0) {
      idx <- d$tail_idx
      polygon(c(d$x[idx], rev(d$x[idx])),
              c(d$y[idx], rep(0, length(idx))),
              col = paste0(col_tail, "25"), border = NA)
    }
    # Density line
    lines(d$x, d$y, col = col_null, lwd = 2)

    # ---- Titles ----
    if (nchar(nm) > 42) nm <- paste0(substr(nm, 1, 39), "...")
    name_line <- if (nzchar(asp)) paste0(nm, " (", asp, ")") else nm
    title(main = name_line, cex.main = cex_main, line = 1, font.main = 1,
          col.main = "#222222")
    title(main = go_id, cex.main = cex_main * 0.85, line = 0.2,
          font.main = 3, col.main = "#666666")

    # ---- GPD threshold line ----
    if (!is.null(gpd))
      abline(v = gpd$u, col = col_thr, lty = 3, lwd = 0.7)

    # ---- Observed DSGE line ----
    abline(v = obs, col = col_obs, lwd = 2.5, lty = 2)

    # ---- Annotation (always show observed value + p-value + FDR) ----
    p_text    <- if (p_val < 0.001) sprintf("p = %.1e", p_val) else sprintf("p = %.3f", p_val)
    padj_text <- if (p_adj < 0.001) sprintf("FDR = %.1e", p_adj) else sprintf("FDR = %.3f", p_adj)

    if (isTRUE(use_std)) {
      leg_labels <- c(sprintf("DSGE_std = %.3f", obs), p_text, padj_text)
    } else {
      leg_labels <- c(sprintf("DSGE = %.3f", obs), p_text, padj_text)
    }
    leg_cols <- c(col_obs, "#333333", "#333333")
    leg_lty  <- c(2, 0, 0)
    leg_lwd  <- c(2.5, NA, NA)

    # ---- Heterogeneity annotation ----
    has_het <- "gini" %in% names(top) && "null_gini_raw" %in% names(result)
    if (has_het) {
      gini_val <- top$gini[i]
      cv_val   <- top$cv[i]
      het_pval <- top$het_p_value[i]
      gini_text <- sprintf("Gini = %.3f", gini_val)
      cv_text   <- if (is.na(cv_val)) "CV = NA" else sprintf("CV = %.3f", cv_val)
      het_text  <- if (het_pval < 0.001) sprintf("het.p = %.1e", het_pval) else sprintf("het.p = %.3f", het_pval)
      dir_text <- if (!is.na(het_pval) && het_pval < 0.05) {
        p_left <- sum(result$null_gini_raw[[key]] <= gini_val) /
                  length(result$null_gini_raw[[key]])
        if (p_left < 0.5) " (het)" else " (uniform)"
      } else ""
      leg_labels <- c(leg_labels, gini_text, cv_text, paste0(het_text, dir_text))
      leg_cols   <- c(leg_cols, "#333333", "#333333", "#333333")
      leg_lty    <- c(leg_lty, 0, 0, 0)
      leg_lwd    <- c(leg_lwd, NA, NA, NA)
    }

    # ---- Legend in the top-left corner ----
    legend("topleft",
           legend = leg_labels,
           col    = leg_cols,
           lty    = leg_lty,
           lwd    = leg_lwd,
           pch    = rep(NA, length(leg_labels)),
           cex    = 0.6, bty = "n", inset = 0.02,
           text.col = "#333333")
  }

  # ---- Overall title ----
  if (!is.null(go_ids)) {
    title_text <- "Selected Pathways"
  } else {
    title_text <- sprintf("Top %d Pathways", n)
  }
  title(main = title_text, outer = TRUE, cex.main = 1.2, line = -0.2,
        font.main = 1, col.main = "#222222")
}


# =========================================================================
# Exported function 5: plot_dsge_volcano -- pathway-level volcano plot
# =========================================================================

#' Gene-level volcano plot for a specific pathway
#'
#' Plots a focused volcano plot showing only the genes belonging to a
#' specified GO pathway. Each point is one gene from the pathway, with
#' log2 fold change on the x-axis and statistical significance (-log10
#' p-value) on the y-axis. This allows visual inspection of the direction,
#' magnitude, and distribution of perturbation within a single pathway.
#'
#' @param de_results Data.frame of differential expression results. Must
#'   contain at least log2FC, p-value, and gene identifier columns.
#' @param dsge_result Optional result from \code{\link{pathway_dsge}()}.
#'   Pass either the full list (where \code{dsge_result$table} is used)
#'   or the result data.frame directly. When provided, the pathway's
#'   \code{dsge_std} and \code{p_adj} are shown in the annotation.
#'   Default \code{NULL}.
#' @param pathway_genes A named list of pathway-gene mappings, as used in
#'   \code{\link{pathway_dsge}()}. Each element is a \code{data.frame}
#'   with a gene identifier column.
#' @param go_id A single character string specifying which pathway (GO ID)
#'   to plot. Must be a name in \code{pathway_genes}.
#' @param logFC_col Column name in \code{de_results} for log2 fold change.
#'   Default \code{"log2FoldChange"}.
#' @param pval_col Column name in \code{de_results} for the p-value
#'   (raw or adjusted). Default \code{"pvalue"}.
#' @param gene_col Column name in \code{de_results} for gene identifiers
#'   that match the pathway mapping. Default \code{"geneName"}.
#' @param gene_id_col Column name in \code{pathway_genes[[go_id]]} that
#'   holds the gene identifiers. Default \code{"db_object_symbol"}.
#' @param threshold p-value significance threshold for the horizontal
#'   reference line. Default \code{0.05}.
#' @param lfc_threshold Numeric vector of logFC thresholds for vertical
#'   reference lines (e.g. \code{c(-1, 1)}). Default \code{NULL}.
#' @param padj_col Optional column name in \code{de_results} for adjusted
#'   p-values (e.g. \code{"padj"}). When provided, points are classified
#'   as adjusted-significant using \code{padj_threshold}, and an
#'   additional "Nominal" category appears in the legend. Default
#'   \code{NULL} (uses \code{threshold} only).
#' @param padj_threshold Adjusted p-value threshold for significance
#'   classification. Default \code{0.05}. Only used when
#'   \code{padj_col} is provided.
#' @param color_up Color for significantly up-regulated genes.
#'   Default \code{"#CC3333"}.
#' @param color_down Color for significantly down-regulated genes.
#'   Default \code{"#3366CC"}.
#' @param color_ns Color for non-significant genes.
#'   Default \code{"#AAAAAA"}.
#' @param color_nom Color for nominally significant genes (raw p-value
#'   below \code{threshold} but adjusted p-value above
#'   \code{padj_threshold}). Default \code{"#CC9933"}. Only used when
#'   \code{padj_col} is provided.
#' @param alpha_sig Transparency for significant points.
#'   Default \code{0.9}.
#' @param alpha_ns Transparency for non-significant points.
#'   Default \code{0.5}.
#' @param cex_point Point size for significant genes. Non-significant
#'   genes are plotted at \code{0.7 * cex_point}. Default \code{1.6}.
#' @param label Whether to label genes with their names. Auto-set to
#'   \code{TRUE} when the pathway has \eqn{\le 80} matched genes,
#'   \code{FALSE} otherwise. Can be forced with \code{TRUE} or
#'   \code{FALSE}.
#' @param label_genes Optional character vector of specific gene names
#'   within the pathway to label. Overrides \code{label}. Useful for
#'   very large pathways where you only want to highlight a few key genes.
#' @param label_sig When \code{TRUE}, only label genes that pass the
#'   \code{threshold}. Default \code{FALSE}. Ignored when
#'   \code{label_genes} is provided.
#' @param cex_label Text size for gene labels. Default \code{0.65}.
#' @param xlab,ylab Axis labels. Auto-generated when \code{NULL}.
#' @param main Plot title. Default \code{NULL}, auto-generated from
#'   GO ID and GO name.
#' @param go_name Optional GO term name for the title. When \code{NULL},
#'   the function looks for a \code{go_name} column in the DSGE result
#'   table or uses only the GO ID.
#' @param ... Additional arguments passed to \code{\link[graphics]{plot}()}.
#'
#' @return Invisibly returns the subset of \code{de_results} for the
#'   pathway genes (data.frame).
#' @importFrom graphics abline legend points text
#' @export
#'
#' @examples
#' \dontrun{
#' res <- read.csv("deseq2_results.csv")
#' gaf <- read_gaf("goa_human.gaf")
#' go  <- read_obo("go.obo")
#' pw  <- get_pathway_genes(gaf, go_names = go, min_size = 15)
#' dsge <- pathway_dsge(pw, res$pvalue, res$baseMean, res$geneName, seed = 42)
#'
#' # T cell receptor complex, with DSGE stats
#' plot_dsge_volcano(res, dsge, pw, go_id = "GO:0042101")
#'
#' # With effect size thresholds, only label significant genes
#' plot_dsge_volcano(res, dsge, pw, go_id = "GO:0032720",
#'   lfc_threshold = c(-1, 1), label_sig = TRUE)
#'
#' # Large pathway, only label specific genes
#' plot_dsge_volcano(res, pw, go_id = "GO:0005737",
#'   label_genes = c("TP53", "MYC", "EGFR"))
#' }
plot_dsge_volcano <- function(de_results,
                               dsge_result    = NULL,
                               pathway_genes,
                               go_id,
                               logFC_col      = "log2FoldChange",
                               pval_col       = "pvalue",
                               padj_col       = NULL,
                               gene_col       = "geneName",
                               gene_id_col    = "db_object_symbol",
                               threshold      = 0.05,
                               padj_threshold = 0.05,
                               lfc_threshold  = NULL,
                               color_up       = "#CC3333",
                               color_down     = "#3366CC",
                               color_nom      = "#CC9933",
                               color_ns       = "#AAAAAA",
                               alpha_sig      = 0.9,
                               alpha_ns       = 0.5,
                               label          = NULL,
                               label_genes    = NULL,
                               label_sig      = FALSE,
                               cex_label      = 0.70,
                               cex_point      = 1.6,
                               xlab           = NULL,
                               ylab           = NULL,
                               main           = NULL,
                               go_name        = NULL,
                               ...) {
  stopifnot(is.data.frame(de_results), nrow(de_results) > 0)
  stopifnot(is.list(pathway_genes), go_id %in% names(pathway_genes))

  needed <- c(logFC_col, pval_col, gene_col)
  missing <- needed[!needed %in% names(de_results)]
  if (length(missing) > 0)
    stop("Columns not found in de_results: ",
         paste(missing, collapse = ", "), call. = FALSE)

  # ---- Extract pathway gene set ----
  pw_df   <- pathway_genes[[go_id]]
  if (!gene_id_col %in% names(pw_df))
    stop("Column '", gene_id_col, "' not found in pathway_genes[['",
         go_id, "]]", call. = FALSE)

  pw_genes <- unique(as.character(pw_df[[gene_id_col]]))
  n_pw     <- length(pw_genes)
  if (n_pw == 0)
    stop("Pathway '", go_id, "' has no genes in the mapping", call. = FALSE)

  # ---- Match pathway genes in DE results ----
  in_pw     <- de_results[[gene_col]] %in% pw_genes
  n_matched <- sum(in_pw)
  if (n_matched == 0)
    stop("None of the ", n_pw, " genes in pathway '", go_id,
         "' were found in de_results", call. = FALSE)

  # ---- Subset to pathway genes only ----
  pw_data <- de_results[in_pw, , drop = FALSE]

  # ---- Core plot data ----
  x_vals <- pw_data[[logFC_col]]
  y_vals <- -log10(pw_data[[pval_col]])

  y_finite <- y_vals[is.finite(y_vals)]
  if (length(y_finite) > 0) {
    y_max <- stats::quantile(y_finite, 0.995, na.rm = TRUE)
    y_vals[!is.finite(y_vals)] <- y_max * 1.1
  }

  is_raw_sig <- pw_data[[pval_col]] <= threshold
  is_adj_sig <- is_raw_sig  # fallback when no padj_col
  if (!is.null(padj_col) && padj_col %in% names(pw_data))
    is_adj_sig <- is_raw_sig & pw_data[[padj_col]] <= padj_threshold
  is_nom_only <- is_raw_sig & !is_adj_sig

  is_up     <- x_vals > 0
  is_down   <- x_vals < 0

  # ---- Point colors (adj-sig-up, adj-sig-down, nom-only, ns) ----
  point_col <- rep(grDevices::adjustcolor(color_ns,  alpha.f = alpha_ns),
                   length(x_vals))
  point_col[is_adj_sig & is_up]   <- grDevices::adjustcolor(color_up,   alpha.f = alpha_sig)
  point_col[is_adj_sig & is_down] <- grDevices::adjustcolor(color_down, alpha.f = alpha_sig)
  point_col[is_nom_only & is_up]   <- grDevices::adjustcolor(color_up,   alpha.f = alpha_sig)
  point_col[is_nom_only & is_down] <- grDevices::adjustcolor(color_down, alpha.f = alpha_sig)

  # ---- Point styling ----
  point_pch <- rep(16, length(x_vals))
  point_pch[is_nom_only] <- 1   # open circle for nominal-only
  point_pch[!is_raw_sig] <- 1   # open circle for ns
  point_cex <- rep(cex_point, length(x_vals))
  point_cex[!is_raw_sig] <- cex_point * 0.7

  # ---- Axis labels ----
  if (is.null(xlab)) xlab <- expression(log[2] ~ "fold change")
  if (is.null(ylab)) ylab <- expression(-log[10](italic(p)))

  # ---- Title ----
  if (is.null(main)) {
    nm <- if (!is.null(go_name) && nchar(go_name) > 0) {
      go_name
    } else if ("go_name" %in% names(pw_data) && nchar(pw_data$go_name[1]) > 0) {
      pw_data$go_name[1]
    } else {
      go_id
    }
    main <- paste0(nm, "  (", n_matched, "/", n_pw, " genes)")
  }

  # ---- Plot (clean, no box) ----
  old_par <- graphics::par(bty = "o", las = 1)
  on.exit(graphics::par(old_par))

  # Determine x-axis limits with symmetric padding
  x_abs_max <- max(abs(x_vals), na.rm = TRUE) * 1.15
  x_abs_max <- max(x_abs_max, if (!is.null(lfc_threshold)) max(abs(lfc_threshold)) * 1.3 else 0)
  y_abs_max <- max(y_vals, na.rm = TRUE) * 1.10

  graphics::plot(x_vals, y_vals,
                 xlim = c(-x_abs_max, x_abs_max),
                 ylim = c(0, y_abs_max),
                 xlab = xlab, ylab = ylab, main = main,
                 col = point_col, cex = point_cex, pch = point_pch,
                 bty = "o", las = 1, ...)

  # ---- Reference lines ----
  graphics::abline(h = -log10(threshold),      col = "#333333", lty = 2, lwd = 0.6)
  if (!is.null(padj_col) && padj_col %in% names(pw_data))
    graphics::abline(h = -log10(padj_threshold), col = "#888888", lty = 3, lwd = 0.6)
  graphics::abline(v = 0, col = "#333333", lty = 3, lwd = 0.4)
  if (!is.null(lfc_threshold)) {
    for (v in lfc_threshold)
      graphics::abline(v = v, col = "#333333", lty = 3, lwd = 0.4)
  }

  # ---- Legend ----
  use_adj <- !is.null(padj_col) && padj_col %in% names(pw_data)
  if (use_adj) {
    n_adj_up   <- sum(is_adj_sig & is_up,   na.rm = TRUE)
    n_adj_down <- sum(is_adj_sig & is_down, na.rm = TRUE)
    n_nom_up   <- sum(is_nom_only & is_up,   na.rm = TRUE)
    n_nom_down <- sum(is_nom_only & is_down, na.rm = TRUE)
    n_ns       <- sum(!is_raw_sig, na.rm = TRUE)

    graphics::legend("topright",
           legend = c(
            sprintf("Sig up (FDR) (%d)", n_adj_up),
            sprintf("Sig down (FDR) (%d)", n_adj_down),
            sprintf("Nominal up (%d)", n_nom_up),
            sprintf("Nominal down (%d)", n_nom_down),
            sprintf("NS (%d)", n_ns)
           ),
           col = c(grDevices::adjustcolor(color_up,   alpha.f = alpha_sig),
                   grDevices::adjustcolor(color_down, alpha.f = alpha_sig),
                   grDevices::adjustcolor(color_up,   alpha.f = alpha_sig),
                   grDevices::adjustcolor(color_down, alpha.f = alpha_sig),
                   grDevices::adjustcolor(color_ns,   alpha.f = alpha_ns)),
           pch = c(16, 16, 1, 1, 1),
           pt.cex = c(cex_point, cex_point, cex_point, cex_point, cex_point * 0.7),
           cex = 0.65, bty = "n", title = "Regulation")
  } else {
    n_up_sig   <- sum(is_raw_sig & is_up,   na.rm = TRUE)
    n_down_sig <- sum(is_raw_sig & is_down, na.rm = TRUE)
    n_ns       <- sum(!is_raw_sig, na.rm = TRUE)

    graphics::legend("topright",
           legend = c(
            sprintf("Up (%d)", n_up_sig),
            sprintf("Down (%d)", n_down_sig),
            sprintf("NS (%d)", n_ns)
           ),
           col = c(grDevices::adjustcolor(color_up,   alpha.f = alpha_sig),
                   grDevices::adjustcolor(color_down, alpha.f = alpha_sig),
                   grDevices::adjustcolor(color_ns,   alpha.f = alpha_ns)),
           pch = c(16, 16, 1),
           pt.cex = c(cex_point, cex_point, cex_point * 0.7),
           cex = 0.65, bty = "n", title = "Regulation")
  }

  # ---- Pathway-level DSGE stats ----
  if (!is.null(dsge_result)) {
    tbl <- if (is.data.frame(dsge_result)) dsge_result else dsge_result$table
    row_idx <- which(tbl$go_id == go_id)
    if (length(row_idx) > 0) {
      r <- tbl[row_idx[1], ]
      dsge_std_val <- if ("dsge_std" %in% names(r)) r$dsge_std else NA
      p_adj_val    <- if ("p_adj" %in% names(r)) r$p_adj else NA

      stat_lines <- character()
      if (!is.null(dsge_std_val) && !is.na(dsge_std_val))
        stat_lines <- c(stat_lines, sprintf("DSGE\u209B\u2099\u2091 = %.2f", dsge_std_val))
      if (!is.null(p_adj_val) && !is.na(p_adj_val)) {
        if (p_adj_val < 0.001)
          stat_lines <- c(stat_lines, sprintf("p.adj = %.1e", p_adj_val))
        else
          stat_lines <- c(stat_lines, sprintf("p.adj = %.4f", p_adj_val))
      }
      # Mean logFC
      stat_lines <- c(stat_lines,
                      sprintf("Mean log\u2082FC = %+.3f", mean(x_vals, na.rm = TRUE)))

      if (length(stat_lines) > 0) {
        graphics::legend("topleft",
               legend = stat_lines,
               cex = 0.65, bty = "n", inset = c(0.02, 0.02),
               text.col = "#333333")
      }
    }
  }

  # ---- Gene labels ----
  # Default: small pathways label all genes, large pathways only significant
  if (is.null(label)) {
    do_label <- n_matched <= 80
    if (!do_label && isTRUE(label_sig)) do_label <- TRUE
  } else {
    do_label <- isTRUE(label)
  }

  if (!is.null(label_genes)) {
    label_idx <- which(pw_data[[gene_col]] %in% label_genes)
  } else if (do_label) {
    if (isTRUE(label_sig)) {
      label_idx <- which(is_raw_sig)
      # When >80 significant genes, label only the top 80 to avoid clutter
      if (length(label_idx) > 80) label_idx <- label_idx[1:80]
    } else {
      label_idx <- seq_len(n_matched)
    }
  } else {
    label_idx <- integer()
  }

  if (length(label_idx) > 0) {
    label_text <- as.character(pw_data[[gene_col]][label_idx])

    # Add white background for readability
    for (i in seq_along(label_idx)) {
      j <- label_idx[i]
      graphics::text(x_vals[j], y_vals[j],
                     labels = label_text[i],
                     cex = cex_label, pos = 3,
                     offset = 0.3,
                     col = "#FFFFFF", xpd = TRUE)
    }
    # Overdraw in dark color for legible text
    for (i in seq_along(label_idx)) {
      j <- label_idx[i]
      graphics::text(x_vals[j], y_vals[j],
                     labels = label_text[i],
                     cex = cex_label, pos = 3,
                     offset = 0.3,
                     col = "#222222", xpd = TRUE)
    }
  }

  invisible(pw_data)
}
