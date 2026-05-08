# =========================================================================
# DSGE: Disruption Score of Gene Expression
# Pathway-level transcriptional perturbation analysis
# =========================================================================
#
# Core idea:
#   1. Convert DESeq2 differential expression p-values to absolute z-scores
#      (|z-score|); larger z means stronger transcriptional perturbation.
#   2. For each GO pathway (gene set), compute the NormZ-weighted mean of
#      member gene z-scores as the pathway DSGE.
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
#   Subset NormZ:        NormZ = Σ(N_i · z_i) / √(Σ N_i²)
#   Pathway DSGE:        DSGE = NormZ / G = Σ(N_i · z_i) / [G · √(Σ N_i²)]
#
#   where N_i is the number of biological replicates for gene i,
#   and G is the number of genes in the subset.
#   When all N_i are equal: DSGE = mean(z_i) / √G
#   Unweighted (n_replicates = NULL): DSGE = mean(z_i)
#
#   Note: The √(Σ N_i²) normalization is applied per subset, not
#   globally. This ensures reasonable scale relationships across
#   pathways of different sizes.
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
# NormZ weighting and normalization are applied later per subset.
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


# ---- Internal helper: compute DSGE on a gene subset (with NormZ weighting) ----
#
# Given a gene index subset, compute DSGE:
#   Unweighted:         DSGE = mean(z_i)
#   NormZ weighted:     DSGE = Σ(N_i · z_i) / [G · √(Σ N_i²)]
#
# When all N_i are equal, the weighted formula reduces to mean(z_i) / √G.
#
# Args:
#   idx    - indices of genes in pool_z
#   pool_z - raw z-scores of the full gene pool
#   pool_N - replicate count vector (NULL = unweighted)
#
# Returns: scalar DSGE
compute_dsge <- function(idx, pool_z, pool_N = NULL) {
  if (is.null(pool_N)) {
    # Unweighted: plain mean of z-scores
    mean(pool_z[idx])
  } else {
    # NormZ weighted: Σ(N·z) / [G · √(Σ N²)]
    N_sub <- pool_N[idx]
    z_sub <- pool_z[idx]
    sum(N_sub * z_sub) / (length(idx) * sqrt(sum(N_sub^2)))
  }
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
#   pool_N - replicate count vector (NULL = unweighted)
#
# Returns: numeric vector of length nb (one DSGE per column)
compute_dsge_batch <- function(mat, pool_z, pool_N = NULL) {
  # vector[matrix] in R drops dimensions; must restore explicitly
  val_z <- pool_z[mat]; dim(val_z) <- dim(mat)

  if (is.null(pool_N)) {
    # Unweighted: column means
    colMeans(val_z)
  } else {
    # NormZ weighted: Σ(N·z) / [sz · √(Σ N²)] per column
    val_N <- pool_N[mat]; dim(val_N) <- dim(mat)
    sz <- nrow(mat)
    colSums(val_N * val_z) / (sz * sqrt(colSums(val_N^2)))
  }
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
fit_gpd_tail <- function(null, tail = 0.90) {
  u <- unname(quantile(null, tail))

  # Too few samples above threshold (< 10): tail fit unreliable, skip GPD
  if (sum(null > u) < 10) return(NULL)

  # Fit GPD using POT::fitgpd with maximum likelihood estimation
  fit <- tryCatch(
    POT::fitgpd(null, u, est = "mle"),
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
# pgpd returns 0 (a genuine p=0).
#
# Args:
#   fit - GPD parameter list from fit_gpd_tail()
#   obs - observed DSGE value
#
# Returns: right-tail p-value (range [0, 1])
eval_gpd_p <- function(fit, obs) {
  # pat = P(X > u), pgpd(*, lower.tail = FALSE) = P(X > obs | X > u)
  p <- fit$pat * evd::pgpd(obs - fit$u, scale  = fit$scale,
                            shape = fit$shape, lower.tail = FALSE)
  # Clamp to [0, 1]; genuine p=0 is allowed
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
# Exported function 1: calc_dsge — genome-wide DSGE
# =========================================================================

#' Compute DSGE (Disruption Score of Gene Expression)
#'
#' Computes the mean absolute z-score for all genes passing filters
#' in a DESeq2 results object. Higher DSGE indicates stronger global
#' transcriptional perturbation.
#'
#' @param pvalue Numeric vector of p-values from DESeq2 results.
#' @param base_mean Numeric vector of baseMean from DESeq2 results
#'   (same length as pvalue). If NULL, no baseMean filtering is applied.
#' @param base_mean_cutoff baseMean filter threshold, default 0.1.
#'   Genes with baseMean at or below this value are excluded as lowly
#'   expressed.
#' @param n_replicates Number of biological replicates per gene. Can be
#'   a single value (all genes identical) or a vector matching pvalue
#'   length. NULL means no NormZ weighting. When all genes share the
#'   same replicate count, weighting reduces to constant scaling by
#'   1/√n_genes.
#'
#' @return A list with elements:
#'   \item{dsge}{scalar, genome-wide DSGE}
#'   \item{n_genes}{integer, number of genes passing filters}
#'   \item{z_scores}{named numeric vector of per-gene raw z-scores}
#' @importFrom stats qnorm quantile setNames
#' @export
#'
#' @examples
#' \dontrun{
#' res <- DESeq2::results(dds)
#' calc_dsge(res$pvalue, res$baseMean)
#' calc_dsge(res$pvalue, res$baseMean, n_replicates = 6)
#' }
calc_dsge <- function(pvalue, base_mean = NULL, base_mean_cutoff = 0.1,
                       n_replicates = NULL) {
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

  # ---- Process replicate counts ----
  if (!is.null(n_replicates)) {
    if (length(n_replicates) == 1L)
      n_replicates <- rep(n_replicates, sum(keep))
    else
      n_replicates <- n_replicates[keep]
  }

  # ---- Compute DSGE (NormZ formula over all retained genes) ----
  dsge_val <- compute_dsge(seq_along(pool_z), pool_z, n_replicates)

  list(dsge = dsge_val, n_genes = sum(keep), z_scores = pool_z)
}


# =========================================================================
# Exported function 2: dsge_perm_test — permutation test for a single gene set
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
#'         observed DSGE using the NormZ formula.
#'   \item Generate null distribution: randomly sample (without
#'         replacement) equally-sized gene sets from the pool,
#'         computing random DSGE via the same NormZ formula, repeated
#'         n_perm times.
#'   \item Empirical right-tail p-value: count(DSGE_null > DSGE_obs) / n_perm.
#' }
#'
#' @param gene_list Character vector of target gene identifiers
#'   (must match values in gene_names).
#' @param pvalue Numeric vector of DESeq2 p-values (all genes).
#' @param base_mean Numeric vector of DESeq2 baseMean values.
#' @param gene_names Character vector of gene names, same length as
#'   pvalue, must be unique.
#' @param base_mean_cutoff baseMean filter threshold, default 0.1.
#' @param n_replicates Number of biological replicates. NULL = unweighted
#'   (DSGE = mean(z)); non-NULL: DSGE = Σ(N·z) / [G · √(Σ N²)].
#' @param n_perm Number of permutations, default 10000.
#' @param seed Optional integer random seed for reproducibility.
#' @param progress Whether to show a progress bar, default TRUE.
#' @param heterogeneity Whether to compute perturbation heterogeneity
#'   indices (Gini, CV) and two-sided p-values. Default \code{FALSE}.
#'   When enabled, also computes Gini and CV null distributions within
#'   the permutation loop, increasing runtime by approximately
#'   30\%-50\%.
#'
#' @return A list with elements:
#'   \item{observed}{observed DSGE value}
#'   \item{n_genes}{number of target genes matched}
#'   \item{null}{permutation null distribution vector}
#'   \item{p_value}{empirical right-tail p-value}
#'   \item{ecdf}{empirical cumulative distribution function of the null}
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
#' }
dsge_perm_test <- function(gene_list, pvalue, base_mean, gene_names,
                           base_mean_cutoff = 0.1,
                           n_replicates     = NULL,
                           n_perm           = 10000L,
                           seed             = NULL,
                           progress         = TRUE,
                           heterogeneity    = FALSE) {
  # ---- Input validation ----
  stopifnot(is.character(gene_list), length(gene_list) > 0)
  stopifnot(length(pvalue) == length(base_mean))
  stopifnot(length(pvalue) == length(gene_names))
  if (anyDuplicated(gene_names))
    stop("'gene_names' must be unique", call. = FALSE)

  # ---- Build filtered gene pool ----
  keep <- !is.na(pvalue) & !is.na(base_mean) & base_mean > base_mean_cutoff
  pool_z <- compute_zscore(pvalue[keep])
  names(pool_z) <- gene_names[keep]
  n_pool <- length(pool_z)
  if (n_pool == 0)
    stop("No genes pass the baseMean > ", base_mean_cutoff, " filter", call. = FALSE)

  # ---- Set up replicate count vector ----
  if (!is.null(n_replicates)) {
    if (length(n_replicates) == 1L)
      pool_N <- rep(n_replicates, n_pool)
    else {
      if (length(n_replicates) != length(pvalue))
        stop("'n_replicates' must match 'pvalue' length", call. = FALSE)
      pool_N <- n_replicates[keep]
    }
  } else {
    pool_N <- NULL
  }

  # ---- Match target genes to gene pool ----
  idx <- match(gene_list, names(pool_z))
  idx <- idx[!is.na(idx)]
  if (length(idx) == 0)
    stop("No genes from 'gene_list' found in the filtered gene pool", call. = FALSE)

  n_target <- length(idx)
  observed <- compute_dsge(idx, pool_z, pool_N)

  if (heterogeneity) {
    observed_gini <- compute_gini(pool_z[idx])
    observed_cv   <- compute_cv(pool_z[idx])
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
    null[i] <- compute_dsge(samp, pool_z, pool_N)
    if (heterogeneity) {
      z_sub <- pool_z[samp]
      null_gini[i] <- compute_gini(z_sub)
      null_cv[i]   <- compute_cv(z_sub)
    }
    if (progress && i %% max(1L, floor(n_perm / 100)) == 0)
      utils::setTxtProgressBar(pb, i)
  }
  if (progress) { utils::setTxtProgressBar(pb, n_perm); close(pb) }

  # ---- Return results ----
  out <- list(observed = observed,
              n_genes  = n_target,
              null     = null,
              p_value  = sum(null > observed) / n_perm,
              ecdf     = stats::ecdf(null))

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
# Exported function 3: pathway_dsge — batch pathway DSGE with
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
#' **Algorithm steps (8 steps)**
#' \enumerate{
#'   \item \strong{Build gene pool}: filter DESeq2 results
#'         (baseMean > cutoff), compute per-gene raw z-scores.
#'   \item \strong{Match pathway genes}: match each pathway's genes to
#'         the gene pool via the \code{gene_id_col} column.
#'   \item \strong{Filter pathways}: retain pathways with
#'         \code{min_size <= n_matched <= max_size}.
#'   \item \strong{Compute observed DSGE}: per pathway using the NormZ
#'         formula: DSGE = Σ(N·z) / [G · √(Σ N²)].
#'   \item \strong{Generate size-grouped null distributions}: each unique
#'         matched-gene count gets its own permutation run (vectorized
#'         batch sampling) using the same NormZ formula, with
#'         simultaneous GPD tail fitting.
#'   \item \strong{Compute p-values}: GPD extrapolation when observed
#'         value falls in the tail (above 90th percentile); empirical
#'         ECDF for the body.
#'   \item \strong{BH FDR correction}: Benjamini-Hochberg multiple
#'         testing correction on all pathway p-values.
#'   \item \strong{Sort and return}: ordered by p_adj ascending.
#' }
#'
#' @param pathway_genes Named list returned by
#'   \code{\link{get_pathway_genes}()}. Each element is a data.frame
#'   containing gene information for one pathway.
#' @param pvalue Numeric vector of DESeq2 p-values (one per gene).
#' @param base_mean Numeric vector of DESeq2 baseMean values.
#' @param gene_names Character vector of gene names (same length as
#'   pvalue), must be unique.
#' @param gene_id_col Column name in pathway_genes data.frames used to
#'   match gene names, default \code{"db_object_symbol"}.
#' @param base_mean_cutoff baseMean filter threshold, default 0.1.
#' @param n_replicates Number of biological replicates. NULL = unweighted
#'   (DSGE = mean(z)); non-NULL: DSGE = Σ(N·z) / [G · √(Σ N²)].
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
pathway_dsge <- function(pathway_genes, pvalue, base_mean, gene_names,
                         gene_id_col       = "db_object_symbol",
                         base_mean_cutoff  = 0.1,
                         n_replicates      = NULL,
                         min_size          = 5L,
                         max_size          = 500L,
                         n_perm            = 10000L,
                         seed              = NULL,
                         return_null       = FALSE,
                         progress          = TRUE,
                         heterogeneity     = FALSE) {
  # ---- Input validation ----
  stopifnot(is.list(pathway_genes), length(pathway_genes) > 0)
  stopifnot(length(pvalue) == length(base_mean))
  stopifnot(length(pvalue) == length(gene_names))
  stopifnot(is.character(gene_id_col), length(gene_id_col) == 1L)
  stopifnot(min_size >= 1L, n_perm >= 1L)

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
    base_mean <- base_mean[uniq]
    gene_names <- gene_names[uniq]
    if (!is.null(n_replicates) && length(n_replicates) > 1L)
      n_replicates <- n_replicates[uniq]
  }

  # Filter: p-value non-missing, baseMean non-missing, baseMean > cutoff
  keep <- !is.na(pvalue) & !is.na(base_mean) & base_mean > base_mean_cutoff
  pool_z <- compute_zscore(pvalue[keep])        # raw z-scores
  names(pool_z) <- gene_names[keep]
  n_pool <- length(pool_z)
  if (n_pool == 0)
    stop("No genes pass the baseMean > ", base_mean_cutoff, " filter", call. = FALSE)

  # Set up replicate count vector (if provided)
  if (!is.null(n_replicates)) {
    if (length(n_replicates) == 1L)
      pool_N <- rep(n_replicates, n_pool)
    else
      pool_N <- n_replicates[keep]
    weighted <- TRUE
  } else {
    pool_N   <- NULL
    weighted <- FALSE
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
  # Step 4: Compute observed DSGE per pathway (and optional Gini, CV)
  # =========================================================================
  # Unweighted: DSGE = mean(z_i)
  # Weighted:   DSGE = Σ(N·z) / [G · √(Σ N²)]
  observed <- vapply(matched, function(idx) compute_dsge(idx, pool_z, pool_N),
                     numeric(1L))
  if (heterogeneity) {
    gini_obs <- vapply(matched, function(idx) compute_gini(pool_z[idx]),
                       numeric(1L))
    cv_obs   <- vapply(matched, function(idx) compute_cv(pool_z[idx]),
                       numeric(1L))
  }

  # =========================================================================
  # Step 5: Generate size-grouped null distributions + GPD tail fitting
  # =========================================================================
  # Pathways sharing the same matched-gene count reuse a single null.
  # Batch vectorized sampling: draw bat sets per iteration, compute
  # all at once with compute_dsge_batch.
  #
  sizes    <- unique(n_matched)
  null_raw <- vector("list", length(sizes))
  null_gpd <- vector("list", length(sizes))
  names(null_raw) <- names(null_gpd) <- as.character(sizes)
  if (heterogeneity) {
    null_gini_raw <- vector("list", length(sizes))
    null_cv_raw   <- vector("list", length(sizes))
    names(null_gini_raw) <- names(null_cv_raw) <- as.character(sizes)
  }

  if (!is.null(seed)) set.seed(seed)

  if (progress) {
    cat("Generating null distributions for", length(sizes),
        "unique pathway sizes\n")
    pb <- utils::txtProgressBar(min = 0, max = length(sizes), style = 3)
  }

  for (s in seq_along(sizes)) {
    sz  <- sizes[s]
    bat <- max(1L, floor(n_pool / sz))
    nul <- numeric(n_perm)
    if (heterogeneity) {
      nul_gini <- numeric(n_perm)
      nul_cv   <- numeric(n_perm)
    }

    for (b in seq(1L, n_perm, by = bat)) {
      nb <- min(bat, n_perm - b + 1L)
      mat <- matrix(sample.int(n_pool, sz * nb, replace = FALSE), nrow = sz)
      nul[b:(b + nb - 1L)] <- compute_dsge_batch(mat, pool_z, pool_N)
      if (heterogeneity) {
        nul_gini[b:(b + nb - 1L)] <- compute_gini_batch(mat, pool_z)
        nul_cv[b:(b + nb - 1L)]   <- compute_cv_batch(mat, pool_z)
      }
    }

    null_raw[[s]] <- nul
    null_gpd[[s]] <- fit_gpd_tail(nul)
    if (heterogeneity) {
      null_gini_raw[[s]] <- nul_gini
      null_cv_raw[[s]]   <- nul_cv
    }
    if (progress) utils::setTxtProgressBar(pb, s)
  }
  if (progress) { utils::setTxtProgressBar(pb, length(sizes)); close(pb) }

  # =========================================================================
  # Step 6: Compute p-value per pathway (GPD tail extrapolation + empirical ECDF)
  # =========================================================================
  # - Observed in tail (above 90th percentile u): GPD extrapolation
  #   p = P(X > u) * P(X > obs | X > u) = pat * pgpd(obs - u)
  # - Observed in body (≤ u): empirical ECDF direct count
  #
  p_val <- numeric(length(observed))
  for (i in seq_along(observed)) {
    key <- as.character(n_matched[i])
    gpd  <- null_gpd[[key]]

    if (!is.null(gpd) && observed[i] > gpd$u) {
      p_val[i] <- eval_gpd_p(gpd, observed[i])
    } else {
      p_val[i] <- sum(null_raw[[key]] > observed[i]) / n_perm
      if (p_val[i] == 0) p_val[i] <- 1 / n_perm
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
  # Step 7: Benjamini-Hochberg FDR correction
  # =========================================================================
  p_adj <- stats::p.adjust(p_val, method = "BH")
  if (heterogeneity)
    het_p_adj <- stats::p.adjust(het_p, method = "BH")

  # =========================================================================
  # Step 8: Assemble results, sorted by p_adj ascending
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
  out <- list(table = result, null_raw = null_raw, null_gpd = null_gpd)
  if (heterogeneity) {
    out$null_gini_raw <- null_gini_raw
    out$null_cv_raw   <- null_cv_raw
  }
  out
}


# =========================================================================
# Exported function 4: plot_dsge — null distribution density plots
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
#'
#' @return No return value; called for its side effect (plotting).
#' @importFrom graphics abline legend lines par rect title
#' @importFrom stats density na.omit
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
#' }
plot_dsge <- function(result, n = 9L,
                       go_ids   = NULL,
                       col_null = "steelblue",
                       col_obs  = "red",
                       col_tail = "#FFA50040",
                       cex_main = 0.85) {
  # ---- Input check ----
  if (!is.list(result) || !all(c("table", "null_raw") %in% names(result)))
    stop("'result' must be from pathway_dsge(..., return_null = TRUE)", call. = FALSE)

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
                 mar = c(3, 3, 2.5, 0.8),
                 mgp = c(1.7, 0.5, 0),
                 oma = c(0, 0, 1.5, 0))
  on.exit(par(old_par))

  for (i in seq_len(n)) {
    key     <- as.character(top$n_matched[i])
    null    <- result$null_raw[[key]]
    obs     <- top$dsge[i]
    nm      <- top$go_name[i]
    go_id   <- top$go_id[i]
    p_val   <- top$p_value[i]
    p_adj   <- top$p_adj[i]
    asp     <- if ("aspect" %in% names(top)) top$aspect[i] else ""

    # ---- Density curve ----
    d <- density(null)
    xlim <- range(c(d$x, obs, null))
    # Add a little space on the right to avoid line touching the edge
    xlim[2] <- xlim[2] + 0.05 * diff(xlim)
    plot(d, xlim = xlim, main = "",
         xlab = "DSGE", ylab = "Density",
         col = col_null, lwd = 2, las = 1, yaxt = "n")

    # ---- Title: GO name + ontology aspect + ID (truncate if too long) ----
    if (nchar(nm) > 40) nm <- paste0(substr(nm, 1, 37), "...")
    title_text <- if (nzchar(asp)) {
      bquote(.(nm) ~ "(" * .(asp) * ")" ~ "\n" ~ italic(.(go_id)))
    } else {
      bquote(.(nm) ~ "\n" ~ italic(.(go_id)))
    }
    title(main = title_text, cex.main = cex_main, line = -0.2)

    # ---- GPD tail region highlight ----
    gpd <- result$null_gpd[[key]]
    if (!is.null(gpd)) {
      tail_x <- seq(gpd$u, max(d$x, obs, gpd$u), length.out = 200)
      usr <- par("usr")
      y_bottom <- usr[3]
      # Semi-transparent rectangle marking the tail region
      rect(gpd$u, y_bottom, usr[2], usr[4],
           col = col_tail, border = NA)
      # Redraw density curve on top
      lines(d, col = col_null, lwd = 2)
      # Threshold vertical line (thin dashed)
      abline(v = gpd$u, col = "#888888", lty = 3, lwd = 0.8)
    }

    # ---- Observed DSGE vertical line ----
    abline(v = obs, col = col_obs, lwd = 2.5, lty = 2)

    # ---- p-value annotation (with heterogeneity info if available) ----
    p_text    <- if (p_val < 0.001) sprintf("p = %.1e", p_val) else sprintf("p = %.3f", p_val)
    padj_text <- if (p_adj < 0.001) sprintf("p.adj = %.1e", p_adj) else sprintf("p.adj = %.3f", p_adj)
    leg_labels <- c(sprintf("DSGE = %.3f", obs), p_text, padj_text)
    leg_cols   <- c(col_obs, "black", "black")
    leg_lty    <- c(2, 0, 0)
    leg_lwd    <- c(2.5, NA, NA)

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
        if (p_left < 0.5) " (heterogeneous)" else " (uniform)"
      } else ""
      leg_labels <- c(leg_labels, gini_text, cv_text, paste0(het_text, dir_text))
      leg_cols   <- c(leg_cols, "black", "black", "black")
      leg_lty    <- c(leg_lty, 0, 0, 0)
      leg_lwd    <- c(leg_lwd, NA, NA, NA)
    }
    legend("topleft",
           legend = leg_labels,
           col    = leg_cols,
           lty    = leg_lty,
           lwd    = leg_lwd,
           pch    = rep(NA, length(leg_labels)),
           cex    = 0.65, bty = "n", inset = 0.02)
  }

  # ---- Overall title ----
  if (!is.null(go_ids)) {
    title_text <- "Selected Pathways — Null Distribution vs. Observed DSGE"
  } else {
    title_text <- sprintf("Top %d Pathways — Null Distribution vs. Observed DSGE", n)
  }
  title(main = title_text, outer = TRUE, cex.main = 1.1, line = -0.5)
}
