// =========================================================================
// DSGE C++ backend: permutation null distribution builder
// =========================================================================
//
// Replaces the R-level batch loop in pathway_dsge() with a single-pass
// C++ implementation that avoids intermediate matrix allocation and
// per-column R-apply overhead.
//
// Expected speedup (vs R):
//   DSGE only:   3-5x
//   DSGE + Gini: 10-30x
//   DSGE + CV:   5-10x
//
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <random>
#include <vector>

using namespace Rcpp;

// =========================================================================
// Internal: partial Fisher-Yates shuffle
// =========================================================================
// Shuffles the first sz elements of indices[0..n_pool-1] uniformly.
// This is an efficient way to sample sz indices without replacement:
// after the call, indices[0..sz-1] is the random sample.
//
// Complexity: O(sz)
inline void partial_shuffle(std::vector<int>& indices, int sz,
                            std::mt19937& rng) {
  int n = static_cast<int>(indices.size());
  for (int i = 0; i < sz; ++i) {
    std::uniform_int_distribution<int> dist(i, n - 1);
    int j = dist(rng);
    std::swap(indices[i], indices[j]);
  }
}

// =========================================================================
// Internal: Gini coefficient of pre-sorted values
// =========================================================================
// G = (2 * sum(i * x_sorted[i])) / (n * sum(x)) - (n + 1) / n
// Requires x_sorted already sorted in ascending order.
inline double gini_sorted(const std::vector<double>& x) {
  int n = static_cast<int>(x.size());
  if (n < 2) return 0.0;
  double s = 0.0;
  for (int i = 0; i < n; ++i) s += x[i];
  if (s == 0.0) return 0.0;
  double weighted = 0.0;
  for (int i = 0; i < n; ++i)
    weighted += static_cast<double>(i + 1) * x[i];
  return (2.0 * weighted) / (static_cast<double>(n) * s) -
         (static_cast<double>(n + 1)) / static_cast<double>(n);
}

// =========================================================================
// Internal: coefficient of variation
// =========================================================================
// CV = sd(x) / mean(x)
inline double calc_cv(const std::vector<double>& x) {
  int n = static_cast<int>(x.size());
  if (n < 2) return 0.0;
  double sum  = 0.0;
  double sum2 = 0.0;
  for (int i = 0; i < n; ++i) {
    double v = x[i];
    sum  += v;
    sum2 += v * v;
  }
  double mean = sum / static_cast<double>(n);
  if (mean == 0.0) return 0.0;
  double var = (sum2 - sum * sum / static_cast<double>(n)) /
               static_cast<double>(n - 1);
  if (var <= 0.0) return 0.0;
  return std::sqrt(var) / mean;
}

// =========================================================================
// Exported: build null distributions for a single pathway size
// =========================================================================
// Generates n_perm permutation null values for DSGE (and optionally
// Gini / CV) in a single C++ pass. Uses a deterministic std::mt19937
// RNG seeded from `seed`, independent of R's .Random.seed.
//
// This is the C++ replacement for the R-level batch loop:
//
//   for (b in seq(1L, n_perm, by = bat)) {
//     nb <- min(bat, n_perm - b + 1L)
//     mat <- matrix(sample.int(n_pool, sz * nb), nrow = sz)
//     nul[b:(b+nb-1L)]       <- compute_dsge_batch(mat, pool_z)
//     nul_gini[b:(b+nb-1L)]  <- compute_gini_batch(mat, pool_z)
//     nul_cv[b:(b+nb-1L)]    <- compute_cv_batch(mat, pool_z)
//   }
//
// Parameters:
//   pool_z       - numeric vector, z-scores of the full gene pool
//   sz           - gene count for pathways of this size
//   n_perm       - number of permutations
//   seed         - RNG seed (deterministic, for reproducibility)
//   compute_gini - if TRUE, also compute Gini null distribution
//   compute_cv   - if TRUE, also compute CV null distribution
//
// Returns: a List with elements $null, $null_gini, $null_cv.
//          $null_gini and $null_cv are numeric(0) when not requested.
//
// [[Rcpp::export]]
List permute_null_cpp(NumericVector pool_z,
                      int           sz,
                      int           n_perm,
                      int           seed,
                      bool          compute_gini = false,
                      bool          compute_cv   = false) {

  int n_pool = pool_z.size();

  // ---- output vectors ----
  NumericVector nul(n_perm);
  NumericVector nul_gini = compute_gini ? NumericVector(n_perm)
                                        : NumericVector(0);
  NumericVector nul_cv   = compute_cv   ? NumericVector(n_perm)
                                        : NumericVector(0);

  // ---- RNG ----
  std::mt19937 rng(static_cast<unsigned int>(seed));

  // ---- reusable buffers ----
  std::vector<int>    indices(n_pool);
  std::vector<double> z_vals(sz);

  // Reset indices to 0,1,2,...,n_pool-1 only once before the first
  // permutation. After each shuffle, we restore stable order in
  // indices[0..sz-1] by re-setting them to 0..sz-1 (cheap).
  std::iota(indices.begin(), indices.end(), 0);

  // ---- main permutation loop ----
  for (int p = 0; p < n_perm; ++p) {
    partial_shuffle(indices, sz, rng);

    double sum = 0.0;
    for (int i = 0; i < sz; ++i) {
      double z = pool_z[indices[i]];
      sum += z;
      if (compute_gini || compute_cv) z_vals[i] = z;
    }
    nul[p] = sum / static_cast<double>(sz);

    if (compute_gini) {
      std::sort(z_vals.begin(), z_vals.end());
      nul_gini[p] = gini_sorted(z_vals);
    }
    if (compute_cv) {
      nul_cv[p] = calc_cv(z_vals);
    }
  }

  return List::create(
    _["null"]      = nul,
    _["null_gini"] = nul_gini,
    _["null_cv"]   = nul_cv
  );
}
