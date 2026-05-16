# =========================================================================
# DSGE Pipeline — End-to-end pathway-level transcriptional perturbation analysis
# =========================================================================
#
# This script demonstrates a complete DSGE workflow, from raw differential
# expression results to ranked pathway significance table. Edit the
# parameters and file paths below to run on your own data.
#
# Two annotation modes are supported:
#   (A) GAF file mode (default): uses GAF + OBO files
#   (B) OrgDb mode: uses a Bioconductor OrgDb object (e.g., org.Hs.eg.db)
#       Set USE_ORGDB = TRUE and provide ORGDB below.
#
# Workflow overview (9 steps):
#   0. Load the DSGE package
#   1. Read differential expression results (CSV)
#   2. Read gene annotations (GAF or OrgDb)
#   3. Read GO term names (OBO or GO.db)
#   4. Build pathway-gene mapping table
#   5. Run pathway DSGE analysis (permutation test + optional GPD tail extrapolation + BH FDR)
#   6. Print results summary
#   7. Print top 20 significant pathways
#   8. Save results to CSV
#   9. Plot: top 9 pathways + selected GO terms (null distribution density)
#
# Input requirements:
#   - Differential expression CSV: must have at least pvalue, baseMean,
#     geneName columns (column names can be adapted; baseMean can be any
#     expression-level measure such as DESeq2 baseMean or Seurat avg_log2FC)
#   - GAF file (mode A): standard GAF 2.2 format (e.g. goa_human.gaf)
#   - OBO file (mode A): standard OBO format (e.g. go.obo)
#   - OrgDb (mode B): Bioconductor OrgDb package (e.g. org.Hs.eg.db)
#
# Output: CSV file under results/ with DSGE and FDR for each pathway

# =========================================================================
# 0. Source local functions (debug/test mode — replaces library(DSGE))
# =========================================================================
# Adjust this path to point to your local DSGE R/ directory
PKG_R_DIR <- "E:/DSGE/R"
source(file.path(PKG_R_DIR, "read_gaf.R"))
source(file.path(PKG_R_DIR, "read_obo.R"))
source(file.path(PKG_R_DIR, "get_pathway_genes.R"))
source(file.path(PKG_R_DIR, "get_pathway_genes_db.R"))
source(file.path(PKG_R_DIR, "dsge.R"))

# Required packages (must be installed; all called via :: so no library() needed):
#   data.table, POT, evd
#   Additional for OrgDb mode: AnnotationDbi, optionally GO.db

# =========================================================================
# Configurable parameters
# =========================================================================

# NormZ weighting: [Deprecated] No longer used; kept for backwards compatibility only
N_REPLICATES  <- NULL

# Pathway size filter (matched gene count)
#   min_size: pathways with fewer genes are not tested
#   max_size: pathways with more genes are excluded (overpowered)
#   Recommended range: [5, 500]
MIN_SIZE      <- 15L
MAX_SIZE      <- 1000L

# Permutation test parameters
#   n_perm: number of permutations per size group; larger = finer p-value resolution
#   10000 is a common balance point; 100000 gives higher resolution
N_PERM        <- 1000000L

# GPD tail extrapolation:
#   TRUE  — fit GPD to the tail of the null distribution for extreme-value
#            p-values (higher resolution for extreme observations).
#            Support-constrained adjustment (arXiv:2602.22975) prevents p = 0.
#   FALSE — use empirical ECDF only (p-values always >= 1/n_perm)
USE_GPD       <- TRUE

# GPD fitting parameters (used when USE_GPD = TRUE):
#   gpd_threshold — quantile threshold for tail fitting (default 0.90).
#                   Lower values give more tail samples (less variance) but
#                   may introduce bias; higher values do the opposite.
#   gpd_method    — estimation method: "mle" (maximum likelihood, default)
#                   or "pwm" (probability-weighted moments, more robust
#                   when n_perm is small or MLE fails to converge)
GPD_THRESHOLD <- 0.90
GPD_METHOD    <- "mle"

# Parallel computing: number of CPU cores for null distribution generation
#   Default 1 (sequential). Use parallel::detectCores() to use all cores.
#   Linux/macOS only (fork-based). Ignored on Windows.
N_CORES       <- 1L

# Random seed for reproducibility
SEED          <- 666L

# GAF column used to match gene names in differential expression results
#   Common options:
#     "db_object_symbol" — gene symbol (e.g.更 CALN1), recommended
#     "db_object_id"     — UniProt ID (e.g. Q16553)
GENE_ID_COL   <- "db_object_symbol"

# Annotation source:
#   USE_ORGDB = FALSE (default) — use GAF + OBO files
#   USE_ORGDB = TRUE           — use a Bioconductor OrgDb object
USE_ORGDB     <- FALSE

# OrgDb object (only used when USE_ORGDB = TRUE)
#   For common model organisms:
#     library(org.Hs.eg.db)  -> org.Hs.eg.db  (human)
#     library(org.Mm.eg.db)  -> org.Mm.eg.db  (mouse)
#     library(org.Dr.eg.db)  -> org.Dr.eg.db  (zebrafish)
#   For other species with a reference genome, use AnnotationHub:
#     hub <- AnnotationHub::AnnotationHub()
#     query(hub, "Ovis aries")
#     ORGDB <- hub[["AH72269"]]
ORGDB         <- NULL

# =========================================================================
# 1. Read differential expression results
# =========================================================================
cat("Loading differential expression results...\n")
res <- read.csv("inst/data_exp/CALN1_GABA_W6_DESeq2_addTPM.csv", stringsAsFactors = FALSE)
cat("  Total genes:", nrow(res), "\n")

# Filter: keep only protein-coding genes with non-empty gene name (".")
# geneType == "protein_coding" excludes lncRNA, pseudogenes, etc.
res <- subset(res, geneType == "protein_coding" & geneName != ".")
cat("  Protein-coding with geneName:", nrow(res), "\n")

# =========================================================================
# 2. / 3. / 4. Read annotations + Build pathway-gene mapping
# =========================================================================
# Two modes: (A) GAF + OBO files, (B) Bioconductor OrgDb object
# =========================================================================

if (isTRUE(USE_ORGDB)) {

  # ---- Mode B: OrgDb object ----
  if (is.null(ORGDB))
    stop("USE_ORGDB = TRUE but ORGDB is NULL. Set ORGDB to an OrgDb object.",
         call. = FALSE)

  cat("\nBuilding pathway-gene map from OrgDb...\n")
  pw <- get_pathway_genes_db(
    orgdb           = ORGDB,
    keytype         = "ENTREZID",
    gene_id_col     = "db_object_id",
    gene_symbol_col = "db_object_symbol",
    min_size        = MIN_SIZE,
    aspect          = NULL,           # all aspects (BP, MF, CC)
    evidence        = NULL,           # all evidence codes
    attach_go_names = TRUE            # fetch GO names via GO.db
  )
  cat("  Pathways (>= ", MIN_SIZE, " genes):", length(pw), "\n")

} else {

  # ---- Mode A: GAF + OBO files ----

  cat("\nLoading GAF annotations...\n")
  gaf <- read_gaf("inst/data_exp/goa_human.gaf/goa_human.gaf")
  cat("  Rows:", nrow(gaf), "\n")

  go <- read_obo("inst/data_exp/go.obo")
  cat("  GO terms:", nrow(go), "\n")

  cat("\nBuilding pathway-gene map...\n")
  pw <- get_pathway_genes(gaf, go_names = go, min_size = MIN_SIZE)
  cat("  Pathways (>= ", MIN_SIZE, " genes in annotation):", length(pw), "\n")

}

# =========================================================================
# 5. Compute pathway DSGE (core step)
# =========================================================================
cat("\nComputing pathway DSGE...\n")
t0 <- Sys.time()

# return_null = TRUE retains null distribution data for plot_dsge()
result <- pathway_dsge(
  pathway_genes = pw,            # output of get_pathway_genes()
  pvalue        = res$pvalue,    # p-value column (e.g. DESeq2, Seurat)
  base_mean     = res$baseMean,  # mean expression column (e.g. DESeq2 baseMean)
  gene_names    = res$geneName,  # gene symbol column (matched to GAF db_object_symbol)
  gene_id_col   = GENE_ID_COL,   # matching column name
  n_replicates  = N_REPLICATES,  # NormZ weighting
  min_size      = MIN_SIZE,      # minimum matched gene count
  max_size      = MAX_SIZE,      # maximum matched gene count
  n_perm        = N_PERM,        # permutation count
  seed          = SEED,          # random seed
  progress      = TRUE,          # show progress bar
  heterogeneity = FALSE,         # compute perturbation heterogeneity (Gini, CV, het_p)
  use_std       = TRUE,           # standardise DSGE: (obs - null_mean) / null_sd
  use_gpd       = USE_GPD,        # GPD tail extrapolation for extreme-value p-values
  gpd_threshold = GPD_THRESHOLD,  # GPD tail quantile threshold
  gpd_method    = GPD_METHOD,     # GPD estimation method ("mle" or "pwm")
  n_cores       = N_CORES,        # parallel cores (1 = sequential, Linux/macOS)
  return_null   = TRUE            # retain null distribution data for plotting
)
result_tbl <- result$table       # extract the results table

t1 <- Sys.time()
cat("\n  Time:", round(difftime(t1, t0, units = "secs"), 1), "seconds\n")
cat("  Pathways tested:", nrow(result_tbl), "\n")

# =========================================================================
# 6. Results summary
# =========================================================================
cat("\n",
    paste(rep("=", 56), collapse = ""), "\n",
    "  RESULTS SUMMARY\n",
    paste(rep("=", 56), collapse = ""), "\n",
    "  FDR < 0.05          : ", sum(result_tbl$p_adj < 0.05), " pathways\n",
    "  FDR < 0.10          : ", sum(result_tbl$p_adj < 0.10), " pathways\n",
    "  FDR < 0.20          : ", sum(result_tbl$p_adj < 0.20), " pathways\n",
    "  DSGE range          : ", paste(round(range(result_tbl$dsge), 4), collapse = " ~ "), "\n",
    "  DSGE_std range       : ", paste(round(range(result_tbl$dsge_std), 4), collapse = " ~ "), "\n")
if ("gini" %in% names(result_tbl)) {
  cat("  Gini range          : ",
      paste(round(range(result_tbl$gini, na.rm = TRUE), 4), collapse = " ~ "), "\n",
      "  Sig. heterogeneity  : ",
      sum(result_tbl$het_p_adj < 0.05, na.rm = TRUE), " pathways (het_p_adj < 0.05)\n")
}
cat("  Median n_matched    : ", median(result_tbl$n_matched), "\n")
if ("aspect" %in% names(result_tbl)) {
  cat("  Aspect breakdown    : BP=", sum(result_tbl$aspect == "BP"),
      ", MF=", sum(result_tbl$aspect == "MF"),
      ", CC=", sum(result_tbl$aspect == "CC"), "\n")
}

# =========================================================================
# 7. Print top 20 most significant pathways
# =========================================================================
cat("\n  TOP 20 (by FDR):\n")
cat(paste(rep("-", 102), collapse = ""), "\n")
show_cols <- c("go_id", "go_name", "aspect", "n_pathway", "n_matched", "dsge",
               "dsge_std", "p_value", "p_adj")
if ("gini" %in% names(result_tbl)) {
  show_cols <- c(show_cols, "gini", "cv", "het_p_value", "het_p_adj")
}
print(result_tbl[1:min(20, nrow(result_tbl)), show_cols, drop = FALSE],
      row.names = FALSE)

# =========================================================================
# 8. Save results to CSV
# =========================================================================
dir.create("results", showWarnings = FALSE)
outfile <- "results/CALN1_GABA_W6_pathway_dsge.csv"
write.csv(result_tbl, outfile, row.names = FALSE)
cat("\nResults saved to:", outfile, "\n")

# =========================================================================
# 9. Plot — null distribution vs. observed DSGE
# =========================================================================
cat("\nPlotting null distributions...\n")
dir.create("results/plots", showWarnings = FALSE)

# 9a. Top 9 pathways — standardised (auto-detected from dsge_std column)
pdf("results/plots/top9_pathways_std.pdf", width = 12, height = 10)
plot_dsge(result, n = 9)
dev.off()
cat("  Top 9 (standardised) plot -> results/plots/top9_pathways_std.pdf\n")

# 9b. Selected GO terms — standardised (auto-detected)
top_go <- head(result_tbl$go_id, 2)
pdf("results/plots/selected_pathways_std.pdf", width = 8, height = 4)
plot_dsge(result, go_ids = top_go)
dev.off()
cat("  Selected GO (standardised) plot -> results/plots/selected_pathways_std.pdf (",
    paste(top_go, collapse = ", "), ")\n")

# 9c. Top 9 pathways — raw scale (explicit use_std = FALSE)
pdf("results/plots/top9_pathways_raw.pdf", width = 12, height = 10)
plot_dsge(result, n = 9, use_std = FALSE)
dev.off()
cat("  Top 9 (raw) plot -> results/plots/top9_pathways_raw.pdf\n")

# =========================================================================
# 10. (Optional) Comparison run with use_gpd = FALSE
# =========================================================================
# Uncomment the block below to run a second analysis with GPD disabled.
# This allows comparison of p-values between GPD tail extrapolation and
# pure empirical ECDF. Requires re-running the full permutation step.
# =========================================================================
# cat("\n\n=== Comparison: use_gpd = FALSE (pure ECDF) ===\n")
# t0 <- Sys.time()
# result_ecdf <- pathway_dsge(
#   pathway_genes = pw,
#   pvalue        = res$pvalue,
#   base_mean     = res$baseMean,
#   gene_names    = res$geneName,
#   gene_id_col   = GENE_ID_COL,
#   min_size      = MIN_SIZE,
#   max_size      = MAX_SIZE,
#   n_perm        = N_PERM,
#   seed          = SEED,
#   progress      = TRUE,
#   heterogeneity = FALSE,
#   use_std       = TRUE,
#   use_gpd       = FALSE,          # disable GPD
#   return_null   = FALSE
# )
# t1 <- Sys.time()
# cat("  Time:", round(difftime(t1, t0, units = "secs"), 1), "seconds\n")
# cat("  Pathways tested:", nrow(result_ecdf), "\n")
# cat("  FDR < 0.05:", sum(result_ecdf$p_adj < 0.05), "pathways\n")
# write.csv(result_ecdf, "results/CALN1_GABA_W6_pathway_dsge_ecdf.csv", row.names = FALSE)
# cat("  Results saved to: results/CALN1_GABA_W6_pathway_dsge_ecdf.csv\n")
# cat("  Compare: sum(p_adj < 0.05) with GPD =",
#     sum(result_tbl$p_adj < 0.05), "vs ECDF =", sum(result_ecdf$p_adj < 0.05), "\n")

cat("\nDone.\n")
