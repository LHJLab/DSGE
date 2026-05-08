# =========================================================================
# DSGE Pipeline — End-to-end pathway-level transcriptional perturbation analysis
# =========================================================================
#
# This script demonstrates a complete DSGE workflow, from raw DESeq2
# results to ranked pathway significance table. Edit the parameters and
# file paths below to run on your own data.
#
# Workflow overview (9 steps):
#   0. Load the DSGE package
#   1. Read DESeq2 differential expression results (CSV)
#   2. Read GAF gene annotation file
#   3. Read OBO ontology file (for GO term names)
#   4. Build pathway-gene mapping table
#   5. Run pathway DSGE analysis (permutation test + GPD tail extrapolation + BH FDR)
#   6. Print results summary
#   7. Print top 20 significant pathways
#   8. Save results to CSV
#   9. Plot: top 9 pathways + selected GO terms (null distribution density)
#
# Input requirements:
#   - DESeq2 results CSV: must have pvalue, baseMean, geneName columns
#   - GAF file: standard GAF 2.2 format (e.g. goa_human.gaf)
#   - OBO file: standard OBO format (e.g. go.obo)
#
# Output: CSV file under results/ with DSGE and FDR for each pathway

# =========================================================================
# 0. Load package
# =========================================================================
library(DSGE)

# =========================================================================
# Configurable parameters
# =========================================================================

# NormZ weighting: biological replicate count per gene
#   NULL    = unweighted (does not affect FDR ranking in standard DESeq2)
#   scalar  = all genes use the same replicate count
#   vector  = per-gene replicate counts (must match DESeq2 row count)
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
N_PERM        <- 10000L

# Random seed for reproducibility
SEED          <- 666L

# GAF column used to match gene names in DESeq2 results
#   Common options:
#     "db_object_symbol" — gene symbol (e.g. CALN1), recommended
#     "db_object_id"     — UniProt ID (e.g. Q16553)
GENE_ID_COL   <- "db_object_symbol"

# =========================================================================
# 1. Read DESeq2 differential expression results
# =========================================================================
cat("Loading DESeq2 results...\n")
res <- read.csv("data_exp/CALN1_GABA_W6_DESeq2_addTPM.csv", stringsAsFactors = FALSE)
cat("  Total genes:", nrow(res), "\n")

# Filter: keep only protein-coding genes with non-empty gene name (".")
# geneType == "protein_coding" excludes lncRNA, pseudogenes, etc.
res <- subset(res, geneType == "protein_coding" & geneName != ".")
cat("  Protein-coding with geneName:", nrow(res), "\n")

# =========================================================================
# 2. Read GAF gene annotations
# =========================================================================
cat("\nLoading GAF annotations...\n")
# read_gaf uses data.table::fread for efficient reading of large files
# (~160 MB, ~2-3 seconds)
gaf <- read_gaf("data_exp/goa_human.gaf/goa_human.gaf")
cat("  Rows:", nrow(gaf), "\n")

# =========================================================================
# 3. Read GO term names
# =========================================================================
cat("\nLoading GO term names...\n")
# read_obo parses the OBO file, extracting id/name/namespace
go <- read_obo("data_exp/go.obo")
cat("  Terms:", nrow(go), "\n")

# =========================================================================
# 4. Build pathway-gene mapping table
# =========================================================================
cat("\nBuilding pathway-gene map...\n")
# min_size here filters out pathways with too few genes in the annotation,
# reducing downstream computation. pathway_dsge() applies a second
# min_size/max_size filter based on actual matched gene counts.
pw <- get_pathway_genes(gaf, go_names = go, min_size = MIN_SIZE)
cat("  Pathways (>= ", MIN_SIZE, " genes in annotation):", length(pw), "\n")

# =========================================================================
# 5. Compute pathway DSGE (core step)
# =========================================================================
cat("\nComputing pathway DSGE...\n")
t0 <- Sys.time()

# return_null = TRUE retains null distribution data for plot_dsge()
result <- pathway_dsge(
  pathway_genes = pw,            # output of get_pathway_genes()
  pvalue        = res$pvalue,    # DESeq2 p-value column
  base_mean     = res$baseMean,  # DESeq2 baseMean column (filter low-expression genes)
  gene_names    = res$geneName,  # DESeq2 gene name column (matched to GAF db_object_symbol)
  gene_id_col   = GENE_ID_COL,   # matching column name
  n_replicates  = N_REPLICATES,  # NormZ weighting
  min_size      = MIN_SIZE,      # minimum matched gene count
  max_size      = MAX_SIZE,      # maximum matched gene count
  n_perm        = N_PERM,        # permutation count
  seed          = SEED,          # random seed
  progress      = TRUE,          # show progress bar
  heterogeneity = FALSE,         # compute perturbation heterogeneity (Gini, CV, het_p)
  return_null   = TRUE           # retain null distribution data for plotting
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
    "  DSGE range          : ", paste(round(range(result_tbl$dsge), 4), collapse = " ~ "), "\n")
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
               "p_value", "p_adj")
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

# 9a. Top 9 pathways
pdf("results/plots/top9_pathways.pdf", width = 12, height = 10)
plot_dsge(result, n = 9)
dev.off()
cat("  Top 9 plot -> results/plots/top9_pathways.pdf\n")

# 9b. Selected GO terms (example: top 2 significant pathways)
top_go <- head(result_tbl$go_id, 2)
pdf("results/plots/selected_pathways.pdf", width = 8, height = 4)
plot_dsge(result, go_ids = top_go)
dev.off()
cat("  Selected GO plot -> results/plots/selected_pathways.pdf (",
    paste(top_go, collapse = ", "), ")\n")

cat("\nDone.\n")
