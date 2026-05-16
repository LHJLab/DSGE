# DSGE: Disruption Score of Gene Expression

Pathway-level transcriptional perturbation analysis. Converts differential expression p-values into absolute z-scores and tests whether GO pathways show significantly stronger transcriptional disruption than expected by chance, using permutation-based null distributions with optional GPD tail extrapolation and Benjamini-Hochberg FDR correction.

## Installation

```r
# install.packages("devtools")   # if not already installed
devtools::install_github("LHJLab/DSGE")
```

## Quick Start

The package ships with an end-to-end pipeline script that you can copy and adapt. The walkthrough below uses the same example data step by step to demonstrate every function and its parameters.

```r
library(DSGE)
```

Our example data: a CALN1-GABA neuron DESeq2 comparison at week 6 (DESeq2 results used as example; any differential expression tool producing p-values works), a human GAF annotation file (~160 MB), and the GO OBO ontology.

### 1. Import differential expression results — `read.csv` (base R)

```r
res <- read.csv("data_exp/CALN1_GABA_W6_DESeq2_addTPM.csv", stringsAsFactors = FALSE)

# Filter to protein-coding genes with non-empty gene symbols
res <- subset(res, geneType == "protein_coding" & geneName != ".")
```

Required columns: `pvalue`, `baseMean`, `geneName` (column names can be adapted;
`baseMean` can be any expression-level measure; pass `base_mean = NULL` to
skip expression filtering).

### 2. Read GAF annotations — `read_gaf()`

Reads GAF 2.2 tab-separated files with `data.table::fread`. Comment lines starting with `!` are auto-skipped.

```r
gaf <- read_gaf("data_exp/goa_human.gaf/goa_human.gaf")

# Inspect the file metadata
head(get_gaf_header("data_exp/goa_human.gaf/goa_human.gaf"))
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `file` | _(required)_ | Path to the GAF file |
| `col_names` | `GAF_COLUMNS` | 17 GAF 2.2 column names; set `NULL` to auto-detect |
| `...` | — | Passed to `data.table::fread` |

### 3. Parse GO terms — `read_obo()`

Extracts `id`, `name`, and `namespace` from each `[Term]` stanza in OBO format.

```r
go <- read_obo("data_exp/go.obo")
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `file` | _(required)_ | Path to the OBO file |

### 4. Build pathway-gene map — `get_pathway_genes()`

Splits the GAF table into a named list by GO term. Each element is a data.frame of genes in that pathway.

```r
pw <- get_pathway_genes(
  gaf,
  genes     = c("db_object_id", "db_object_symbol"),   # which columns identify genes
  unique    = TRUE,                                     # deduplicate within each term
  min_size  = 5,                                        # drop pathways with < 5 genes
  qualifier = NULL,                                     # NULL = keep all; e.g. c("enables", "involved_in")
  evidence  = NULL,                                     # NULL = keep all; e.g. c("IDA", "IMP")
  aspect    = NULL,                                     # NULL = all; "P" = BP, "F" = MF, "C" = CC
  go_names  = go                                        # attach go_name + go_namespace (BP/MF/CC)
)
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `gaf_data` | _(required)_ | Output of `read_gaf()` |
| `genes` | `c("db_object_id", "db_object_symbol")` | Columns kept for downstream matching |
| `unique` | `TRUE` | Remove duplicate gene entries per term |
| `min_size` | `5` | Discard pathways below this gene count |
| `qualifier` | `NULL` | Filter by GAF qualifier column (e.g. `"enables"`) |
| `evidence` | `NULL` | Filter by evidence code (e.g. `c("IDA", "IPI")`) |
| `aspect` | `NULL` | Filter by ontology: `"P"`, `"F"`, `"C"` |
| `go_names` | `NULL` | Output of `read_obo()` — adds `go_name`, `go_namespace` columns |

### 5. (Alternative) Build pathway-gene map from OrgDb — `get_pathway_genes_db()`

An alternative to steps 2–4 for users who prefer Bioconductor's OrgDb
packages over GAF + OBO files. Produces the same output format.

#### Common model organisms

```r
library(org.Hs.eg.db)    # human
library(org.Mm.eg.db)    # mouse
library(org.Dr.eg.db)    # zebrafish
library(org.Rn.eg.db)    # rat
library(org.Dm.eg.db)    # fruit fly
library(org.Ce.eg.db)    # C. elegans
library(org.Sc.sgd.db)   # yeast
library(org.At.tair.db)  # Arabidopsis

pw <- get_pathway_genes_db(org.Hs.eg.db)
```

#### Non-model organisms via AnnotationHub

```r
library(AnnotationHub)
hub <- AnnotationHub()
query(hub, "Ovis aries")               # search for sheep
sheep_orgdb <- hub[["AH72269"]]        # load the OrgDb

pw <- get_pathway_genes_db(sheep_orgdb)
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `orgdb` | _(required)_ | An `OrgDb` object (e.g., `org.Hs.eg.db`) |
| `keytype` | `"ENTREZID"` | Key type for gene IDs in the OrgDb |
| `gene_id_col` | `"db_object_id"` | Gene ID column name in output |
| `gene_symbol_col` | `"db_object_symbol"` | Gene symbol column name in output |
| `min_size` | `5` | Drop pathways below this gene count |
| `aspect` | `NULL` | Ontology filter: `"BP"`, `"MF"`, `"CC"`, or `NULL` (all) |
| `evidence` | `NULL` | Evidence code filter (e.g., `"IDA"`); `NULL` = all |
| `attach_go_names` | `TRUE` | Fetch GO term names via `GO.db` |

### 6. Pathway DSGE analysis — `pathway_dsge()`

The core function. Computes DSGE for every pathway, generates size-grouped permutation null distributions, fits GPD to the upper tail, and applies Benjamini-Hochberg FDR correction.

```r
result <- pathway_dsge(
  pathway_genes    = pw,                      # from get_pathway_genes()
  pvalue           = res$pvalue,              # p-value column (any DE tool)
  base_mean        = res$baseMean,            # mean expression column
  gene_names       = res$geneName,            # gene symbols (match GAF gene_id_col)
  gene_id_col      = "db_object_symbol",      # column in pw to match gene_names against
  base_mean_cutoff = 0.1,                     # exclude genes with baseMean <= 0.1
  min_size         = 5,                       # drop pathways with fewer matched genes
  max_size         = 500,                     # drop pathways with more matched genes (Inf to keep all)
  n_perm           = 10000,                   # permutations per size group
  seed             = 42,                      # random seed for reproducibility
  return_null      = TRUE,                    # keep null distributions for plot_dsge()
  progress         = TRUE,                    # show progress bars
  heterogeneity    = FALSE,                   # compute Gini, CV, and het_p (adds ~30% runtime)
  use_std          = TRUE,                    # compute standardised DSGE (observed vs null)
  use_gpd          = TRUE,                     # GPD tail extrapolation for extreme-value p-values
  gpd_threshold    = 0.99,                     # GPD tail quantile threshold
  gpd_method       = "mle",                    # GPD estimation method
  n_cores          = 1                         # parallel cores (Linux/macOS only)
)

result_tbl <- result$table    # data.frame, sorted by p_adj ascending
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `pathway_genes` | _(required)_ | Named list from `get_pathway_genes()` |
| `pvalue` | _(required)_ | p-value vector from differential expression analysis (DESeq2, edgeR, Seurat, etc.) |
| `base_mean` | _(required)_ | Mean expression vector (e.g., DESeq2 baseMean, Seurat avg_log2FC); pass NULL to skip filtering |
| `gene_names` | _(required)_ | Gene symbols, must be unique |
| `gene_id_col` | `"db_object_symbol"` | Column in pathway data.frames to match `gene_names` |
| `base_mean_cutoff` | `0.1` | Exclude genes with baseMean at or below this value |
| `n_replicates` | `NULL` | **[Deprecated]** No longer used; will be removed in a future version |
| `min_size` | `5` | Minimum matched genes per pathway |
| `max_size` | `500` | Maximum matched genes (set `Inf` to disable) |
| `n_perm` | `10000` | Permutations per size group |
| `seed` | `NULL` | Random seed for reproducibility |
| `return_null` | `FALSE` | If `TRUE`, return list with null distributions (needed for `plot_dsge`) |
| `progress` | `TRUE` | Show progress bars during computation |
| `heterogeneity` | `FALSE` | If `TRUE`, also compute Gini, CV, and heterogeneity p-values |
| `use_std` | `TRUE` | If `TRUE`, compute `(observed - mean(null)) / sd(null)` and include `dsge_std` column |
| `use_gpd` | `TRUE` | If `TRUE`, use GPD tail extrapolation with support-constrained adjustment (avoids p=0). If `FALSE`, always use empirical ECDF (p-values always >= 1/n_perm) |
| `gpd_threshold` | `0.99` | Tail quantile threshold for GPD fitting. Lower = more tail samples (less variance, more bias); higher = fewer samples (more variance, less bias) |
| `gpd_method` | `"mle"` | GPD estimation method passed to `POT::fitgpd`. Default `"mle"`. Also: `"mple"`, `"moments"`, `"pwmu"`, `"pwmb"`, `"mdpd"`, `"med"`, `"pickands"`, `"lme"`, `"mgf"` |
| `n_cores` | `1` | Number of CPU cores for parallel null generation (Linux/macOS only, uses `parallel::mclapply`). Set to `parallel::detectCores()` to use all cores |

Result columns: `go_id`, `go_name`, `aspect`, `n_pathway`, `n_matched`, `dsge`, `dsge_std`, `p_value`, `p_adj`. When `heterogeneity = TRUE`: also `gini`, `cv`, `het_p_value`, `het_p_adj`.

### 7. Inspect and save results

```r
# How many pathways are significant?
sum(result_tbl$p_adj < 0.05)

# Ontology breakdown
table(result_tbl$aspect)    # BP / MF / CC

# Top pathways
head(result_tbl[, c("go_id", "go_name", "aspect", "n_matched", "dsge", "dsge_std", "p_adj")])

# Save
write.csv(result_tbl, "pathway_dsge_results.csv", row.names = FALSE)
```

### 8. Plot null distributions — `plot_dsge()`

Draws the null density curve for each selected pathway with the observed DSGE marked as a red dashed line. GPD tail region shown in orange. Set `use_std = TRUE` to plot on the standardised scale. Requires `pathway_dsge(..., return_null = TRUE)`.

```r
# Top 9 by significance
plot_dsge(result, n = 9)

# Specific GO terms
plot_dsge(result, go_ids = c("GO:0007156", "GO:0007268"))

# Standardised scale (requires dsge_std = TRUE in pathway_dsge)
plot_dsge(result, n = 9, use_std = TRUE)

# Save to PDF
pdf("top9.pdf", width = 12, height = 10)
plot_dsge(result, n = 9)
dev.off()
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `result` | _(required)_ | Output of `pathway_dsge(..., return_null = TRUE)` |
| `n` | `9` | Number of top pathways to plot (by p_adj); ignored when `go_ids` is set |
| `go_ids` | `NULL` | Character vector of GO terms to plot, e.g. `c("GO:0007156")` |
| `col_null` | `"steelblue"` | Null distribution density curve colour |
| `col_obs` | `"red"` | Observed DSGE vertical line colour |
| `col_tail` | `"#FFA50040"` | GPD tail region colour (semi-transparent) |
| `cex_main` | `0.85` | Title font scaling |
| `use_std` | auto | If `TRUE`, plot with standardised null and `dsge_std` on x-axis; defaults to `TRUE` when `dsge_std` column exists |

### 9. Test a single gene set — `dsge_perm_test()`

For testing one custom gene list without running the full pathway pipeline.

```r
my_genes <- c("CALN1", "GAD1", "GAD2", "SLC32A1", "SLC17A6")

test <- dsge_perm_test(
  gene_list        = my_genes,
  pvalue           = res$pvalue,
  base_mean        = res$baseMean,
  gene_names       = res$geneName,
  base_mean_cutoff = 0.1,
  n_perm           = 10000,
  seed             = 42,
  progress         = TRUE,
  heterogeneity    = FALSE,
  use_std          = TRUE,
  use_gpd          = TRUE
)
test$observed   # DSGE of the gene set
test$p_value    # empirical right-tail p-value
test$dsge_std   # standardised DSGE (if use_std = TRUE)
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `gene_list` | _(required)_ | Character vector of gene symbols |
| `pvalue` | _(required)_ | p-value vector from differential expression analysis |
| `base_mean` | _(required)_ | Mean expression vector; pass NULL to skip filtering |
| `gene_names` | _(required)_ | Gene symbols matching those in `gene_list` |
| `base_mean_cutoff` | `0.1` | Exclude genes with baseMean at or below this value |
| `n_perm` | `10000` | Number of permutations |
| `seed` | `NULL` | Random seed |
| `progress` | `TRUE` | Show progress bar |
| `heterogeneity` | `FALSE` | If `TRUE`, also compute Gini, CV, and het_p |
| `use_std` | `TRUE` | If `TRUE`, return `dsge_std = (observed - mean(null)) / sd(null)` |
| `use_gpd` | `TRUE` | If `TRUE`, use GPD tail extrapolation with support-constrained adjustment (avoids p=0). If `FALSE`, empirical ECDF only (p always >= 1/n_perm) |
| `gpd_threshold` | `0.99` | Tail quantile threshold for GPD fitting |
| `gpd_method` | `"mle"` | GPD estimation method passed to `POT::fitgpd`. Default `"mle"`. Also: `"mple"`, `"moments"`, `"pwmu"`, `"pwmb"`, `"mdpd"`, `"med"`, `"pickands"`, `"lme"`, `"mgf"` |

### 10. Genome-wide summary — `calc_dsge()`

A single scalar DSGE for the whole transcriptome, plus the per-gene z-score pool.

```r
dsge_res <- calc_dsge(
  pvalue           = res$pvalue,
  base_mean        = res$baseMean,
  base_mean_cutoff = 0.1
)
dsge_res$dsge       # scalar: global perturbation strength
dsge_res$n_genes    # number of genes after filtering
hist(dsge_res$z_scores, breaks = 100)   # per-gene z-score distribution
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `pvalue` | _(required)_ | p-value vector from differential expression analysis |
| `base_mean` | `NULL` | Mean expression vector (e.g., DESeq2 baseMean); NULL skips filtering |
| `base_mean_cutoff` | `0.1` | Exclude genes with mean expression at or below this value |

## Key Concepts

**DSGE formula.** `z_i = |Φ⁻¹(1 - p_i/2)|`; pathway DSGE = `mean(z_i)` (unweighted mean of absolute z-scores).

**Size-grouped permutation.** Pathways sharing the same matched-gene count reuse one null distribution, reducing computation from `K × n_perm` to `|sizes| × n_perm`.

**GPD tail extrapolation.** When `use_gpd = TRUE` (default), observed DSGE values above the `gpd_threshold` null percentile (default 0.99) get p-values from a fitted Generalized Pareto Distribution instead of direct counting — providing higher resolution for extreme observations. A **support-constrained adjustment** (Peschel et al. 2025, arXiv:2602.22975) is applied when the fitted GPD would otherwise produce p = 0 (due to a finite upper bound), ensuring a valid non-zero p-value while minimal deviation from the MLE. When set to `FALSE`, pure empirical ECDF is used (p-values always >= 1/n_perm, no risk of p = 0).

**Perturbation heterogeneity.** Optional Gini + CV with two-sided permutation test. Low heterogeneity = uniform pathway-wide perturbation; high heterogeneity = selective targeting of a few key genes (driver vs. passenger pathway).

## Input Data Sources

| File | Format | Where to get it |
|------|--------|-----------------|
| Differential expression results | CSV/table with `pvalue`, `baseMean`, `geneName` columns (column names can be adapted) | Your own DE analysis (DESeq2, edgeR, Seurat, limma, etc.) |
| GAF annotations (mode A) | GAF 2.2 (tab-separated, 17 cols) | [GOA Human](https://ftp.ebi.ac.uk/pub/databases/GO/goa/HUMAN/) |
| OBO ontology (mode A) | OBO 1.2/1.4 | [Gene Ontology Downloads](https://geneontology.org/docs/download-ontology/) |
| OrgDb (mode B) | Bioconductor OrgDb package | `BiocManager::install("org.Hs.eg.db")` or `AnnotationHub` |

## License

MIT
