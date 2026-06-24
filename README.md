[![R-CMD-check](https://github.com/LHJLab/DSGE/actions/workflows/R-CMD-check.yaml/badge.svg?branch=v1.3.0)](https://github.com/LHJLab/DSGE/actions/workflows/R-CMD-check.yaml)
[![R package](https://github.com/LHJLab/DSGE/actions/workflows/R-package.yaml/badge.svg?branch=v1.3.0)](https://github.com/LHJLab/DSGE/actions/workflows/R-package.yaml)
[![lint](https://github.com/LHJLab/DSGE/actions/workflows/lint.yaml/badge.svg?branch=v1.3.0)](https://github.com/LHJLab/DSGE/actions/workflows/lint.yaml)

# DSGE: Disruption Score of Gene Expression

Gene set-level transcriptional perturbation analysis. Converts differential expression p-values into per-gene z-scores and tests whether each gene set shows stronger perturbation than expected under a size-matched permutation null. Uses GPD tail extrapolation for extreme-value p-values.

## Installation

```r
# install.packages("devtools")
devtools::install_github("LHJLab/DSGE")
```

## Quick Start

```r
library(DSGEr)
library(org.Hs.eg.db)

# Build pathway-gene map from Bioconductor OrgDb
pw <- get_pathway_genes_db(org.Hs.eg.db, min_size = 10)

# Read DE results (any tool — DESeq2, edgeR, limma, Seurat, etc.)
res <- read.csv("your_de_results.csv")
# Required columns: pvalue, gene symbol
# Optional: baseMean/AveExpr (for expression filtering), log2FoldChange (for direction)

# Run pathway analysis
result <- pathway_dsge(pw, pvalue = res$pvalue, base_mean = res$AveExpr,
                        gene_names = res$gene, gene_id_col = "db_object_symbol",
                        n_perm = 100000, n_cores = 4,
                        directional = TRUE, direction_vec = res$log2FoldChange,
                        return_null = TRUE)

# Significant pathways
head(result$table[result$table$p_adj < 0.05, c("go_id", "go_name", "dsge_std", "p_adj")])

# Plot null distribution for selected pathways
plot_dsge(result, go_ids = c("GO:0007264", "GO:0018108"))
```

Key parameters: `min_size`, `max_size`, `n_perm`, `use_gpd`, `directional` (with `direction_vec`), `nds_top_frac` , `n_cores`.

## Documentation

For detailed documentation with worked examples, parameter references, and diagnostic plots, see the [DSGE Wiki](https://github.com/LHJLab/DSGE/wiki).

## License

MIT
