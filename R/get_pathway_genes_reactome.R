# =========================================================================
# get_pathway_genes_reactome: Build Reactome pathway-gene mapping
# =========================================================================
#
# Alternative to get_pathway_genes_db() that extracts Reactome pathway
# annotations using reactome.db (reactomePATHID2EXTID) instead of GO
# terms. Gene symbols are resolved from an OrgDb object.
#
# Reactome pathway IDs follow the pattern "R-HSA-<number>" (e.g.
# R-HSA-177929). Pathway names are fetched locally from the reactome.db
# Bioconductor package using batched queries (one select() call, not
# per-pathway).
#
# Typical usage:
#   library(org.Hs.eg.db)
#   pw <- get_pathway_genes_reactome(org.Hs.eg.db)

#' Build Reactome pathway-gene mapping from a Bioconductor OrgDb object
#'
#' Extracts Reactome pathway to gene mappings, using \code{reactome.db}
#' for pathway-to-gene associations and an \code{OrgDb} for gene symbol
#' resolution. Returns a named list ready to pass to
#' \code{\link{pathway_dsge}()}.
#'
#' Reactome pathway IDs follow the pattern \code{R-HSA-<number>}
#' (e.g. \code{R-HSA-177929}). Pathway names are fetched locally from
#' the \code{reactome.db} Bioconductor package in a single batched query.
#'
#' @param orgdb An \code{OrgDb} object (e.g. \code{org.Hs.eg.db}) used
#'   to resolve Entrez IDs to gene symbols. Only human (R-HSA) pathways
#'   are included.
#' @param keytype Key type used to query the \code{OrgDb}. Default
#'   \code{"ENTREZID"}. Must be one of the key types returned by
#'   \code{keytypes(orgdb)}.
#' @param gene_id_col Name for the gene ID column in the output
#'   data.frames. Default \code{"reactome_gene_id"}.
#' @param gene_symbol_col Name for the gene symbol column in the output
#'   data.frames. Default \code{"reactome_gene_symbol"}.
#' @param min_size Minimum gene count per pathway; pathways below this
#'   are discarded. Default \code{5}. Pass \code{NULL} to keep all.
#' @param attach_path_names Logical. Whether to fetch pathway names via
#'   \code{reactome.db}. Default \code{TRUE}. When \code{reactome.db} is
#'   not installed, falls back to \code{NA} names.
#' @param species_prefix Reactome species prefix for filtering pathways.
#'   Default \code{"R-HSA"} (Homo sapiens). Use \code{NULL} to keep all
#'   species.
#'
#' @return A named list where each element is a \code{data.frame} with
#'   columns \code{reactome_name} (if attached), \code{gene_id_col}, and
#'   \code{gene_symbol_col}. The list names are Reactome pathway IDs
#'   (e.g. \code{"R-HSA-177929"}). S3 class \code{"reactome_pathway"}
#'   is set on list elements for source auto-detection by
#'   \code{\link{pathway_dsge}()}.
#'
#' @note Requires \code{reactome.db} and \code{AnnotationDbi}. Install
#'   with: \code{BiocManager::install(c("AnnotationDbi", "reactome.db"))}.
#'
#' @seealso \code{\link{get_pathway_genes_db}} for the GO-based equivalent.
#'   \code{\link{get_pathway_genes_kegg}} for the KEGG equivalent.
#' @importFrom methods is
#' @export
#'
#' @examples
#' # Requires additional Bioconductor database packages
#' \dontrun{
#' library(org.Hs.eg.db)
#' pw <- get_pathway_genes_reactome(org.Hs.eg.db)
#' str(pw[1:2])
#' }
get_pathway_genes_reactome <- function(orgdb,
                                        keytype           = "ENTREZID",
                                        gene_id_col       = "reactome_gene_id",
                                        gene_symbol_col   = "reactome_gene_symbol",
                                        min_size          = 5L,
                                        attach_path_names = TRUE,
                                        species_prefix    = "R-HSA") {

  # ---- Dependency checks ----
  if (!requireNamespace("reactome.db", quietly = TRUE))
    stop("Package 'reactome.db' is required. ",
         "Install with: BiocManager::install('reactome.db')",
         call. = FALSE)
  if (!requireNamespace("AnnotationDbi", quietly = TRUE))
    stop("Package 'AnnotationDbi' is required. ",
         "Install with: BiocManager::install('AnnotationDbi')",
         call. = FALSE)

  # ---- Input validation ----
  if (!is(orgdb, "OrgDb"))
    stop("'orgdb' must be an OrgDb object (e.g. org.Hs.eg.db)", call. = FALSE)

  # ---- Get pathway-to-gene mappings from reactome.db ----
  # reactomePATHID2EXTID is an AnnDbBimap (S4), not an environment.
  # toTable() converts it to a data.frame (path_id, gene_id).
  gs_data <- suppressWarnings(
    AnnotationDbi::toTable(reactome.db::reactomePATHID2EXTID)
  )

  # Rename columns positionally: 1st = pathway ID, 2nd = Entrez ID
  colnames(gs_data)[1:2] <- c("reactome_id", "entrez_id")

  # Filter by species prefix on pathway ID
  if (!is.null(species_prefix)) {
    gs_data <- gs_data[grepl(paste0("^", species_prefix), gs_data$reactome_id), ,
                       drop = FALSE]
  }

  if (nrow(gs_data) == 0L)
    stop("No Reactome pathways found for prefix '", species_prefix, "'",
         call. = FALSE)

  # ---- Resolve gene symbols from OrgDb ----
  all_entrez <- unique(gs_data$entrez_id)

  sym_data <- tryCatch(
    suppressWarnings(
      AnnotationDbi::select(orgdb,
                             keys    = all_entrez,
                             columns = "SYMBOL",
                             keytype = keytype)
    ),
    error = function(e) {
      stop("Failed to query gene symbols from OrgDb: ", e$message,
           call. = FALSE)
    }
  )

  colnames(sym_data)[colnames(sym_data) == keytype] <- "entrez_id"

  # Merge symbols onto pathway data
  gs_data <- merge(gs_data, sym_data[, c("entrez_id", "SYMBOL")],
                   by = "entrez_id", all.x = TRUE)
  gs_data <- gs_data[!is.na(gs_data$SYMBOL), , drop = FALSE]

  if (nrow(gs_data) == 0L)
    stop("No Reactome pathways with resolvable gene symbols found",
         call. = FALSE)

  # ---- Rename columns ----
  colnames(gs_data)[colnames(gs_data) == "entrez_id"] <- gene_id_col
  colnames(gs_data)[colnames(gs_data) == "SYMBOL"]     <- gene_symbol_col

  keep_cols <- intersect(c(gene_id_col, gene_symbol_col, "reactome_id"),
                          colnames(gs_data))
  gs_data <- gs_data[, keep_cols, drop = FALSE]

  # ---- Deduplicate within each Reactome pathway ----
  gs_data <- gs_data[!duplicated(gs_data[, c("reactome_id", gene_symbol_col)]), ,
                     drop = FALSE]

  # ---- Split by Reactome pathway ID ----
  result <- split(gs_data[, c(gene_id_col, gene_symbol_col), drop = FALSE],
                  gs_data$reactome_id)
  result <- result[order(names(result))]

  # ---- Attach pathway names via reactome.db (batched) ----
  if (isTRUE(attach_path_names)) {
    if (!requireNamespace("reactome.db", quietly = TRUE)) {
      warning("Package 'reactome.db' not available; pathway names set to NA",
              call. = FALSE)
    } else {
      all_ids <- names(result)
      name_df <- tryCatch(
        suppressWarnings(
          AnnotationDbi::select(reactome.db::reactome.db,
                                 keys    = all_ids,
                                 keytype = "PATHID",
                                 columns = "PATHNAME")
        ),
        error = function(e) {
          warning("Reactome name lookup failed: ", e$message,
                  call. = FALSE)
          NULL
        }
      )

      if (!is.null(name_df) && nrow(name_df) > 0) {
        # Build lookup: PATHID -> PATHNAME
        # PATHNAME may contain "Homo sapiens: " prefix, remove it
        name_df$PATHNAME <- sub("^[^:]+:\\s*", "", name_df$PATHNAME)
        name_lookup <- stats::setNames(name_df$PATHNAME, name_df$PATHID)

        for (nm in names(result)) {
          name_val <- if (nm %in% names(name_lookup)) name_lookup[[nm]] else NA_character_
          result[[nm]]$reactome_name <- name_val
        }
      }
    }
  }

  # ---- Attach S3 class for source detection ----
  result <- lapply(result, function(df) {
    # Move reactome_name to the front
    if ("reactome_name" %in% names(df)) {
      df <- df[, c("reactome_name", setdiff(names(df), "reactome_name")),
                drop = FALSE]
    }
    structure(df, class = unique(c("reactome_pathway", class(df))))
  })

  # ---- Filter by minimum gene count ----
  if (!is.null(min_size)) {
    min_size <- as.integer(min_size)
    if (is.na(min_size) || min_size < 1L)
      stop("'min_size' must be a positive integer", call. = FALSE)
    keep <- vapply(result, nrow, integer(1L)) >= min_size
    result <- result[keep]
  }

  result
}
