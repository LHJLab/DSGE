# =========================================================================
# get_pathway_genes_db: Build pathway-gene mapping from an OrgDb object
# =========================================================================
#
# Alternative to get_pathway_genes() for users who have an OrgDb-style
# gene annotation database (from Bioconductor's org.*.eg.db packages or
# AnnotationHub). Produces the same named-list output that
# pathway_dsge() expects, without requiring GAF + OBO files.
#
# Key difference from get_pathway_genes():
#   - Input:  OrgDb object (programmatic access)
#   - Input:  GAF data.frame (file-based)
#
# Typical usage:
#   library(org.Hs.eg.db)
#   pw <- get_pathway_genes_db(org.Hs.eg.db)

#' Build pathway-gene mapping from a Bioconductor OrgDb object
#'
#' Extracts GO term to gene mappings from an \code{OrgDb} annotation
#' package (e.g. \code{org.Hs.eg.db}) or an \code{AnnotationHub} record.
#' Returns a named list in the same format as
#' \code{\link{get_pathway_genes}()}, ready to pass to
#' \code{\link{pathway_dsge}()}.
#'
#' This function provides an alternative to \code{get_pathway_genes()}
#' for users who prefer Bioconductor's annotation infrastructure over
#' GAF + OBO files.
#'
#' @param orgdb An \code{OrgDb} object, e.g. \code{org.Hs.eg.db} for
#'   human, \code{org.Mm.eg.db} for mouse, or a similar object retrieved
#'   from \code{AnnotationHub}.
#' @param keytype Key type used to query the \code{OrgDb}. Default
#'   \code{"ENTREZID"}. Must be one of the key types returned by
#'   \code{keytypes(orgdb)}.
#' @param gene_id_col Name for the gene ID column in the output
#'   data.frames. Default \code{"db_object_id"} (matching the
#'   \code{get_pathway_genes()} convention; the values are Entrez IDs
#'   by default).
#' @param gene_symbol_col Name for the gene symbol column in the output
#'   data.frames. Default \code{"db_object_symbol"} (matching the
#'   \code{get_pathway_genes()} convention).
#' @param min_size Minimum gene count per pathway; pathways below this
#'   are discarded. Default \code{5}. Pass \code{NULL} to keep all.
#' @param aspect Ontology aspect filter. \code{NULL} (default) returns
#'   all. One or more of \code{"BP"} (Biological Process),
#'   \code{"MF"} (Molecular Function), \code{"CC"} (Cellular Component).
#' @param evidence Evidence code filter (e.g. \code{"IDA"},
#'   \code{"IEA"}). \code{NULL} (default) keeps all. Pass a character
#'   vector to keep only annotations with those evidence codes.
#' @param use_goall Logical. If \code{TRUE}, query the \code{GOALL}
#'   column instead of \code{GO}, propagating gene annotations to all
#'   ancestor terms along the GO directed acyclic graph. This produces
#'   a broader pathway set consistent with \code{clusterProfiler}'s
#'   \code{gseGO}/\code{enrichGO} default behaviour. Default
#'   \code{FALSE} uses direct annotations only.
#' @param attach_go_names Logical. Whether to fetch GO term names and
#'   ontology classifications via \code{GO.db}. Default \code{TRUE}.
#'   Requires the \code{GO.db} package to be installed.
#'
#' @return A named list where each element is a \code{data.frame} with
#'   columns \code{go_name} (if attached), \code{go_namespace}
#'   (if attached), \code{gene_id_col}, and \code{gene_symbol_col}.
#'   The list names are GO term IDs. This matches the output format of
#'   \code{\link{get_pathway_genes}()} and can be passed directly to
#'   \code{\link{pathway_dsge}()}.
#'
#' @note This function requires the \code{AnnotationDbi} package. The
#'   \code{GO.db} package is required when \code{attach_go_names = TRUE}
#'   (the default). Both are Bioconductor packages; install with
#'   \code{BiocManager::install(c("AnnotationDbi", "GO.db"))}.
#'
#' @seealso \code{\link{get_pathway_genes}} for the GAF-based equivalent.
#' @importFrom methods is
#' @export
#'
#' @examples
#' # Requires additional Bioconductor database packages
#' \donttest{
#' if (require("org.Hs.eg.db", quietly = TRUE)) {
#'   pw <- get_pathway_genes_db(org.Hs.eg.db)
#'   str(pw[1:2])
#' }
#' }
get_pathway_genes_db <- function(orgdb,
                                  keytype          = "ENTREZID",
                                  gene_id_col      = "db_object_id",
                                  gene_symbol_col  = "db_object_symbol",
                                  min_size         = 5L,
                                  aspect           = NULL,
                                  evidence         = NULL,
                                  attach_go_names  = TRUE,
                                  use_goall        = FALSE) {

  # ---- Dependency checks ----
  if (!requireNamespace("AnnotationDbi", quietly = TRUE))
    stop("Package 'AnnotationDbi' is required. ",
         "Install with: BiocManager::install('AnnotationDbi')",
         call. = FALSE)

  # ---- Input validation ----
  if (!is(orgdb, "OrgDb"))
    stop("'orgdb' must be an OrgDb object (e.g. org.Hs.eg.db)", call. = FALSE)

  # ---- Get all GO annotations from the OrgDb ----
  # When use_goall=FALSE: select() with "GO" returns keytype, SYMBOL,
  #   GO, EVIDENCE, ONTOLOGY.
  # When use_goall=TRUE:  select() with "GOALL" returns keytype,
  #   SYMBOL, GOALL, EVIDENCEALL, ONTOLOGYALL — gene annotations are
  #   propagated to all ancestor terms along the GO DAG, producing a
  #   broader pathway set consistent with clusterProfiler's default.
  go_col  <- if (use_goall) "GOALL"       else "GO"
  evi_col <- if (use_goall) "EVIDENCEALL" else "EVIDENCE"
  ont_col <- if (use_goall) "ONTOLOGYALL" else "ONTOLOGY"

  all_keys <- AnnotationDbi::keys(orgdb, keytype = keytype)

  go_data <- tryCatch(
    suppressWarnings(
      AnnotationDbi::select(orgdb,
                             keys    = all_keys,
                             columns = c(go_col, "SYMBOL", evi_col),
                             keytype = keytype)
    ),
    error = function(e) {
      go_data <- suppressWarnings(
        AnnotationDbi::select(orgdb,
                               keys    = all_keys,
                               columns = c(go_col, "SYMBOL"),
                               keytype = keytype)
      )
      go_data[[evi_col]] <- NA_character_
      go_data[[ont_col]] <- NA_character_
      go_data
    }
  )

  # ---- Clean up ----
  go_data <- go_data[!is.na(go_data[[go_col]]), , drop = FALSE]
  go_data <- go_data[!is.na(go_data$SYMBOL), , drop = FALSE]

  if (nrow(go_data) == 0L)
    stop("No GO annotations found for the given OrgDb and keytype",
         call. = FALSE)

  # ---- Filter by evidence code ----
  if (!is.null(evidence)) {
    go_data <- go_data[go_data[[evi_col]] %in% evidence, , drop = FALSE]
    if (nrow(go_data) == 0L)
      stop("No rows remain after evidence code filter", call. = FALSE)
  }

  # ---- Filter by ontology aspect ----
  if (!is.null(aspect)) {
    aspect <- match.arg(aspect, c("BP", "MF", "CC"), several.ok = TRUE)
    go_data <- go_data[!is.na(go_data[[ont_col]]) &
                         go_data[[ont_col]] %in% aspect, , drop = FALSE]
    if (nrow(go_data) == 0L)
      stop("No rows remain after aspect filter", call. = FALSE)
  }

  # ---- Rename columns to match get_pathway_genes() convention ----
  key_col <- keytype
  cols <- colnames(go_data)
  cols[cols == key_col]       <- gene_id_col
  cols[cols == "SYMBOL"]      <- gene_symbol_col
  cols[cols == go_col]        <- "go_id"
  colnames(go_data) <- cols

  keep_cols <- intersect(c(gene_id_col, gene_symbol_col, "go_id"),
                          colnames(go_data))
  go_data <- go_data[, keep_cols, drop = FALSE]

  # ---- Deduplicate within each GO term ----
  # Use gene_symbol_col for dedup (more robust than gene_id since
  # different IDs may map to the same symbol)
  go_data <- go_data[!duplicated(go_data[, c("go_id", gene_symbol_col)]), ,
                     drop = FALSE]

  # ---- Split by GO term ----
  result <- split(go_data[, c(gene_id_col, gene_symbol_col), drop = FALSE],
                  go_data$go_id)
  result <- result[order(names(result))]

  # ---- Attach GO term names and ontology from GO.db ----
  if (isTRUE(attach_go_names)) {
    if (!requireNamespace("GO.db", quietly = TRUE))
      stop("Package 'GO.db' is required when attach_go_names = TRUE. ",
           "Install with: BiocManager::install('GO.db')",
           call. = FALSE)

    go_ids <- names(result)
    go_info <- tryCatch(
      suppressWarnings(
        AnnotationDbi::select(GO.db::GO.db,
                               keys    = go_ids,
                               columns = c("TERM", "ONTOLOGY"))
      ),
      error = function(e) data.frame(GOID = go_ids,
                                      TERM = NA_character_,
                                      ONTOLOGY = NA_character_,
                                      stringsAsFactors = FALSE)
    )

    lookup_name <- stats::setNames(go_info$TERM, go_info$GOID)
    lookup_ns   <- stats::setNames(go_info$ONTOLOGY, go_info$GOID)

    for (go in names(result)) {
      nm <- if (go %in% names(lookup_name)) lookup_name[[go]] else NA_character_
      ns <- if (go %in% names(lookup_ns)) lookup_ns[[go]] else NA_character_
      result[[go]]$go_name      <- nm
      result[[go]]$go_namespace <- ns
      # Move go_name and go_namespace to the front
      result[[go]] <- result[[go]][, c("go_name", "go_namespace",
                                        setdiff(names(result[[go]]),
                                                c("go_name", "go_namespace"))),
                                    drop = FALSE]
    }
  }

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
