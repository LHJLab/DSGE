# =========================================================================
# get_pathway_genes_kegg: Build KEGG pathway-gene mapping from OrgDb
# =========================================================================
#
# Alternative to get_pathway_genes_db() that extracts KEGG pathway
# annotations (from the PATH column of an OrgDb) instead of GO terms.
# Returns a named list in the same structure that pathway_dsge() expects.
#
# KEGG pathway IDs embed an organism code prefix (e.g. "hsa00010" for
# human, "mmu00010" for mouse). The organism is detected automatically
# from the OrgDb via AnnotationDbi::species() and a hardcoded mapping
# table covering 15 common model organisms.
#
# Pathway names are fetched online via KEGGREST::keggList().
# When KEGGREST is not available or the network is down, names are set
# to NA and the function continues with a warning.
#
# Typical usage:
#   library(org.Hs.eg.db)
#   pw <- get_pathway_genes_kegg(org.Hs.eg.db)

# ---- KEGG organism code lookup table ----
# Maps common scientific names to KEGG 3-letter organism codes.
# AnnotationDbi::species() returns the scientific name, which we use
# as the lookup key.
KEGG_ORG_CODES <- c(
  "Homo sapiens"             = "hsa",
  "Mus musculus"             = "mmu",
  "Rattus norvegicus"        = "rno",
  "Danio rerio"              = "dre",
  "Drosophila melanogaster"  = "dme",
  "Caenorhabditis elegans"   = "cel",
  "Saccharomyces cerevisiae" = "sce",
  "Arabidopsis thaliana"     = "ath",
  "Escherichia coli"         = "eco",
  "Sus scrofa"               = "ssc",
  "Bos taurus"               = "bta",
  "Canis lupus familiaris"   = "cfa",
  "Macaca mulatta"           = "mcc",
  "Gallus gallus"            = "gga",
  "Xenopus tropicalis"       = "xtr"
)

#' Build KEGG pathway-gene mapping from a Bioconductor OrgDb object
#'
#' Extracts KEGG pathway to gene mappings from an \code{OrgDb} annotation
#' package (e.g. \code{org.Hs.eg.db}) or an \code{AnnotationHub} record.
#' Returns a named list in the same format as
#' \code{\link{get_pathway_genes}()}, ready to pass to
#' \code{\link{pathway_dsge}()}.
#'
#' KEGG pathway IDs are stored in the \code{PATH} column of OrgDb objects.
#' IDs include an organism prefix (e.g. \code{"hsa00010"} for human,
#' \code{"mmu00010"} for mouse). The organism code is auto-detected from
#' the OrgDb species and a built-in lookup table of 15 common model
#' organisms.
#'
#' Pathway names are fetched online via \code{KEGGREST::keggList()}.
#' When the \code{KEGGREST} package is not available or network access
#' fails, names are set to \code{NA} with a warning.
#'
#' @param orgdb An \code{OrgDb} object, e.g. \code{org.Hs.eg.db} for
#'   human, \code{org.Mm.eg.db} for mouse, or a similar object retrieved
#'   from \code{AnnotationHub}.
#' @param keytype Key type used to query the \code{OrgDb}. Default
#'   \code{"ENTREZID"}. Must be one of the key types returned by
#'   \code{keytypes(orgdb)}.
#' @param gene_id_col Name for the gene ID column in the output
#'   data.frames. Default \code{"kegg_gene_id"}.
#' @param gene_symbol_col Name for the gene symbol column in the output
#'   data.frames. Default \code{"kegg_gene_symbol"}.
#' @param min_size Minimum gene count per pathway; pathways below this
#'   are discarded. Default \code{5}. Pass \code{NULL} to keep all.
#' @param attach_path_names Logical. Whether to fetch pathway names via
#'   \code{KEGGREST}. Default \code{TRUE}. When \code{KEGGREST} is not
#'   available or network access fails, falls back to \code{NA} names.
#'
#' @return A named list where each element is a \code{data.frame} with
#'   columns \code{kegg_name} (if attached), \code{organism_code},
#'   \code{gene_id_col}, and \code{gene_symbol_col}. The list names are
#'   full KEGG pathway IDs including the organism prefix (e.g.
#'   \code{"hsa00010"}). S3 class \code{"kegg_pathway"} is set on list
#'   elements for source auto-detection by \code{\link{pathway_dsge}()}.
#'
#' @note Requires the \code{AnnotationDbi} package. The \code{KEGGREST}
#'   package is required only when \code{attach_path_names = TRUE}.
#'   Install with:
#'   \code{BiocManager::install(c("AnnotationDbi", "KEGGREST"))}.
#'
#' @seealso \code{\link{get_pathway_genes_db}} for the GO-based equivalent.
#'   \code{\link{get_pathway_genes_reactome}} for the Reactome equivalent.
#' @importFrom methods is
#' @export
#'
#' @examples
#' \donttest{
#' library(org.Hs.eg.db)
#'
#' # All KEGG pathways, minimum 5 genes
#' pw <- get_pathway_genes_kegg(org.Hs.eg.db)
#'
#' # Without fetching pathway names (offline)
#' pw <- get_pathway_genes_kegg(org.Hs.eg.db, attach_path_names = FALSE)
#' }
get_pathway_genes_kegg <- function(orgdb,
                                    keytype           = "ENTREZID",
                                    gene_id_col       = "kegg_gene_id",
                                    gene_symbol_col   = "kegg_gene_symbol",
                                    min_size          = 5L,
                                    attach_path_names = TRUE) {

  # ---- Dependency checks ----
  if (!requireNamespace("AnnotationDbi", quietly = TRUE))
    stop("Package 'AnnotationDbi' is required. ",
         "Install with: BiocManager::install('AnnotationDbi')",
         call. = FALSE)

  # ---- Input validation ----
  if (!is(orgdb, "OrgDb"))
    stop("'orgdb' must be an OrgDb object (e.g. org.Hs.eg.db)", call. = FALSE)

  # ---- Extract organism code ----
  # AnnotationDbi::species() returns a plain character string
  sp <- AnnotationDbi::species(orgdb)
  org_code <- KEGG_ORG_CODES[[sp]]
  if (is.null(org_code))
    stop("Unrecognised organism '", sp, "' for KEGG pathway lookup. ",
         "Supported organisms: ", paste(names(KEGG_ORG_CODES), collapse = ", "),
         call. = FALSE)

  # ---- Get all KEGG pathway annotations from the OrgDb ----
  all_keys <- AnnotationDbi::keys(orgdb, keytype = keytype)

  gs_data <- tryCatch(
    suppressWarnings(
      AnnotationDbi::select(orgdb,
                             keys    = all_keys,
                             columns = c("PATH", "SYMBOL"),
                             keytype = keytype)
    ),
    error = function(e) {
      stop("Failed to query KEGG pathway annotations: ", e$message,
           call. = FALSE)
    }
  )

  # ---- Clean up ----
  gs_data <- gs_data[!is.na(gs_data$PATH), , drop = FALSE]
  gs_data <- gs_data[!is.na(gs_data$SYMBOL), , drop = FALSE]

  # Prepend organism code to PATH IDs (OrgDb stores e.g. "00010",
  # KEGGREST returns "hsa00010" — make them consistent)
  gs_data$PATH <- paste0(org_code, gs_data$PATH)

  if (nrow(gs_data) == 0L)
    stop("No KEGG pathway annotations found for the given OrgDb and keytype",
         call. = FALSE)

  # ---- Rename columns ----
  key_col <- keytype
  cols <- colnames(gs_data)
  cols[cols == key_col]   <- gene_id_col
  cols[cols == "SYMBOL"]  <- gene_symbol_col
  cols[cols == "PATH"]    <- "kegg_id"
  colnames(gs_data) <- cols

  keep_cols <- intersect(c(gene_id_col, gene_symbol_col, "kegg_id"),
                          colnames(gs_data))
  gs_data <- gs_data[, keep_cols, drop = FALSE]

  # ---- Deduplicate within each KEGG pathway ----
  gs_data <- gs_data[!duplicated(gs_data[, c("kegg_id", gene_symbol_col)]), ,
                     drop = FALSE]

  # ---- Split by KEGG pathway ID ----
  # IMPORTANT: keep full IDs (e.g. "hsa00010"), do NOT strip prefix
  result <- split(gs_data[, c(gene_id_col, gene_symbol_col), drop = FALSE],
                  gs_data$kegg_id)
  result <- result[order(names(result))]

  # ---- Attach organism code ----
  for (nm in names(result)) {
    result[[nm]]$organism_code <- org_code
  }

  # ---- Attach pathway names via KEGGREST ----
  if (isTRUE(attach_path_names)) {
    if (!requireNamespace("KEGGREST", quietly = TRUE)) {
      warning("Package 'KEGGREST' not available; pathway names set to NA",
              call. = FALSE)
    } else {
      raw <- tryCatch(
        KEGGREST::keggList("pathway", org_code),
        error = function(e) {
          warning("KEGG name lookup failed: ", e$message,
                  call. = FALSE)
          NULL
        }
      )

      if (!is.null(raw)) {
        # KEGG names: vector, names=full IDs (e.g. "path:hsa00010"), values=descriptions
        # Strip "path:" prefix to match internal IDs (e.g. "hsa00010")
        id_stripped <- sub("^path:", "", names(raw))

        # Build name lookup: stripped ID -> pathway description
        name_lookup <- stats::setNames(raw, id_stripped)

        for (nm in names(result)) {
          # Full ID e.g. "hsa00010"
          name_val <- if (nm %in% names(name_lookup)) {
            name_lookup[[nm]]
          } else {
            NA_character_
          }
          result[[nm]]$kegg_name <- name_val
        }
      }
    }
  }

  # ---- Attach S3 class for source detection ----
  result <- lapply(result, function(df) {
    # Move organism_code and kegg_name to the front
    front <- c("kegg_name", "organism_code")
    front <- front[front %in% names(df)]
    rest  <- setdiff(names(df), front)
    df    <- df[, c(front, rest), drop = FALSE]
    structure(df, class = unique(c("kegg_pathway", class(df))))
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
