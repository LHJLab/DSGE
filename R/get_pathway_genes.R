# =========================================================================
# get_pathway_genes: Extract pathway-gene mapping table from GAF data
# =========================================================================
#
# Splits the flat GAF data returned by read_gaf() into a named list by
# GO term. Each list element is a data.frame containing the gene
# information for that pathway.
#
# Main features:
#   1. Split by GO term (split by go_id)
#   2. Filter by qualifier, evidence code, and ontology aspect
#   3. Optionally remove duplicate gene entries within a pathway
#   4. Optionally attach GO term names (from read_obo output)
#   5. Filter by pathway gene count (min_size)
#
# The returned data structure is the direct input to pathway_dsge().

#' Extract gene-pathway (GO term) associations from GAF data
#'
#' Splits the parsed GAF data.frame returned by
#' \code{\link{read_gaf}()} into a named list by GO term. Supports
#' filtering by qualifier, evidence code, ontology aspect, and minimum
#' pathway size. Optionally attaches GO term names from an OBO file.
#'
#' @param gaf_data A \code{data.frame} returned by \code{\link{read_gaf}()}.
#' @param genes Character vector of column names identifying genes.
#'   Default \code{c("db_object_id", "db_object_symbol")}.
#'   The first is typically a UniProt ID, the second a gene symbol.
#' @param unique Logical. Whether to remove duplicate gene entries within
#'   a GO term. Default \code{TRUE} (the same gene may be annotated to
#'   the same term via multiple evidence lines).
#' @param min_size Minimum gene count threshold; pathways with fewer
#'   genes are discarded. Default \code{5}. Set to 5 to ensure sufficient
#'   statistical power for the permutation test.
#' @param qualifier Qualifier filter, e.g. \code{"enables"},
#'   \code{"involved_in"}, \code{"located_in"}. \code{NULL} (default)
#'   returns all. A common combination is \code{c("enables",
#'   "involved_in")} (excludes NOT annotations and located_in
#'   annotations).
#' @param evidence Evidence code filter, e.g. \code{"IDA"} (direct
#'   experimental assay), \code{"IEA"} (electronic annotation).
#'   \code{NULL} (default) returns all. For high-confidence manual
#'   annotations only, use e.g. \code{c("IDA", "IPI", "IMP", "IGI",
#'   "IEP")}.
#' @param aspect Ontology aspect filter: \code{"F"} (Molecular Function),
#'   \code{"P"} (Biological Process), \code{"C"} (Cellular Component).
#'   \code{NULL} (default) returns all three.
#' @param go_names A \code{data.frame} returned by
#'   \code{\link{read_obo}()}. When provided, each pathway data.frame
#'   gains \code{go_name} and \code{go_namespace} columns (placed in
#'   the first two positions). \code{go_namespace} uses abbreviated
#'   forms: \code{"BP"}, \code{"MF"}, \code{"CC"}.
#'
#' @return A named \code{list} where each element name is a GO term ID
#'   (e.g. \code{"GO:0005515"}) and the value is a \code{data.frame} of
#'   associated genes. The list is sorted alphabetically by GO ID.
#' @export
#'
#' @examples
#' \dontrun{
#' gaf <- read_gaf("goa_human.gaf")
#' go  <- read_obo("go.obo")
#'
#' # Get all pathways (with GO names, at least 5 genes)
#' pathway_genes <- get_pathway_genes(gaf, go_names = go, min_size = 5)
#' pathway_genes[["GO:0005515"]][1:5, ]
#'
#' # Only experimentally validated biological process annotations
#' pathway_genes <- get_pathway_genes(gaf, evidence = "IDA", aspect = "P")
#' }
get_pathway_genes <- function(gaf_data,
                               genes     = c("db_object_id", "db_object_symbol"),
                               unique    = TRUE,
                               min_size  = 5L,
                               qualifier = NULL,
                               evidence  = NULL,
                               aspect    = NULL,
                               go_names  = NULL) {

  # ---- Input validation ----
  if (!is.data.frame(gaf_data))
    stop("'gaf_data' must be a data.frame returned by read_gaf()", call. = FALSE)

  # Check that required columns exist (genes columns + go_id column)
  required <- c(genes, "go_id")
  miss <- setdiff(required, names(gaf_data))
  if (length(miss) > 0)
    stop("Required column(s) not found: ", paste(miss, collapse = ", "),
         call. = FALSE)

  # ---- Filter by qualifier / evidence / aspect ----

  # qualifier filter: column 4, e.g. "enables", "NOT|enables", etc.
  # NOTE: the NOT qualifier indicates the gene product does NOT have
  #       the specified function (negative annotation)
  dat <- gaf_data
  if (!is.null(qualifier)) {
    if (!"qualifier" %in% names(dat))
      stop("Column 'qualifier' not found", call. = FALSE)
    dat <- dat[dat$qualifier %in% qualifier, , drop = FALSE]
  }

  # evidence_code filter: column 7
  # IEA (electronic annotation) has the widest coverage but lower
  # confidence; IDA/IMP etc. are experimentally validated, higher
  # quality but lower coverage
  if (!is.null(evidence)) {
    if (!"evidence_code" %in% names(dat))
      stop("Column 'evidence_code' not found", call. = FALSE)
    dat <- dat[dat$evidence_code %in% evidence, , drop = FALSE]
  }

  # aspect filter: column 9, F/P/C
  # F = Molecular Function, P = Biological Process, C = Cellular Component
  if (!is.null(aspect)) {
    aspect <- match.arg(aspect, c("F", "P", "C"), several.ok = TRUE)
    if (!"aspect" %in% names(dat))
      stop("Column 'aspect' not found", call. = FALSE)
    dat <- dat[dat$aspect %in% aspect, , drop = FALSE]
  }

  # ---- Split by GO term ----
  # split() groups by go_id factor; each group retains only the columns
  # specified by 'genes'. Result is a named list: names = GO IDs,
  # values = data.frames.
  result <- split(dat[, genes, drop = FALSE], dat$go_id)

  # ---- Deduplicate ----
  # The same gene may be annotated to the same GO term via multiple
  # evidence sources (e.g. two publications with different PMIDs).
  # When unique = TRUE, each gene appears at most once per term.
  if (isTRUE(unique))
    result <- lapply(result, function(df) df[!duplicated(df), , drop = FALSE])

  # Sort by GO ID alphabetically
  result <- result[order(names(result))]

  # ---- Attach GO term names ----
  # Look up id-to-name mapping from read_obo() output and add go_name
  # and go_namespace columns at the front of each data.frame for
  # downstream readability
  if (!is.null(go_names)) {
    if (!is.data.frame(go_names) || !all(c("id", "name") %in% names(go_names)))
      stop("'go_names' must be a data.frame with 'id' and 'name' columns",
           call. = FALSE)

    # Build id -> name / id -> namespace lookup tables
    lookup_name <- setNames(go_names$name, go_names$id)
    if ("namespace" %in% names(go_names)) {
      lookup_ns <- setNames(go_names$namespace, go_names$id)
    } else {
      lookup_ns <- NULL
    }

    for (go in names(result)) {
      result[[go]]$go_name <- if (go %in% names(lookup_name)) lookup_name[[go]] else NA_character_
      if (!is.null(lookup_ns)) {
        # Map namespace to short form: biological_process -> BP, etc.
        ns <- lookup_ns[[go]]
        if (!is.na(ns)) {
          ns <- switch(ns,
                       biological_process = "BP",
                       molecular_function = "MF",
                       cellular_component = "CC",
                       ns)
        }
        result[[go]]$go_namespace <- if (!is.na(ns)) ns else NA_character_
      }
      # Move go_name, go_namespace to the front
      nc <- names(result[[go]])
      front <- c("go_name", "go_namespace")
      front <- front[front %in% nc]
      result[[go]] <- result[[go]][, c(front, setdiff(nc, front)), drop = FALSE]
    }
  }

  # ---- Filter by minimum gene count ----
  # Pathways with too few genes lack statistical power for the
  # permutation test; 5 is a common minimum. pathway_dsge() further
  # applies its own min_size/max_size filter downstream.
  if (!is.null(min_size)) {
    min_size <- as.integer(min_size)
    if (is.na(min_size) || min_size < 1L)
      stop("'min_size' must be a positive integer", call. = FALSE)
    keep <- vapply(result, nrow, integer(1L)) >= min_size
    result <- result[keep]
  }

  result
}
