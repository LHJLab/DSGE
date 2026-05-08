# =========================================================================
# read_gaf: Efficiently read GO Annotation File (GAF 2.x)
# =========================================================================
#
# GAF (Gene Ontology Annotation File) is the standard annotation format
# defined by the GO Consortium, recording associations between gene
# products (proteins/genes) and GO terms. This package supports
# GAF 2.2 (17 tab-delimited standard columns).
#
# GAF file characteristics:
#   - Files can be very large (human annotations > 800k rows, ~160 MB)
#   - Lines starting with "!" are comment/metadata headers at the top
#   - Data rows have 17 columns separated by tabs
#   - Some fields may be empty
#
# Why data.table::fread instead of read.table?
#   - fread() is a C-level multi-threaded reader, 10-50x faster than base R
#   - For a 160 MB GAF file: fread ~2-3 s, read.table ~30-60 s
#
# Dependency: data.table (fread)

# =========================================================================
# GAF 2.2 standard column names
# =========================================================================
# 17 columns in order, with meaning and examples in comments
GAF_COLUMNS <- c(
  "db",                    #  1 Database source, e.g. UniProtKB
  "db_object_id",          #  2 Unique gene product identifier in database, e.g. Q16553
  "db_object_symbol",      #  3 Gene symbol or name, e.g. CALN1 (key column for matching DESeq2 geneName)
  "qualifier",             #  4 Qualifier: enables/involved_in/located_in or NOT|enables, etc.
  "go_id",                 #  5 GO identifier, e.g. GO:0005515 (protein binding)
  "db_reference",          #  6 Reference, e.g. PMID:33961781 or GO_REF:0000002
  "evidence_code",         #  7 Evidence code: IEA (electronic)/IPI (protein interaction)/IDA (direct assay)/IBA (biological inference), etc.
  "with_from",             #  8 with/from field, auxiliary evidence info
  "aspect",                #  9 Ontology aspect: F (Molecular Function)/P (Biological Process)/C (Cellular Component)
  "db_object_name",        # 10 Full name of the gene product
  "db_object_synonym",     # 11 Gene synonym list (pipe-delimited)
  "db_object_type",        # 12 Object type, usually protein
  "taxon",                 # 13 Taxonomic identifier, e.g. taxon:9606 (human)
  "date",                  # 14 Annotation date, format YYYYMMDD
  "assigned_by",           # 15 Annotation source database
  "annotation_extension",  # 16 Annotation extension info
  "gene_product_form_id"   # 17 Isoform or specific gene product identifier
)

# =========================================================================
# Exported function 1: get_gaf_header — extract GAF comment header
# =========================================================================

#' Extract header lines from a GAF file
#'
#' Lines starting with \code{!} in GAF files contain metadata (database
#' version, date, column definitions, etc.). This function extracts them
#' for inspecting file version and provenance.
#'
#' @param file Path to the GAF file.
#'
#' @return Character vector, one string per header line (with leading
#'   \code{!} removed).
#' @export
#'
#' @examples
#' \dontrun{
#' get_gaf_header("goa_human.gaf")
#' }
get_gaf_header <- function(file) {
  # Comment headers are typically within the first 500 lines; no need
  # to read the entire file
  header_lines <- readLines(file, n = 500)
  comment_idx <- grep("^!", header_lines)
  if (length(comment_idx) == 0) return(character(0))
  # Strip the leading "!" to return clean metadata text
  sub("^!", "", header_lines[comment_idx])
}

# =========================================================================
# Exported function 2: read_gaf — read and parse a GAF file into a data.frame
# =========================================================================

#' Read a Gene Ontology Annotation File (GAF)
#'
#' Efficiently reads a GAF 2.x format gene ontology annotation file using
#' \code{data.table::fread()}. Comment header lines starting with
#' \code{!} are automatically skipped. Returns a data.frame with all 17
#' standard GAF columns.
#'
#' @param file Path to the GAF file.
#' @param col_names Character vector of column names. Defaults to the
#'   standard 17 GAF 2.2 column names. Pass \code{NULL} to keep
#'   auto-detected column names (if the file has a header row).
#' @param ... Additional arguments passed to \code{data.table::fread}.
#'
#' @return A \code{data.frame} containing the GAF annotation data.
#' @importFrom data.table fread
#' @export
#'
#' @examples
#' \dontrun{
#' gaf <- read_gaf("goa_human.gaf")
#' head(gaf)
#' }
read_gaf <- function(file, col_names = GAF_COLUMNS, ...) {
  # ---- File existence check ----
  if (!file.exists(file)) {
    stop("File does not exist: ", file, call. = FALSE)
  }

  # ---- Count comment header lines ----
  # The GAF spec requires all comment lines (starting with "!") to precede
  # data lines, so scanning the first 500 lines suffices to determine the
  # number of lines to skip. Typical GAF files have 20-50 header lines.
  first_lines <- readLines(file, n = 500L)
  n_skip <- sum(grepl("^!", first_lines))
  rm(first_lines)  # free memory promptly

  # ---- Read data with data.table::fread ----
  # sep = "\t"     tab-delimited
  # quote = ""     disable quote parsing (some GAF fields contain quotes
  #                that can cause misalignment)
  # header = FALSE GAF has no column header row; names supplied via col_names
  # skip = n_skip  skip comment header lines
  # fill = TRUE    allow rows with fewer than 17 columns (fill with NA for
  #                robustness)
  # na.strings = "" treat empty strings as NA
  dt <- data.table::fread(
    file      = file,
    sep       = "\t",
    quote     = "",
    header    = FALSE,
    skip      = n_skip,
    col.names = col_names,
    fill      = TRUE,
    na.strings = "",
    ...
  )

  # Convert to plain data.frame (data.table inherits from data.frame but
  # some functions behave differently; explicit conversion ensures
  # compatibility)
  as.data.frame(dt)
}
