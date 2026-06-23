# =========================================================================
# read_obo: Parse Gene Ontology OBO format files
# =========================================================================
#
# OBO (Open Biomedical Ontologies) is the standard exchange format for
# the Gene Ontology. An OBO file contains multiple [Term] stanzas, each
# defining a GO term's basic information: id, name, namespace, etc.
#
# This function extracts only the three core fields, ignoring
# is_a/relationship/alt_id and other relationship data, because DSGE
# analysis only needs GO term names for result readability.
#
# OBO file structure example:
#   format-version: 1.4
#   ...
#   [Term]
#   id: GO:0000001
#   name: mitochondrion inheritance
#   namespace: biological_process
#   ...
#   [Term]
#   id: GO:0000002
#   name: mitochondrial genome maintenance
#   namespace: biological_process
#   ...
#   [Typedef]
#   id: part_of
#   ...                          <- typedef stanzas are skipped
#
# Dependency: none; uses only base R readLines and regex matching.

#' Read GO term names from an OBO file
#'
#' Parses a Gene Ontology OBO format file (version 1.2 / 1.4), extracting
#' the \code{id}, \code{name}, and \code{namespace} from each [Term]
#' stanza. [Typedef] stanzas are skipped.
#'
#' @param file Path to the OBO file (e.g. \code{"go-basic.obo"}).
#'
#' @return A \code{data.frame} with columns:
#'   \item{id}{GO identifier (e.g. \code{"GO:0005515"})}
#'   \item{name}{Human-readable term name (e.g. \code{"protein binding"})}
#'   \item{namespace}{Ontology classification: \code{"molecular_function"},
#'     \code{"biological_process"}, or \code{"cellular_component"}}
#' @export
#'
#' @examples
#' # Create a temporary OBO file for demonstration
#' obo_file <- tempfile(fileext = ".obo")
#' writeLines(c(
#'   "format-version: 1.2",
#'   "",
#'   "[Term]",
#'   "id: GO:0003674",
#'   "name: molecular_function",
#'   "namespace: molecular_function",
#'   "",
#'   "[Term]",
#'   "id: GO:0005575",
#'   "name: cellular_component",
#'   "namespace: cellular_component",
#'   "",
#'   "[Term]",
#'   "id: GO:0008150",
#'   "name: biological_process",
#'   "namespace: biological_process"
#' ), obo_file)
#' go_names <- read_obo(obo_file)
#' head(go_names)
#' \donttest{
#' go_names <- read_obo("go-basic.obo")
#' head(go_names)
#' }
read_obo <- function(file) {
  # ---- Input validation ----
  if (!file.exists(file))
    stop("File does not exist: ", file, call. = FALSE)

  # ---- Read all lines into memory ----
  # OBO files are typically in the tens of MB; reading entirely is feasible
  lines <- readLines(file, warn = FALSE)

  # ---- Locate all [Term] stanza start lines ----
  start <- which(lines == "[Term]")
  if (length(start) == 0)
    stop("No [Term] stanzas found in file", call. = FALSE)

  # ---- Determine end line for each [Term] stanza ----
  # Strategy: find all stanza boundary lines (starting with "["); each
  # [Term] stanza extends from its own start line to the line before the
  # next stanza (or end of file)
  bounds <- which(grepl("^\\[", lines))
  end <- integer(length(start))
  for (i in seq_along(start)) {
    pos <- match(start[i], bounds)       # position of current [Term] in boundary list
    end[i] <- if (pos < length(bounds)) bounds[pos + 1] - 1 else length(lines)
  }

  # ---- Extract id, name, namespace from each [Term] stanza ----
  ids  <- character(length(start))   # GO identifier
  nms  <- character(length(start))   # GO term name
  nsps <- character(length(start))   # namespace

  for (i in seq_along(start)) {
    block <- lines[start[i]:end[i]]   # all lines of the current [Term] stanza

    # Extract "id:" line, strip "id: " prefix to get the raw ID
    id_line <- grep("^id:", block, value = TRUE)[1]
    ids[i]  <- if (!is.na(id_line)) sub("^id:\\s*", "", id_line) else NA_character_

    # Extract "name:" line, strip "name: " prefix to get the raw name
    nm_line <- grep("^name:", block, value = TRUE)[1]
    nms[i]  <- if (!is.na(nm_line)) sub("^name:\\s*", "", nm_line) else NA_character_

    # Extract "namespace:" line, strip "namespace: " prefix to get the ontology category
    ns_line <- grep("^namespace:", block, value = TRUE)[1]
    nsps[i] <- if (!is.na(ns_line)) sub("^namespace:\\s*", "", ns_line) else NA_character_
  }

  # ---- Assemble data.frame, filter out invalid rows without id ----
  df <- data.frame(id = ids, name = nms, namespace = nsps,
                   stringsAsFactors = FALSE)
  df[!is.na(df$id), , drop = FALSE]
}
