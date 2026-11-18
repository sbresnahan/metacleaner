# Internal helpers -----------------------------------------------------------

#' Taxonomic ranks used by metacleaner
#'
#' The seven ranks that metacleaner reports and that can be used for
#' `level` and `require`. `"superkingdom"` and `"kingdom"` are accepted as
#' synonyms of `"domain"` (NCBI renamed this rank in 2025).
#'
#' @return A character vector of rank names, from broadest to narrowest.
#' @examples
#' taxonomic_ranks()
#' @export
taxonomic_ranks <- function() {
  c("domain", "phylum", "class", "order", "family", "genus", "species")
}

normalise_rank <- function(x, arg = "level") {
  if (!is.character(x) || length(x) != 1L || is.na(x)) {
    stop(sprintf("`%s` must be a single rank name, one of: %s.", arg,
                 paste(taxonomic_ranks(), collapse = ", ")), call. = FALSE)
  }
  y <- tolower(trimws(x))
  if (y %in% c("superkingdom", "kingdom", "realm")) y <- "domain"
  if (!y %in% taxonomic_ranks()) {
    stop(sprintf("`%s` = \"%s\" is not a recognized rank. Use one of: %s.",
                 arg, x, paste(taxonomic_ranks(), collapse = ", ")), call. = FALSE)
  }
  y
}

msg <- function(verbose, ...) {
  if (isTRUE(verbose)) message(...)
  invisible(NULL)
}

check_number <- function(x, arg, lower = -Inf, upper = Inf) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x < lower || x > upper) {
    stop(sprintf("`%s` must be a single number between %s and %s.", arg, lower, upper),
         call. = FALSE)
  }
  invisible(x)
}

check_flag <- function(x, arg) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(sprintf("`%s` must be TRUE or FALSE.", arg), call. = FALSE)
  }
  invisible(x)
}

check_threads <- function(threads) {
  if (identical(threads, "auto")) {
    n <- suppressWarnings(parallel::detectCores(logical = TRUE))
    return(max(1L, as.integer(n) - 1L, na.rm = TRUE))
  }
  check_number(threads, "threads", 1, 1024)
  as.integer(threads)
}

# Accession without version suffix, used for matching IDs between sources
base_accession <- function(x) sub("\\.[0-9]+$", "", x)

# Match ids against a key vector, first exactly and then ignoring versions
match_ids <- function(x, table) {
  m <- match(x, table)
  miss <- is.na(m)
  if (any(miss)) {
    m[miss] <- match(base_accession(x[miss]), base_accession(table))
  }
  m
}

# Turn a file path or a named character vector into a named character vector
as_sequences <- function(x, arg) {
  if (is.null(x)) return(NULL)
  if (is.character(x) && length(x) == 1L && is.null(names(x))) {
    if (!file.exists(x)) {
      stop(sprintf("The file given for `%s` does not exist: %s", arg, x), call. = FALSE)
    }
    x <- read_fasta(x)
  }
  if (!is.character(x) || is.null(names(x))) {
    stop(sprintf(paste0("`%s` must be the path to a FASTA file or a named character ",
                        "vector of sequences (for example the output of read_fasta())."), arg),
         call. = FALSE)
  }
  if (length(x) == 0L) {
    stop(sprintf("`%s` contains no sequences.", arg), call. = FALSE)
  }
  dup <- duplicated(names(x))
  if (any(dup)) {
    warning(sprintf(paste0("%d sequence(s) in `%s` have duplicated IDs; only the first ",
                           "sequence with each ID is used."), sum(dup), arg), call. = FALSE)
    x <- x[!dup]
  }
  x
}

collapse_lineage <- function(lin) {
  if (nrow(lin) == 0L) return(character(0))
  m <- as.matrix(lin[, taxonomic_ranks(), drop = FALSE])
  m[is.na(m)] <- "NA"
  gsub(" ", "_", apply(m, 1L, paste, collapse = "|"), fixed = TRUE)
}
