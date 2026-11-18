#' Search query sequences against a database of sequences
#'
#' Finds, for every query sequence, the database sequences that it matches
#' with at least `identity` percent identity over at least `coverage` percent
#' of the query length (both strands are searched). This is the search step
#' used inside [clean_reference()]; most users do not need to call it
#' directly.
#'
#' @section Native engine:
#' The default engine is written in C++ and needs no external software. The
#' database is indexed by k-mers (words of length `k`). For each query a
#' q-gram filter computes how many k-mers any qualifying match must share
#' with the query; only database sequences that pass this filter are
#' aligned (banded local alignment with megablast-like scores: match +1,
#' mismatch -2, gap 2.5). Because the filter is derived from the thresholds,
#' no qualifying match is missed when the thresholds are high (identity of
#' about 90% or more). With 100% identity and 100% coverage (the default) the
#' search is an exact containment test and is very fast.
#'
#' Identity is the number of identical positions divided by the alignment
#' length (as BLAST's `pident`). Upper/lower case and U/T are treated as the
#' same; an ambiguity code (e.g. `N`, `R`, `Y`) is identical only to the same
#' code. BLAST instead scores ambiguity codes as mismatches, so a sequence
#' containing them can be 100% identical to its own copy with the native
#' engine but not with BLAST. Coverage is the fraction of the query
#' covered by the best alignment to that database sequence. BLAST's `qcovs`
#' adds together several separate alignments to the same sequence, so
#' coverage values close to the threshold can differ slightly between the
#' engines.
#'
#' @section BLAST engine:
#' With `engine = "blast"`, BLAST+ (`makeblastdb` and `blastn`, default
#' megablast settings) must be installed and on the system `PATH`; see
#' [blast_available()]. Identity is `pident` of the best-scoring HSP and
#' coverage is `qcovs`, as in metacleaner 1.x.
#'
#' @param query Query sequences: a FASTA file path or a named character vector.
#' @param database Database sequences: a FASTA file path or a named character
#'   vector.
#' @param identity Minimum percent identity (0-100). Default 100.
#' @param coverage Minimum percent query coverage (0-100). Default 100.
#' @param max_hits Maximum number of qualifying hits reported per query
#'   (best first). Default 100.
#' @param best_hit If `TRUE` (default), queries without any qualifying hit
#'   also report their best-scoring non-qualifying hit (with `pass = FALSE`),
#'   when one exists.
#' @param exclude_self If `TRUE`, a query is not allowed to match a database
#'   sequence with the same ID. Default `FALSE`.
#' @param engine `"native"` (default) or `"blast"`.
#' @param threads Number of CPU threads to use, or `"auto"` to use all
#'   available cores but one. Default 1. The search is split across threads
#'   on all operating systems (including macOS); results are identical
#'   whatever the number of threads.
#' @param k Word length for the native engine index (8-13). Default 12.
#'   Smaller values are more sensitive at low identity but slower.
#' @param verbose Print progress messages. Default `TRUE`.
#'
#' @return A data frame with one row per reported hit and columns `query`,
#'   `subject` (database ID), `identity`, `coverage`, `score` (alignment score
#'   for the native engine, bit score for BLAST), `strand`, `qstart`, `qend`,
#'   `sstart`, `send`, `pass` (whether the hit meets both thresholds) and
#'   `n_pass` (total number of qualifying hits for that query, which can
#'   exceed `max_hits`).
#' @seealso [clean_reference()]
#' @examples
#' q <- c(q1 = "ACGTACGTTAGCATGCATCGATCGATGCTAGCTAGG")
#' db <- c(s1 = "TTTTACGTACGTTAGCATGCATCGATCGATGCTAGCTAGGTTTT",
#'         s2 = "ACGTACGTTAGCATGCATCGTTCGATGCTAGCTAGG")
#' search_sequences(q, db, verbose = FALSE)
#' search_sequences(q, db, identity = 95, coverage = 95, verbose = FALSE)
#' @export
search_sequences <- function(query, database, identity = 100, coverage = 100,
                             max_hits = 100, best_hit = TRUE, exclude_self = FALSE,
                             engine = c("native", "blast"), threads = 1, k = 12,
                             verbose = TRUE) {
  engine <- match.arg(engine)
  query <- as_sequences(query, "query")
  database <- as_sequences(database, "database")
  check_number(identity, "identity", 0, 100)
  check_number(coverage, "coverage", 0, 100)
  check_number(max_hits, "max_hits", 1)
  check_flag(best_hit, "best_hit")
  check_flag(exclude_self, "exclude_self")
  check_number(k, "k", 8, 13)
  threads <- check_threads(threads)
  if (engine == "blast") {
    return(search_blast(query, database, identity, coverage, max_hits, best_hit,
                        exclude_self, threads, verbose))
  }
  exclude <- if (exclude_self) match(names(query), names(database)) else integer(0)
  norm <- function(x) chartr("U", "T", toupper(unname(x)))
  res <- native_search_cpp(
    query = norm(query), subject = norm(database),
    identity = identity, coverage = coverage, k = as.integer(k),
    max_hits = as.integer(min(max_hits, .Machine$integer.max)),
    best_fail = best_hit, exclude = as.integer(exclude), threads = threads,
    occ_skip = 1000L, batch_size = 2000L, verbose = isTRUE(verbose))
  h <- res$hits
  # report minus-strand subject coordinates like BLAST (sstart > send)
  minus <- h$strand == "-"
  if (any(minus)) {
    tmp <- h$sstart[minus]
    h$sstart[minus] <- h$send[minus]
    h$send[minus] <- tmp
  }
  format_hits(h, names(query), names(database), res$n_pass)
}

format_hits <- function(h, qnames, snames, n_pass) {
  h$n_pass <- n_pass[h$query]
  h$query <- qnames[h$query]
  h$subject <- snames[h$subject]
  rownames(h) <- NULL
  h
}

empty_hits <- function() {
  data.frame(query = character(0), subject = character(0), identity = numeric(0),
             coverage = numeric(0), score = numeric(0), strand = character(0),
             qstart = integer(0), qend = integer(0), sstart = integer(0),
             send = integer(0), pass = logical(0), n_pass = integer(0),
             stringsAsFactors = FALSE)
}
