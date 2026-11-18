#' Clean a barcode reference database
#'
#' The main metacleaner function. It removes sequences that are probably
#' mislabeled or contaminants from a DNA barcode reference database (for
#' example a plant ITS2 database built from NCBI GenBank records). For each
#' query sequence it:
#'
#' 1. searches a database of unwanted sequences (`contaminants`, e.g.
#'    fungal ITS) and flags queries that match one;
#' 2. searches a database of wanted sequences (`reference`, e.g. plant ITS)
#'    and flags queries without a match;
#' 3. compares the taxonomy of each query with the taxonomy of its matches at
#'    the chosen `level` (genus by default) and flags queries whose matches
#'    belong to a different or to several taxa.
#'
#' A match ("hit") counts only if it reaches both `identity` and `coverage`.
#' The full decision rules are listed in [classify_queries()].
#'
#' @param query The database to clean: a FASTA file path or a named character
#'   vector (see [read_fasta()]). Sequence names must be accessions that can be
#'   found in `taxonomy`.
#' @param reference Sequences of the wanted group (e.g. all plant ITS
#'   sequences from NCBI): FASTA file path or named character vector. It may
#'   contain the query sequences themselves.
#' @param contaminants Optional sequences of unwanted groups (e.g. fungal
#'   ITS): FASTA file path or named character vector. `NULL` (default) skips
#'   this step.
#' @param taxonomy Where to find the taxonomy of query and reference
#'   sequences. One of:
#'   * a taxonomy file path (see [read_taxonomy()] for the formats),
#'   * `"ncbi"` to look up accessions online (see [taxonomy_from_ncbi()]),
#'   * the path to a taxonomizr database ending in `.sql`
#'     (see [taxonomy_from_taxonomizr()]),
#'   * a lineage table (data frame) such as returned by [read_taxonomy()],
#'   * or a list / vector of several of these, tried in order: e.g.
#'     `c("query.tax", "ncbi")` uses the file first and NCBI for anything
#'     not in the file.
#' @param identity Minimum percent identity for a hit (0-100). Default 100.
#' @param coverage Minimum percent of the query covered by a hit (0-100).
#'   Default 100.
#' @param level Taxonomic rank at which a query and its hits must agree: one
#'   of `"domain"`, `"phylum"`, `"class"`, `"order"`, `"family"`, `"genus"`
#'   (default) or `"species"`.
#' @param require Optional extra filter, e.g. `c(class = "Magnoliopsida")`:
#'   queries with any hit outside this taxon are flagged. (The old
#'   command-line form `"class,Magnoliopsida"` is also accepted.)
#' @param rescue If `TRUE` (default, as in metacleaner 1.x), a query with no
#'   hit at the thresholds is kept when its best lower-scoring hit belongs to
#'   the same taxon at `level`. Set to `FALSE` for stricter cleaning.
#' @param exclude_self If `TRUE`, a query may not match the reference sequence
#'   with the same accession. Default `FALSE`.
#' @param engine Search engine: `"native"` (default, built in) or `"blast"`
#'   (requires BLAST+; see [blast_available()]).
#' @param threads Number of CPU threads to use, or `"auto"` to use all
#'   available cores but one. Default 1. The search is split across threads
#'   on all operating systems (including macOS); results are identical
#'   whatever the number of threads.
#' @param max_hits Maximum hits per query used in the taxonomy checks
#'   (best-scoring first). Default 100.
#' @param outdir Optional output folder. If given, the results are also
#'   written there with [write_results()].
#' @param prefix File name prefix for the output files. Defaults to the query
#'   file name without extension, or `"metacleaner"`.
#' @param verbose Print progress messages. Default `TRUE`.
#'
#' @return An object of class `metacleaner_result`: a list with
#'   * `clean`: the kept sequences (named character vector),
#'   * `flagged`: a data frame with one row per removed sequence and the
#'     reason, the hit that triggered it and both taxonomies,
#'   * `status`: one row per query sequence (kept or not),
#'   * `taxonomy`: lineage table of the kept sequences,
#'   * `hits`: the search results (`reference` and `contaminants`),
#'   * `settings`: the settings used.
#'
#'   Print it for a summary; use [write_results()] to save it.
#' @seealso [write_results()], [read_fasta()], [read_taxonomy()],
#'   [build_reference_db()]
#' @examples
#' ex <- function(f) system.file("extdata", f, package = "metacleaner")
#' res <- clean_reference(
#'   query        = ex("query.fasta"),
#'   reference    = ex("reference.fasta"),
#'   contaminants = ex("contaminants.fasta"),
#'   taxonomy     = ex("taxonomy.tax"),
#'   level        = "genus"
#' )
#' res
#' head(res$flagged)
#' \donttest{
#' # save the cleaned database and the report
#' write_results(res, outdir = tempdir(), prefix = "example")
#' }
#' @export
clean_reference <- function(query, reference, contaminants = NULL, taxonomy,
                            identity = 100, coverage = 100, level = "genus",
                            require = NULL, rescue = TRUE, exclude_self = FALSE,
                            engine = c("native", "blast"), threads = 1, max_hits = 100,
                            outdir = NULL, prefix = NULL, verbose = TRUE) {
  t0 <- Sys.time()
  engine <- match.arg(engine)
  if (missing(taxonomy)) {
    stop("Please say where to find the taxonomy, e.g. taxonomy = \"my_database.tax\" ",
         "or taxonomy = \"ncbi\". See ?clean_reference.", call. = FALSE)
  }
  level <- normalise_rank(level, "level")
  check_require(require)
  check_number(identity, "identity", 0, 100)
  check_number(coverage, "coverage", 0, 100)
  if (is.null(prefix)) {
    prefix <- if (is.character(query) && length(query) == 1L && is.null(names(query))) {
      sub("\\.(fa|fasta|fna|fas)(\\.gz)?$", "", basename(query), ignore.case = TRUE)
    } else "metacleaner"
  }
  sources <- resolve_taxonomy_sources(taxonomy)

  msg(verbose, "metacleaner: reading sequences")
  q <- as_sequences(query, "query")
  ref <- as_sequences(reference, "reference")
  con <- as_sequences(contaminants, "contaminants")
  msg(verbose, "  ", length(q), " query, ", length(ref), " reference",
      if (!is.null(con)) paste0(", ", length(con), " contaminant") else "", " sequences")

  ch <- NULL
  if (!is.null(con)) {
    msg(verbose, "Step 1/3: searching the contaminant database")
    ch <- search_sequences(q, con, identity = identity, coverage = coverage, max_hits = 1,
                           best_hit = FALSE, engine = engine, threads = threads,
                           verbose = verbose)
  }
  bad <- if (is.null(ch)) character(0) else unique(ch$query[ch$pass])
  q2 <- q[!names(q) %in% bad]
  msg(verbose, "Step 2/3: searching the reference database")
  rh <- if (length(q2)) {
    search_sequences(q2, ref, identity = identity, coverage = coverage, max_hits = max_hits,
                     best_hit = rescue, exclude_self = exclude_self, engine = engine,
                     threads = threads, verbose = verbose)
  } else empty_hits()
  msg(verbose, "Step 3/3: comparing taxonomy at ", level, " level")
  cl <- classify_queries(names(q), rh, ch, taxonomy = sources, level = level,
                         require = require, rescue = rescue, verbose = verbose)
  st <- cl$status
  lin <- cl$lineages
  keep <- st$id[!st$flagged]
  fl <- st[st$flagged, , drop = FALSE]
  flagged <- data.frame(
    id = fl$id,
    query_taxonomy = lineage_strings(fl$id, lin),
    top_hit = fl$top_hit,
    hit_taxonomy = lineage_strings(fl$top_hit, lin),
    reason = fl$reason,
    database = fl$database,
    identity = fl$identity,
    coverage = fl$coverage,
    n_hits = fl$n_hits,
    hit_taxa = fl$hit_taxa,
    stringsAsFactors = FALSE)
  flagged <- flagged[order(flagged$id), , drop = FALSE]
  rownames(flagged) <- NULL
  tax_keep <- lin[match(keep, lin$id), , drop = FALSE]
  tax_keep$id <- keep
  rownames(tax_keep) <- NULL

  res <- structure(list(
    clean = q[keep],
    flagged = flagged,
    status = st,
    taxonomy = tax_keep,
    hits = list(reference = rh, contaminants = ch),
    settings = list(identity = identity, coverage = coverage, level = level,
                    require = require, rescue = rescue, exclude_self = exclude_self,
                    engine = engine, max_hits = max_hits, n_query = length(q),
                    n_reference = length(ref),
                    n_contaminants = if (is.null(con)) 0L else length(con),
                    prefix = prefix, runtime = difftime(Sys.time(), t0, units = "secs"),
                    version = as.character(utils::packageVersion("metacleaner")))
  ), class = "metacleaner_result")
  msg(verbose, sprintf("Done: kept %d of %d sequences, flagged %d (%.1f s).",
                       length(keep), length(q), nrow(flagged),
                       as.numeric(res$settings$runtime)))
  if (!is.null(outdir)) write_results(res, outdir = outdir, prefix = prefix, verbose = verbose)
  res
}

lineage_strings <- function(ids, lin) {
  out <- rep(NA_character_, length(ids))
  m <- match(ids, lin$id)
  ok <- !is.na(m)
  if (any(ok)) out[ok] <- collapse_lineage(lin[m[ok], , drop = FALSE])
  out
}

#' @export
print.metacleaner_result <- function(x, ...) {
  s <- x$settings
  cat("<metacleaner result>\n")
  cat(sprintf("  query sequences:   %d\n", s$n_query))
  cat(sprintf("  kept:              %d\n", length(x$clean)))
  cat(sprintf("  flagged:           %d\n", nrow(x$flagged)))
  if (nrow(x$flagged)) {
    tab <- sort(table(x$flagged$reason), decreasing = TRUE)
    for (i in seq_along(tab)) cat(sprintf("    %-40s %d\n", names(tab)[i], tab[[i]]))
  }
  cat(sprintf("  settings: identity >= %s%%, coverage >= %s%%, level = %s, engine = %s%s\n",
              s$identity, s$coverage, s$level, s$engine,
              if (!is.null(s$require)) paste0(", require ", names(s$require), " = ", s$require) else ""))
  cat("Use write_results() to save the cleaned database and report.\n")
  invisible(x)
}

#' @export
summary.metacleaner_result <- function(object, ...) {
  s <- object$settings
  reasons <- if (nrow(object$flagged)) {
    as.data.frame(table(reason = object$flagged$reason), stringsAsFactors = FALSE)
  } else data.frame(reason = character(0), Freq = integer(0))
  names(reasons)[2] <- "n"
  reasons <- reasons[order(-reasons$n), , drop = FALSE]
  rownames(reasons) <- NULL
  out <- list(n_query = s$n_query, n_kept = length(object$clean),
              n_flagged = nrow(object$flagged), reasons = reasons, settings = s)
  class(out) <- "summary.metacleaner_result"
  out
}

#' @export
print.summary.metacleaner_result <- function(x, ...) {
  cat(sprintf("metacleaner: kept %d of %d sequences (%.1f%%); flagged %d\n",
              x$n_kept, x$n_query, 100 * x$n_kept / max(1, x$n_query), x$n_flagged))
  if (nrow(x$reasons)) print(x$reasons, row.names = FALSE)
  invisible(x)
}
