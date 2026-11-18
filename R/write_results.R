#' Save metacleaner results to files
#'
#' Writes the cleaned database and a report of the removed sequences:
#'
#' * `<prefix>_clean.fasta`: the kept sequences (one line per sequence);
#' * `<prefix>_clean.tax`: taxonomy of the kept sequences;
#' * `<prefix>_flagged.tsv`: one row per removed sequence with its taxonomy,
#'   the hit that caused it to be flagged, that hit's taxonomy and the reason;
#' * `<prefix>_summary.txt`: settings and counts per reason.
#'
#' @param result The output of [clean_reference()].
#' @param outdir Folder to write to (created if needed), e.g. `"cleaned"`.
#' @param prefix File name prefix. Default: the prefix stored in `result`.
#' @param tax_format Format of the `_clean.tax` file: `"taxonomizr"`
#'   (default; the comma-separated format used by metacleaner 1.x and as
#'   MetaCurator input: `taxid,accession,taxid,domain,...,species`) or
#'   `"lineage"` (tab-separated `accession` and a `k__;p__;...;s__` lineage
#'   string, as used by QIIME-style classifiers).
#' @param overwrite Overwrite existing files? Default `TRUE`.
#' @param verbose Print the names of the written files. Default `TRUE`.
#' @return The paths of the written files, invisibly.
#' @examples
#' ex <- function(f) system.file("extdata", f, package = "metacleaner")
#' res <- clean_reference(ex("query.fasta"), ex("reference.fasta"),
#'                        ex("contaminants.fasta"), taxonomy = ex("taxonomy.tax"),
#'                        verbose = FALSE)
#' files <- write_results(res, outdir = tempdir(), prefix = "example")
#' files
#' @export
write_results <- function(result, outdir, prefix = NULL,
                          tax_format = c("taxonomizr", "lineage"), overwrite = TRUE,
                          verbose = TRUE) {
  if (!inherits(result, "metacleaner_result")) {
    stop("`result` must be the output of clean_reference().", call. = FALSE)
  }
  tax_format <- match.arg(tax_format)
  if (missing(outdir) || !is.character(outdir) || length(outdir) != 1L) {
    stop("Please give the folder to save the results in, e.g. outdir = \"cleaned\".",
         call. = FALSE)
  }
  if (is.null(prefix)) prefix <- result$settings$prefix
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  f <- file.path(outdir, paste0(prefix, c("_clean.fasta", "_clean.tax", "_flagged.tsv",
                                          "_summary.txt")))
  names(f) <- c("fasta", "taxonomy", "flagged", "summary")
  if (!overwrite && any(file.exists(f))) {
    stop("Output files already exist in ", outdir, "; use overwrite = TRUE or another prefix.",
         call. = FALSE)
  }
  write_fasta(result$clean, f[["fasta"]])
  tx <- result$taxonomy
  if (tax_format == "taxonomizr") {
    d <- data.frame(tx$taxid, tx$id, tx$taxid, tx[, taxonomic_ranks()],
                    stringsAsFactors = FALSE)
    utils::write.table(d, f[["taxonomy"]], sep = ",", row.names = FALSE, col.names = FALSE)
  } else {
    pre <- c("k__", "p__", "c__", "o__", "f__", "g__", "s__")
    m <- as.matrix(tx[, taxonomic_ranks(), drop = FALSE])
    m[is.na(m)] <- ""
    lin <- if (nrow(m)) {
      apply(m, 1L, function(v) paste0(pre, v, collapse = ";"))
    } else character(0)
    utils::write.table(data.frame(tx$id, lin), f[["taxonomy"]], sep = "\t", quote = FALSE,
                       row.names = FALSE, col.names = FALSE)
  }
  utils::write.table(result$flagged, f[["flagged"]], sep = "\t", quote = FALSE,
                     row.names = FALSE, na = "NA")
  s <- result$settings
  sm <- summary(result)
  lines <- c(
    sprintf("metacleaner %s  (%s)", s$version, format(Sys.time(), "%Y-%m-%d %H:%M")),
    "",
    sprintf("query sequences:       %d", s$n_query),
    sprintf("reference sequences:   %d", s$n_reference),
    sprintf("contaminant sequences: %d", s$n_contaminants),
    sprintf("identity >= %s%%, coverage >= %s%%, level = %s, rescue = %s, exclude_self = %s",
            s$identity, s$coverage, s$level, s$rescue, s$exclude_self),
    if (!is.null(s$require)) sprintf("require: %s = %s", names(s$require), s$require),
    sprintf("engine: %s", s$engine),
    "",
    sprintf("kept:    %d", sm$n_kept),
    sprintf("flagged: %d", sm$n_flagged),
    if (nrow(sm$reasons)) sprintf("  %-40s %d", sm$reasons$reason, sm$reasons$n)
  )
  writeLines(lines, f[["summary"]])
  msg(verbose, "Results written to ", normalizePath(outdir), ":\n  ",
      paste(basename(f), collapse = "\n  "))
  invisible(f)
}
