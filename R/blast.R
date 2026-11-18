#' Check whether BLAST+ is installed
#'
#' metacleaner does not need BLAST+, but it can use it as an alternative
#' search engine (`engine = "blast"`). This function reports whether the
#' `makeblastdb` and `blastn` programs can be found on the system `PATH`.
#'
#' @param verbose If `TRUE`, print where the programs were found or how to
#'   install them.
#' @return `TRUE` or `FALSE`, invisibly when `verbose = TRUE`.
#' @examples
#' blast_available()
#' @export
blast_available <- function(verbose = FALSE) {
  paths <- Sys.which(c("makeblastdb", "blastn"))
  ok <- all(nzchar(paths))
  if (isTRUE(verbose)) {
    if (ok) {
      message("BLAST+ found:\n  ", paste(paths, collapse = "\n  "))
    } else {
      message("BLAST+ was not found. The default engine = \"native\" does not need it.\n",
              "To use engine = \"blast\", install BLAST+ from ",
              "https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ ",
              "and make sure 'blastn' and 'makeblastdb' are on your PATH.")
    }
    return(invisible(ok))
  }
  ok
}

run_tool <- function(cmd, args) {
  out <- suppressWarnings(system2(Sys.which(cmd), args, stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  if (!is.null(status) && status != 0) {
    stop(sprintf("%s failed (exit status %s):\n%s", cmd, status,
                 paste(utils::tail(out, 20), collapse = "\n")), call. = FALSE)
  }
  invisible(out)
}

search_blast <- function(query, database, identity, coverage, max_hits, best_hit,
                         exclude_self, threads, verbose) {
  if (!blast_available()) {
    stop("engine = \"blast\" needs BLAST+, which was not found. ",
         "Run blast_available(verbose = TRUE) for help, or use engine = \"native\".",
         call. = FALSE)
  }
  td <- tempfile("metacleaner_blast_")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  qf <- file.path(td, "query.fasta")
  sf <- file.path(td, "db.fasta")
  db <- file.path(td, "db")
  out <- file.path(td, "hits.tsv")
  # simple local IDs avoid problems with unusual accession formats
  # (IDs such as "s10000" would be parsed as GenBank accessions by
  # -parse_seqids, so plain numbered IDs are used without that option)
  write_fasta(stats_set_names(query, paste0("mcq_", seq_along(query))), qf)
  write_fasta(stats_set_names(database, paste0("mcs_", seq_along(database))), sf)
  msg(verbose, "  building BLAST database (", length(database), " sequences)")
  run_tool("makeblastdb", c("-in", shQuote(sf), "-dbtype", "nucl", "-out", shQuote(db)))
  msg(verbose, "  running blastn (", length(query), " queries, ", threads, " thread(s))")
  ntarget <- max(as.integer(min(max_hits, 1e6)), 10L)
  run_tool("blastn", c("-query", shQuote(qf), "-db", shQuote(db), "-out", shQuote(out),
                       "-outfmt", shQuote("6 qseqid sseqid pident qcovs length qstart qend sstart send bitscore sstrand"),
                       "-max_target_seqs", ntarget, "-num_threads", threads))
  if (!file.exists(out) || file.size(out) == 0) return(empty_hits())
  h <- utils::read.delim(out, header = FALSE, stringsAsFactors = FALSE,
                         col.names = c("qseqid", "sseqid", "pident", "qcovs", "length",
                                       "qstart", "qend", "sstart", "send", "bitscore",
                                       "sstrand"))
  qi <- as.integer(sub("^.*mcq_([0-9]+).*$", "\\1", h$qseqid))
  si <- as.integer(sub("^.*mcs_([0-9]+).*$", "\\1", h$sseqid))
  if (anyNA(qi) || anyNA(si)) stop("Could not parse BLAST output IDs.", call. = FALSE)
  hits <- data.frame(query = qi, subject = si, identity = h$pident, coverage = h$qcovs,
                     score = h$bitscore,
                     strand = ifelse(h$sstrand == "minus", "-", "+"),
                     qstart = h$qstart, qend = h$qend, sstart = h$sstart, send = h$send,
                     stringsAsFactors = FALSE)
  if (exclude_self) {
    hits <- hits[names(query)[hits$query] != names(database)[hits$subject], , drop = FALSE]
  }
  # best HSP per query-subject pair
  hits <- hits[order(hits$query, -hits$score, -hits$identity, hits$subject), , drop = FALSE]
  hits <- hits[!duplicated(hits[, c("query", "subject")]), , drop = FALSE]
  hits$pass <- hits$identity + 1e-9 >= identity & hits$coverage + 1e-9 >= coverage
  n_pass <- tabulate(hits$query[hits$pass], nbins = length(query))
  keep_pass <- hits$pass & stats_ave_seq(hits$query, hits$pass) <= max_hits
  keep <- keep_pass
  if (best_hit) {
    nopass <- n_pass[hits$query] == 0L
    first <- !duplicated(hits$query)
    keep <- keep | (nopass & first)
  }
  hits <- hits[keep, , drop = FALSE]
  format_hits(hits, names(query), names(database), n_pass)
}

stats_set_names <- function(x, nm) {
  names(x) <- nm
  x
}

# running count of TRUE values within each group (rows already ordered)
stats_ave_seq <- function(g, flag) {
  out <- integer(length(g))
  if (length(g) == 0L) return(out)
  cs <- stats::ave(as.integer(flag), g, FUN = cumsum)
  as.integer(cs)
}
