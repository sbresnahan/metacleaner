#' Download sequences from NCBI GenBank
#'
#' Downloads all NCBI Nucleotide records that match a search term (written
#' exactly as you would type it on the NCBI website) to a FASTA file, in
#' batches. Requires an internet connection and the suggested package
#' \pkg{rentrez}.
#'
#' Large downloads can take a long time. If a download is interrupted, run
#' the same command again with `resume = TRUE` to continue where it stopped.
#' An NCBI API key (free, see <https://account.ncbi.nlm.nih.gov/>)
#' makes downloads faster.
#'
#' @param term NCBI Nucleotide search term, e.g.
#'   `'(ITS2) AND Embryophyta[Organism] AND ("0"[SLEN] : "10000"[SLEN])'`.
#' @param file Output FASTA file.
#' @param batch_size Records per request (max 10000). Default 500.
#' @param max_records Maximum number of records to download. Default: all.
#' @param resume If `TRUE` and `file` exists, continue an interrupted download
#'   by appending to it. If `FALSE` (default) an existing file is an error.
#' @param api_key Optional NCBI API key. Defaults to the `ENTREZ_KEY`
#'   environment variable.
#' @param verbose Print progress. Default `TRUE`.
#' @return The number of records in `file`, invisibly.
#' @seealso [build_reference_db()]
#' @examples
#' \donttest{
#' if (requireNamespace("rentrez", quietly = TRUE)) {
#'   out <- tempfile(fileext = ".fasta")
#'   try(download_ncbi_sequences("Malus[Organism] AND ITS2", out, max_records = 20))
#' }
#' }
#' @export
download_ncbi_sequences <- function(term, file, batch_size = 500, max_records = Inf,
                                    resume = FALSE, api_key = Sys.getenv("ENTREZ_KEY"),
                                    verbose = TRUE) {
  need_pkg("rentrez", "download sequences from NCBI")
  check_number(batch_size, "batch_size", 1, 10000)
  check_number(max_records, "max_records", 1)
  check_flag(resume, "resume")
  if (nzchar(api_key)) rentrez::set_entrez_key(api_key)
  done <- 0L
  if (file.exists(file)) {
    if (!resume) {
      stop("The file ", file, " already exists. Use resume = TRUE to continue a download, ",
           "or choose another file name.", call. = FALSE)
    }
    done <- count_fasta_records(file)
  }
  s <- with_retry(function() {
    rentrez::entrez_search(db = "nuccore", term = term, use_history = TRUE, retmax = 0)
  })
  total <- min(as.numeric(s$count), max_records)
  msg(verbose, "NCBI search found ", s$count, " record(s); downloading ", total,
      if (done > 0) paste0(" (", done, " already in file)") else "")
  if (total <= done) return(invisible(done))
  con <- file(file, open = if (done > 0) "a" else "w")
  on.exit(close(con))
  start <- done
  while (start < total) {
    n <- as.integer(min(batch_size, total - start))
    txt <- with_retry(function() {
      rentrez::entrez_fetch(db = "nuccore", web_history = s$web_history, rettype = "fasta",
                            retmode = "text", retstart = start, retmax = n)
    })
    txt <- sub("\n+$", "", txt)
    if (nzchar(txt)) writeLines(txt, con)
    flush(con)
    start <- start + n
    msg(verbose, sprintf("  %d / %d", start, as.integer(total)))
    ncbi_pause()
  }
  invisible(count_fasta_records(file))
}

count_fasta_records <- function(file) {
  n <- 0L
  con <- file(file, "r")
  on.exit(close(con))
  repeat {
    l <- readLines(con, n = 100000L, warn = FALSE)
    if (length(l) == 0L) break
    n <- n + sum(startsWith(l, ">"))
  }
  n
}

#' Build a reference (or contaminant) database from NCBI
#'
#' A one-step replacement for the shell commands previously needed to build
#' metacleaner databases: downloads all NCBI Nucleotide records that match
#' `term`, simplifies the headers to accessions, removes empty and
#' (optionally) duplicate sequences and writes one sequence per line.
#'
#' @inheritParams download_ncbi_sequences
#' @param file Output FASTA file for the finished database.
#' @param dereplicate Remove duplicate sequences (see [dereplicate()])?
#'   Default `FALSE`. Recommended for contaminant databases; **not**
#'   recommended for the reference database, because identical sequences from
#'   different genera are what reveals ambiguous or mislabeled entries.
#' @param min_length Remove sequences shorter than this. Default 0.
#' @param keep_download Keep the raw download next to `file` (with suffix
#'   `.raw.fasta`)? Default `FALSE`.
#' @return The path to `file`, invisibly.
#' @examples
#' \dontrun{
#' # Plant ITS1/ITS2 reference database (as used by Crone et al. 2023)
#' build_reference_db(
#'   '((((ITS1) OR 5.8S) OR 28S) OR ITS2) AND Embryophyta[Organism] AND ("0"[SLEN] : "10000"[SLEN])',
#'   "ITS_plants.fasta")
#' # Non-plant (contaminant) database
#' build_reference_db(
#'   '((((ITS1) OR 5.8S) OR 28S) OR ITS2) NOT Embryophyta[Organism] AND ("0"[SLEN] : "10000"[SLEN])',
#'   "ITS_nonplants.fasta", dereplicate = TRUE)
#' }
#' @export
build_reference_db <- function(term, file, dereplicate = FALSE, min_length = 0,
                               max_records = Inf, keep_download = FALSE,
                               api_key = Sys.getenv("ENTREZ_KEY"), verbose = TRUE) {
  check_flag(dereplicate, "dereplicate")
  check_number(min_length, "min_length", 0)
  raw_file <- paste0(sub("\\.(fa|fasta|fna|fas)$", "", file), ".raw.fasta")
  download_ncbi_sequences(term, raw_file, max_records = max_records,
                          resume = file.exists(raw_file), api_key = api_key,
                          verbose = verbose)
  seqs <- read_fasta(raw_file)
  n0 <- length(seqs)
  seqs <- seqs[!duplicated(names(seqs))]
  seqs <- seqs[nchar(seqs) >= min_length]
  if (length(seqs) == 0) {
    stop("None of the ", n0, " downloaded sequences is at least ", min_length,
         " bases long, so nothing was written. Lower min_length or check the search term. ",
         "The download was kept in ", raw_file, ".", call. = FALSE)
  }
  if (dereplicate) seqs <- dereplicate(seqs)
  write_fasta(seqs, file)
  if (!keep_download) unlink(raw_file)
  msg(verbose, sprintf("Wrote %d of %d downloaded sequences to %s", length(seqs), n0, file))
  invisible(file)
}
