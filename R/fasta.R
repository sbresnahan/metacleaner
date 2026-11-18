#' Read sequences from a FASTA file
#'
#' Reads a FASTA file (plain text or gzip-compressed, with sequences on one
#' line or wrapped over many lines) into a named character vector. By default
#' each sequence is named by its accession, i.e. the first word of its header
#' line, with old-style NCBI headers such as `gi|12345|gb|AB123456.1|` reduced
#' to `AB123456.1`.
#'
#' @param file Path to a FASTA file. Files ending in `.gz` are decompressed
#'   automatically.
#' @param ids How sequences are named: `"accession"` (default) uses the
#'   accession (first word of the header); `"header"` keeps the full header
#'   line.
#' @param drop_version If `TRUE`, remove the version suffix from accessions
#'   (`AB123456.1` becomes `AB123456`). Default `FALSE`.
#'
#' @return A named character vector: the names are sequence IDs and the values
#'   are the sequences in upper case.
#' @seealso [write_fasta()], [dereplicate()]
#' @examples
#' f <- system.file("extdata", "query.fasta", package = "metacleaner")
#' seqs <- read_fasta(f)
#' length(seqs)
#' head(names(seqs))
#' @export
read_fasta <- function(file, ids = c("accession", "header"), drop_version = FALSE) {
  ids <- match.arg(ids)
  check_flag(drop_version, "drop_version")
  if (!is.character(file) || length(file) != 1L) {
    stop("`file` must be the path to one FASTA file.", call. = FALSE)
  }
  if (!file.exists(file)) stop("File not found: ", file, call. = FALSE)
  raw <- read_raw(file)
  parsed <- parse_fasta_raw(raw)
  seqs <- parsed$sequence
  nm <- if (ids == "accession") extract_accession(parsed$header) else parsed$header
  if (drop_version) nm <- base_accession(nm)
  names(seqs) <- nm
  empty <- !nzchar(seqs)
  if (any(empty)) {
    warning(sprintf("%d empty sequence(s) were removed from %s.", sum(empty), basename(file)),
            call. = FALSE)
    seqs <- seqs[!empty]
  }
  seqs
}

read_raw <- function(file) {
  is_gz <- grepl("\\.gz$", file, ignore.case = TRUE)
  if (!is_gz) {
    return(readBin(file, what = "raw", n = file.size(file)))
  }
  con <- gzfile(file, open = "rb")
  on.exit(close(con))
  chunks <- list()
  repeat {
    r <- readBin(con, what = "raw", n = 64L * 1024L * 1024L)
    if (length(r) == 0L) break
    chunks[[length(chunks) + 1L]] <- r
  }
  if (length(chunks) == 0L) raw(0) else do.call(c, chunks)
}

#' Write sequences to a FASTA file
#'
#' @param x A named character vector of sequences (for example from
#'   [read_fasta()]). The names are written as header lines.
#' @param file Output path. If it ends in `.gz` the file is gzip-compressed.
#' @param width Maximum number of bases per line. The default, `0`, writes
#'   each sequence on a single line, which is what most metabarcoding tools
#'   expect.
#'
#' @return The output path, invisibly.
#' @examples
#' seqs <- c(seq1 = "ACGTACGT", seq2 = "TTGACCA")
#' out <- tempfile(fileext = ".fasta")
#' write_fasta(seqs, out)
#' readLines(out)
#' @export
write_fasta <- function(x, file, width = 0) {
  if (!is.character(x) || (length(x) > 0L && is.null(names(x)))) {
    stop("`x` must be a named character vector of sequences.", call. = FALSE)
  }
  check_number(width, "width", 0)
  if (width > 0 && length(x) > 0L) {
    w <- as.integer(width)
    x <- vapply(x, function(s) {
      n <- nchar(s)
      if (n <= w) return(s)
      st <- seq.int(1L, n, by = w)
      paste(substring(s, st, pmin(st + w - 1L, n)), collapse = "\n")
    }, character(1), USE.NAMES = TRUE)
  }
  con <- if (grepl("\\.gz$", file, ignore.case = TRUE)) gzfile(file, "w") else file(file, "w")
  on.exit(close(con))
  if (length(x) > 0L) writeLines(paste0(">", names(x), "\n", x), con, useBytes = TRUE)
  invisible(file)
}

#' Extract accessions from FASTA header lines
#'
#' Returns the first word of each header. Old-style NCBI headers
#' (`gi|12345|gb|AB123456.1|`) and pipe-delimited headers (`gb|AB123456.1|`)
#' are reduced to the accession.
#'
#' @param header Character vector of FASTA header lines (with or without the
#'   leading `>`).
#' @return A character vector of accessions.
#' @examples
#' extract_accession(c(">MN123456.1 Malus domestica ITS2",
#'                     "gi|12345|gb|AB123456.1|"))
#' @export
extract_accession <- function(header) {
  h <- sub("^>", "", trimws(header))
  h <- sub("[[:space:]].*$", "", h)
  piped <- grepl("|", h, fixed = TRUE)
  if (any(piped)) {
    h[piped] <- vapply(strsplit(h[piped], "|", fixed = TRUE), function(tok) {
      tok <- tok[nzchar(tok)]
      if (length(tok) == 0L) return("")
      if (tok[1] == "gi" && length(tok) >= 4L) return(tok[4])
      if (length(tok) >= 2L) return(tok[2])
      tok[1]
    }, character(1))
  }
  h
}

#' Reverse complement DNA sequences
#'
#' @param x Character vector of DNA sequences (IUPAC ambiguity codes are
#'   complemented correctly).
#' @return A character vector of the same length, with names kept.
#' @examples
#' reverse_complement(c(a = "AACGTN", b = "RYKM"))
#' @export
reverse_complement <- function(x) {
  out <- revcomp_cpp(toupper(as.character(x)))
  names(out) <- names(x)
  out
}

#' Remove duplicate sequences
#'
#' Keeps the first copy of every distinct sequence, like `seqkit rmdup -s`.
#'
#' Be careful when de-replicating a *reference* ("good") database that is
#' used by [clean_reference()]: if two genera share an identical sequence,
#' keeping only one copy hides exactly the ambiguity that the
#' "matches multiple genus" check is designed to detect. De-replicating the
#' contaminant database is harmless.
#'
#' @param x A named character vector of sequences, or a path to a FASTA file.
#' @param both_strands If `TRUE` (default), a sequence and its reverse
#'   complement count as duplicates.
#' @return A named character vector with duplicates removed.
#' @examples
#' seqs <- c(a = "ACGTT", b = "ACGTT", c = "AACGT", d = "GGGCC")
#' dereplicate(seqs)
#' dereplicate(seqs, both_strands = FALSE)
#' @export
dereplicate <- function(x, both_strands = TRUE) {
  x <- as_sequences(x, "x")
  check_flag(both_strands, "both_strands")
  key <- toupper(x)
  if (both_strands) {
    rc <- revcomp_cpp(unname(key))
    key <- ifelse(key <= rc, key, rc)
  }
  x[!duplicated(key)]
}
