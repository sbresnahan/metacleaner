#' Read a taxonomy file
#'
#' Reads a file that gives the taxonomy (lineage) of each sequence and
#' returns it as a standard lineage table. Three layouts are recognized
#' automatically:
#'
#' * **taxonomizr / MetaCurator input format** (also written by
#'   [write_results()]): comma-separated, no header, columns `taxid`,
#'   `accession`, `taxid`, then domain, phylum, class, order, family, genus and
#'   species. Example:
#'   `"49992","Z37472.1","49992","Eukaryota","Streptophyta",NA,"Lamiales",...`
#' * **lineage strings** (QIIME, Metaxa2 and similar): two tab-separated
#'   columns, the ID and a `;`-separated lineage, optionally with rank prefixes
#'   such as `k__`, `d__`, `p__`, `c__`, `o__`, `f__`, `g__`, `s__`.
#' * **a table with a header row** (tab- or comma-separated) with an ID column
#'   (`id`, `accession`, `seqid` or `qseqid`) and columns named after ranks
#'   (`superkingdom`/`kingdom`/`domain`, `phylum`, ..., `species`).
#'
#' @param file Path to the taxonomy file.
#' @param format `"auto"` (default), `"taxonomizr"`, `"lineage"` or `"table"`.
#' @return A lineage table: a data frame with columns `id`, `taxid` and the
#'   seven ranks returned by [taxonomic_ranks()]; missing values are `NA`.
#' @examples
#' f <- system.file("extdata", "taxonomy.tax", package = "metacleaner")
#' tax <- read_taxonomy(f)
#' head(tax)
#' @export
read_taxonomy <- function(file, format = c("auto", "taxonomizr", "lineage", "table")) {
  format <- match.arg(format)
  if (!file.exists(file)) stop("Taxonomy file not found: ", file, call. = FALSE)
  first <- readLines(file, n = 5L, warn = FALSE)
  first <- first[nzchar(first)]
  if (length(first) == 0L) return(empty_lineages())
  if (format == "auto") format <- guess_tax_format(first)
  out <- switch(format,
    taxonomizr = read_tax_taxonomizr(file),
    lineage = read_tax_lineage(file),
    table = read_tax_table(file)
  )
  out <- out[!is.na(out$id) & nzchar(out$id), , drop = FALSE]
  out <- out[!duplicated(out$id), , drop = FALSE]
  rownames(out) <- NULL
  out
}

guess_tax_format <- function(lines) {
  l1 <- lines[1]
  if (grepl("\t", l1)) {
    f <- strsplit(l1, "\t", fixed = TRUE)[[1]]
    if (length(f) == 2L && (grepl(";", f[2]) || !tolower(f[1]) %in% c("id", "accession"))) {
      return("lineage")
    }
    return("table")
  }
  f <- utils::read.csv(text = l1, header = FALSE, stringsAsFactors = FALSE,
                       colClasses = "character")
  if (ncol(f) == 10L && grepl("^ *[0-9]+ *$|^NA$", f[[1]])) return("taxonomizr")
  "table"
}

empty_lineages <- function() {
  out <- data.frame(id = character(0), taxid = character(0), stringsAsFactors = FALSE)
  for (r in taxonomic_ranks()) out[[r]] <- character(0)
  out
}

make_lineages <- function(id, taxid, ranks) {
  out <- data.frame(id = as.character(id), taxid = as.character(taxid),
                    stringsAsFactors = FALSE)
  for (r in taxonomic_ranks()) {
    v <- if (!is.null(ranks[[r]])) as.character(ranks[[r]]) else rep(NA_character_, length(id))
    v <- trimws(v)
    v[v %in% c("", "NA", "na", "N/A", "unclassified", "Unclassified")] <- NA_character_
    out[[r]] <- v
  }
  out$taxid <- trimws(out$taxid)
  out$taxid[out$taxid %in% c("", "NA")] <- NA_character_
  out
}

read_tax_taxonomizr <- function(file) {
  d <- utils::read.csv(file, header = FALSE, stringsAsFactors = FALSE,
                       colClasses = "character", na.strings = c("NA", ""))
  if (ncol(d) < 10L) {
    stop("Expected 10 comma-separated columns (taxid, accession, taxid, 7 ranks) in ",
         file, call. = FALSE)
  }
  ranks <- stats::setNames(as.list(d[, 4:10]), taxonomic_ranks())
  make_lineages(d[[2]], d[[1]], ranks)
}

read_tax_lineage <- function(file) {
  d <- utils::read.delim(file, header = FALSE, stringsAsFactors = FALSE,
                         colClasses = "character", quote = "", comment.char = "")
  if (tolower(d[1, 1]) %in% c("id", "accession", "feature id", "seqid")) d <- d[-1, , drop = FALSE]
  parts <- strsplit(d[[2]], ";", fixed = TRUE)
  prefix_map <- c(k = "domain", d = "domain", p = "phylum", c = "class", o = "order",
                  f = "family", g = "genus", s = "species")
  ranks <- lapply(stats::setNames(taxonomic_ranks(), taxonomic_ranks()),
                  function(r) rep(NA_character_, nrow(d)))
  for (i in seq_along(parts)) {
    p <- trimws(parts[[i]])
    p <- p[nzchar(p)]
    if (length(p) == 0L) next
    has_prefix <- grepl("^[a-z]__", p)
    if (all(has_prefix)) {
      key <- prefix_map[substr(p, 1L, 1L)]
      val <- sub("^[a-z]__", "", p)
      ok <- !is.na(key)
      for (j in which(ok)) ranks[[key[j]]][i] <- val[j]
    } else {
      p <- sub("^[a-z]__", "", p)
      n <- min(length(p), 7L)
      for (j in seq_len(n)) ranks[[j]][i] <- p[j]
    }
  }
  make_lineages(d[[1]], NA_character_, ranks)
}

read_tax_table <- function(file) {
  l1 <- readLines(file, n = 1L, warn = FALSE)
  sep <- if (grepl("\t", l1)) "\t" else ","
  d <- utils::read.table(file, header = TRUE, sep = sep, stringsAsFactors = FALSE,
                         colClasses = "character", quote = "\"", comment.char = "",
                         check.names = FALSE, na.strings = c("NA", ""))
  nm <- tolower(trimws(names(d)))
  id_col <- which(nm %in% c("id", "accession", "seqid", "qseqid", "sequence_id"))[1]
  if (is.na(id_col)) {
    stop("Could not find an ID column (named id, accession or seqid) in ", file, call. = FALSE)
  }
  nm[nm %in% c("superkingdom", "kingdom")] <- "domain"
  ranks <- list()
  for (r in taxonomic_ranks()) {
    j <- which(nm == r)[1]
    if (!is.na(j)) ranks[[r]] <- d[[j]]
  }
  if (length(ranks) == 0L) {
    stop("No rank columns (domain/superkingdom, phylum, ..., species) found in ", file,
         call. = FALSE)
  }
  tj <- which(nm %in% c("taxid", "taxonomy_id", "tax_id"))[1]
  make_lineages(d[[id_col]], if (is.na(tj)) NA_character_ else d[[tj]], ranks)
}

#' Look up taxonomy for accessions from NCBI
#'
#' Looks up the lineage of each accession with the NCBI E-utilities web
#' service (no local database needed). Only the accessions you ask for are
#' downloaded, in batches. Requires an internet connection and the
#' suggested packages \pkg{rentrez} and \pkg{xml2}.
#'
#' NCBI limits anonymous users to 3 requests per second; with a free API key
#' (see <https://account.ncbi.nlm.nih.gov/>) the limit is 10 per second.
#' Use `cache` to store results in a file so that repeated runs do not
#' download the same information again.
#'
#' @param ids Character vector of nucleotide accessions.
#' @param api_key Optional NCBI API key. Defaults to the `ENTREZ_KEY`
#'   environment variable if set.
#' @param cache Optional path to an `.rds` file used to store and reuse
#'   lookups between runs.
#' @param batch_size Number of IDs per web request. Default 200.
#' @param verbose Print progress messages. Default `TRUE`.
#' @return A lineage table (see [read_taxonomy()]) with one row per ID that
#'   could be resolved.
#' @examples
#' \donttest{
#' if (requireNamespace("rentrez", quietly = TRUE) &&
#'     requireNamespace("xml2", quietly = TRUE)) {
#'   try(taxonomy_from_ncbi(c("MW545606.1", "MZ509176.1")))
#' }
#' }
#' @export
taxonomy_from_ncbi <- function(ids, api_key = Sys.getenv("ENTREZ_KEY"), cache = NULL,
                               batch_size = 200, verbose = TRUE) {
  need_pkg("rentrez", "look up taxonomy from NCBI")
  need_pkg("xml2", "look up taxonomy from NCBI")
  check_number(batch_size, "batch_size", 1, 500)
  ids <- unique(as.character(ids[!is.na(ids) & nzchar(ids)]))
  if (nzchar(api_key)) rentrez::set_entrez_key(api_key)
  cached <- empty_lineages()
  if (!is.null(cache) && file.exists(cache)) cached <- readRDS(cache)
  todo <- ids[is.na(match_ids(ids, cached$id))]
  if (length(todo) > 0L) {
    msg(verbose, "  looking up ", length(todo), " accession(s) at NCBI")
    # work in chunks and update the cache after each one, so an interrupted
    # lookup can be resumed without starting again
    chunks <- split(todo, ceiling(seq_along(todo) / (batch_size * 10)))
    for (ci in seq_along(chunks)) {
      acc2tax <- ncbi_accession_taxids(chunks[[ci]], batch_size, FALSE)
      taxids <- unique(stats::na.omit(acc2tax$taxid))
      lin <- ncbi_taxid_lineages(taxids, batch_size, FALSE)
      m <- match(acc2tax$taxid, lin$taxid)
      new <- make_lineages(acc2tax$id, acc2tax$taxid,
                           lapply(stats::setNames(taxonomic_ranks(), taxonomic_ranks()),
                                  function(r) lin[[r]][m]))
      cached <- rbind(cached, new)
      cached <- cached[!duplicated(cached$id), , drop = FALSE]
      if (!is.null(cache)) saveRDS(cached[!is.na(cached$taxid), , drop = FALSE], cache)
      msg(verbose, sprintf("    %d / %d", min(ci * batch_size * 10, length(todo)), length(todo)))
    }
  }
  out <- cached[match_ids(ids, cached$id), , drop = FALSE]
  out$id <- ids
  out <- out[!is.na(out$taxid) | rowSums(!is.na(out[, taxonomic_ranks()])) > 0, , drop = FALSE]
  rownames(out) <- NULL
  out
}

need_pkg <- function(pkg, what) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("To %s, please install the '%s' package first: install.packages(\"%s\")",
                 what, pkg, pkg), call. = FALSE)
  }
}

with_retry <- function(expr_fun, tries = 4L) {
  for (i in seq_len(tries)) {
    res <- tryCatch(expr_fun(), error = function(e) e)
    if (!inherits(res, "error")) return(res)
    if (i < tries) Sys.sleep(2^i)
  }
  stop("NCBI request failed after ", tries, " attempts: ", conditionMessage(res), call. = FALSE)
}

ncbi_pause <- function() {
  Sys.sleep(if (nzchar(Sys.getenv("ENTREZ_KEY"))) 0.11 else 0.35)
}

ncbi_accession_taxids <- function(ids, batch_size, verbose) {
  batches <- split(ids, ceiling(seq_along(ids) / batch_size))
  res <- vector("list", length(batches))
  for (b in seq_along(batches)) {
    x <- batches[[b]]
    sm <- tryCatch(with_retry(function() {
      rentrez::entrez_summary(db = "nuccore", id = x, always_return_list = TRUE)
    }), error = function(e) NULL)
    if (is.null(sm)) {
      res[[b]] <- data.frame(id = x, taxid = NA_character_, stringsAsFactors = FALSE)
      next
    }
    acc <- vapply(sm, function(s) if (is.null(s$accessionversion)) NA_character_ else as.character(s$accessionversion), "")
    tx <- vapply(sm, function(s) if (is.null(s$taxid)) NA_character_ else as.character(s$taxid), "")
    m <- match_ids(x, acc)
    res[[b]] <- data.frame(id = x, taxid = unname(tx[m]), stringsAsFactors = FALSE)
    ncbi_pause()
    if (verbose && length(batches) > 1L) message("\r    accessions: batch ", b, "/", length(batches), appendLF = FALSE)
  }
  if (verbose && length(batches) > 1L) message("")
  do.call(rbind, res)
}

ncbi_taxid_lineages <- function(taxids, batch_size, verbose) {
  out <- list()
  if (length(taxids) == 0L) {
    return(make_lineages(character(0), character(0), list()))
  }
  batches <- split(taxids, ceiling(seq_along(taxids) / batch_size))
  for (b in seq_along(batches)) {
    xml <- with_retry(function() {
      rentrez::entrez_fetch(db = "taxonomy", id = batches[[b]], rettype = "xml")
    })
    out[[b]] <- parse_taxonomy_xml(xml)
    ncbi_pause()
    if (verbose && length(batches) > 1L) message("\r    taxa: batch ", b, "/", length(batches), appendLF = FALSE)
  }
  if (verbose && length(batches) > 1L) message("")
  do.call(rbind, out)
}

parse_taxonomy_xml <- function(xml) {
  doc <- xml2::read_xml(xml)
  taxa <- xml2::xml_find_all(doc, "/TaxaSet/Taxon")
  rows <- lapply(taxa, function(t) {
    tid <- xml2::xml_text(xml2::xml_find_first(t, "./TaxId"))
    aka <- xml2::xml_text(xml2::xml_find_all(t, "./AkaTaxIds/TaxId"))
    nm <- xml2::xml_text(xml2::xml_find_first(t, "./ScientificName"))
    rk <- xml2::xml_text(xml2::xml_find_first(t, "./Rank"))
    lin <- xml2::xml_find_all(t, "./LineageEx/Taxon")
    lr <- c(xml2::xml_text(xml2::xml_find_first(lin, "./Rank")), rk)
    ln <- c(xml2::xml_text(xml2::xml_find_first(lin, "./ScientificName")), nm)
    if (!"domain" %in% lr) lr[lr == "superkingdom"] <- "domain"
    ids <- c(tid, aka)
    d <- data.frame(taxid = ids, stringsAsFactors = FALSE)
    for (r in taxonomic_ranks()) {
      j <- which(lr == r)
      d[[r]] <- rep(if (length(j) == 0L) NA_character_ else ln[j[1]], length(ids))
    }
    d
  })
  if (length(rows) == 0L) {
    return(make_lineages(character(0), character(0), list()))
  }
  d <- do.call(rbind, rows)
  d
}

#' Look up taxonomy for accessions in a local taxonomizr database
#'
#' Uses a local SQLite database created by
#' [taxonomizr::prepareDatabase()] (usually called `accessionTaxa.sql`).
#' Building that database downloads tens of gigabytes from NCBI, so for most
#' users [taxonomy_from_ncbi()] or a taxonomy file ([read_taxonomy()]) is
#' easier.
#'
#' @param ids Character vector of nucleotide accessions.
#' @param sqlfile Path to the taxonomizr SQLite database.
#' @return A lineage table (see [read_taxonomy()]).
#' @examples
#' \dontrun{
#' taxonomy_from_taxonomizr(c("MW545606.1"), "accessionTaxa.sql")
#' }
#' @export
taxonomy_from_taxonomizr <- function(ids, sqlfile) {
  need_pkg("taxonomizr", "use a taxonomizr database")
  if (!file.exists(sqlfile)) stop("taxonomizr database not found: ", sqlfile, call. = FALSE)
  ids <- unique(as.character(ids[!is.na(ids) & nzchar(ids)]))
  if (length(ids) == 0L) return(empty_lineages())
  taxid <- taxonomizr::accessionToTaxa(base_accession(ids), sqlfile, version = "base")
  tx <- tryCatch(
    taxonomizr::getTaxonomy(taxid, sqlfile, desiredTaxa = taxonomic_ranks()),
    error = function(e) NULL)
  if (is.null(tx) || all(is.na(tx[, 1]))) {
    # older taxonomizr databases call the first rank "superkingdom"
    tx2 <- tryCatch(taxonomizr::getTaxonomy(
      taxid, sqlfile,
      desiredTaxa = c("superkingdom", taxonomic_ranks()[-1])), error = function(e) NULL)
    if (!is.null(tx2)) tx <- tx2
  }
  if (is.null(tx)) stop("Could not read taxonomy from ", sqlfile, call. = FALSE)
  ranks <- stats::setNames(lapply(seq_len(7L), function(j) tx[, j]), taxonomic_ranks())
  out <- make_lineages(ids, taxid, ranks)
  out[!is.na(out$taxid), , drop = FALSE]
}

# Resolve lineages for ids from one or more taxonomy sources -------------------

resolve_taxonomy_sources <- function(taxonomy) {
  if (is.data.frame(taxonomy)) return(list(taxonomy))
  if (is.character(taxonomy)) taxonomy <- as.list(taxonomy)
  if (!is.list(taxonomy) || length(taxonomy) == 0L) {
    stop(paste0("`taxonomy` must be a taxonomy file, \"ncbi\", a taxonomizr .sql ",
                "database, a lineage table, or a list of these."), call. = FALSE)
  }
  lapply(taxonomy, function(s) {
    if (is.data.frame(s)) {
      if (!all(c("id", taxonomic_ranks()) %in% names(s))) {
        stop("A taxonomy table must have columns: id, ",
             paste(taxonomic_ranks(), collapse = ", "), call. = FALSE)
      }
      return(s)
    }
    if (!is.character(s) || length(s) != 1L) stop("Invalid taxonomy source.", call. = FALSE)
    if (tolower(s) == "ncbi") return("ncbi")
    if (!file.exists(s)) stop("Taxonomy file not found: ", s, call. = FALSE)
    if (grepl("\\.(sql|sqlite|sqlite3|db)$", s, ignore.case = TRUE)) {
      return(structure(s, class = "taxonomizr_db"))
    }
    read_taxonomy(s)
  })
}

get_lineages <- function(ids, sources, verbose = TRUE, ncbi_cache = NULL) {
  ids <- unique(ids[!is.na(ids)])
  out <- empty_lineages()
  todo <- ids
  for (s in sources) {
    if (length(todo) == 0L) break
    if (is.data.frame(s)) {
      m <- match_ids(todo, s$id)
      got <- s[m[!is.na(m)], c("id", "taxid", taxonomic_ranks()), drop = FALSE]
      got$id <- todo[!is.na(m)]
    } else if (inherits(s, "taxonomizr_db")) {
      msg(verbose, "  looking up ", length(todo), " accession(s) in taxonomizr database")
      got <- taxonomy_from_taxonomizr(todo, unclass(s))
    } else if (identical(s, "ncbi")) {
      got <- taxonomy_from_ncbi(todo, cache = ncbi_cache, verbose = verbose)
    }
    out <- rbind(out, got[, names(out), drop = FALSE])
    todo <- setdiff(todo, got$id)
  }
  rownames(out) <- NULL
  out
}
