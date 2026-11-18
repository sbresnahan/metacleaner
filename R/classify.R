#' Decide which query sequences to flag
#'
#' Applies the metacleaner decision rules to search results. This is the
#' classification step used inside [clean_reference()]; call it directly only
#' if you produced the hit tables yourself (for example with
#' [search_sequences()]).
#'
#' @section Decision rules:
#' Each query is checked in this order; the first rule that applies gives the
#' reason it is flagged.
#' 1. **matches contaminant database**: the query has a qualifying hit (at
#'    `identity` and `coverage`) in the contaminant database.
#' 2. **no match in reference database**: the query has no qualifying hit in
#'    the reference database. If `rescue = TRUE` and the query's best
#'    non-qualifying hit belongs to the same taxon as the query at `level`,
#'    the query is kept and that hit is used in the checks below.
#' 3. **hit outside required <rank>**: (only if `require` is set) at least one
#'    qualifying hit is not in the required taxon, or its taxonomy is unknown.
#' 4. **no taxonomy for query**: the query's own taxonomy at `level` is
#'    unknown.
#' 5. **matches multiple <level>**: the qualifying hits belong to more than one
#'    taxon at `level` (a hit with unknown taxonomy counts as a separate
#'    taxon).
#' 6. **matches different <level>**: all qualifying hits belong to one taxon
#'    that differs from the query's.
#'
#' Queries not flagged by any rule are kept.
#'
#' @param query_ids Character vector of all query IDs.
#' @param reference_hits Data frame of hits against the reference database,
#'   with at least columns `query`, `subject`, `identity`, `coverage`, `score`
#'   and `pass` (as returned by [search_sequences()]).
#' @param contaminant_hits Optional data frame of hits against the contaminant
#'   database (same columns), or `NULL`.
#' @param taxonomy Where to find taxonomy; see [clean_reference()].
#' @param level Rank at which query and hit taxonomy must agree. Default
#'   `"genus"`.
#' @param require Optional named character vector of length one, e.g.
#'   `c(class = "Magnoliopsida")`: every qualifying hit must belong to this
#'   taxon.
#' @param rescue Keep queries without a qualifying reference hit when their
#'   best hit agrees with their taxonomy at `level` (see rule 2). Default
#'   `TRUE`, as in metacleaner 1.x.
#' @param verbose Print progress messages.
#'
#' @return A list with elements `status` (one row per query: `id`, `flagged`,
#'   `reason`, `database`, `top_hit`, `identity`, `coverage`, `n_hits`,
#'   `hit_taxa`) and `lineages` (the lineage table used).
#' @examples
#' hits <- data.frame(query = c("q1", "q2", "q2"), subject = c("r1", "r2", "r3"),
#'                    identity = 100, coverage = 100, score = 1, pass = TRUE)
#' tax <- data.frame(id = c("q1", "q2", "r1", "r2", "r3"), taxid = NA,
#'                   domain = "Eukaryota", phylum = NA, class = NA, order = NA,
#'                   family = NA, genus = c("Malus", "Malus", "Malus", "Malus", "Pyrus"),
#'                   species = NA)
#' res <- classify_queries(c("q1", "q2", "q3"), hits, taxonomy = tax, verbose = FALSE)
#' res$status[, c("id", "flagged", "reason")]
#' @export
classify_queries <- function(query_ids, reference_hits, contaminant_hits = NULL,
                             taxonomy, level = "genus", require = NULL, rescue = TRUE,
                             verbose = TRUE) {
  level <- normalise_rank(level, "level")
  check_flag(rescue, "rescue")
  req <- check_require(require)
  sources <- resolve_taxonomy_sources(taxonomy)
  ids <- unique(as.character(query_ids))
  n <- length(ids)
  st <- data.frame(id = ids, flagged = FALSE, reason = NA_character_,
                   database = NA_character_, top_hit = NA_character_,
                   identity = NA_real_, coverage = NA_real_, n_hits = 0L,
                   hit_taxa = NA_character_, stringsAsFactors = FALSE)
  rh <- order_hits(reference_hits)
  ch <- if (is.null(contaminant_hits)) NULL else order_hits(contaminant_hits)

  set_flag <- function(st, qid, reason, database, hits = NULL) {
    i <- match(qid, st$id)
    st$flagged[i] <- TRUE
    st$reason[i] <- reason
    st$database[i] <- database
    if (!is.null(hits)) {
      j <- match(qid, hits$query)
      st$top_hit[i] <- hits$subject[j]
      st$identity[i] <- hits$identity[j]
      st$coverage[i] <- hits$coverage[j]
    }
    st
  }

  # 1. contaminants
  if (!is.null(ch)) {
    cp <- ch[ch$pass & ch$query %in% ids, , drop = FALSE]
    bad <- unique(cp$query)
    if (length(bad)) {
      st <- set_flag(st, bad, "matches contaminant database", "contaminants", cp)
      st$n_hits[match(bad, st$id)] <- as.integer(table(cp$query)[bad])
    }
    msg(verbose, "  ", length(bad), " sequence(s) match the contaminant database")
  }
  remaining <- st$id[!st$flagged]

  # 2. reference database
  rh <- rh[rh$query %in% remaining, , drop = FALSE]
  rp <- rh[rh$pass, , drop = FALSE]
  nomatch <- setdiff(remaining, rp$query)
  rf <- rh[!rh$pass & rh$query %in% nomatch, , drop = FALSE]
  rf <- rf[!duplicated(rf$query), , drop = FALSE]

  # taxonomy for everything we need
  need <- unique(c(remaining, rp$subject, rf$subject,
                   if (!is.null(ch)) st$top_hit[!is.na(st$top_hit)]))
  msg(verbose, "  retrieving taxonomy for ", length(need), " sequence(s)")
  lin <- get_lineages(need, sources, verbose = verbose)
  tax_at <- function(x, rank = level) lin[[rank]][match(x, lin$id)]
  missing_q <- sum(is.na(match(remaining, lin$id)))
  if (missing_q > 0L) {
    msg(verbose, "  note: no taxonomy found for ", missing_q, " query sequence(s)")
  }

  rescued <- character(0)
  if (rescue && nrow(rf)) {
    qt <- tax_at(rf$query)
    stx <- tax_at(rf$subject)
    ok <- !is.na(qt) & !is.na(stx) & qt == stx
    rescued <- rf$query[ok]
    if (length(rescued)) rp <- rbind(rp, rf[ok, , drop = FALSE])
  }
  nm <- setdiff(nomatch, rescued)
  if (length(nm)) st <- set_flag(st, nm, "no match in reference database", "reference", rf)
  msg(verbose, "  ", length(nm), " sequence(s) have no qualifying match in the reference database",
      if (rescue) paste0(" (", length(rescued), " more rescued by their best hit)") else "")

  rp <- order_hits(rp)
  cand <- setdiff(remaining, nm)
  cnt <- table(rp$query)
  st$n_hits[match(cand, st$id)] <- as.integer(cnt[cand])
  st$n_hits[is.na(st$n_hits)] <- 0L
  hit_tax <- tax_at(rp$subject)
  hit_tax_lab <- ifelse(is.na(hit_tax), "<unknown>", hit_tax)
  taxa_by_q <- split(hit_tax_lab, factor(rp$query, levels = cand))
  st$hit_taxa[match(cand, st$id)] <- vapply(taxa_by_q, function(v) {
    paste(unique(v), collapse = "; ")
  }, "")

  # 3. required taxon
  if (!is.null(req)) {
    rt <- tax_at(rp$subject, req$rank)
    off <- is.na(rt) | rt != req$value
    offq <- unique(rp$query[off])
    if (length(offq)) {
      st <- set_flag(st, offq, sprintf("hit outside required %s (%s)", req$rank, req$value),
                     "reference", rp[off, , drop = FALSE])
    }
    msg(verbose, "  ", length(offq), " sequence(s) have hits outside ", req$rank, " ", req$value)
    cand <- setdiff(cand, offq)
  }

  # 4. query taxonomy
  qtax <- tax_at(cand)
  notax <- cand[is.na(qtax)]
  if (length(notax)) {
    st <- set_flag(st, notax, sprintf("no taxonomy for query at %s", level), "reference",
                   rp[rp$query %in% notax, , drop = FALSE])
  }
  cand <- cand[!is.na(qtax)]
  qtax <- qtax[!is.na(qtax)]

  # 5./6. agreement at level
  if (length(cand)) {
    sub <- rp[rp$query %in% cand, , drop = FALSE]
    sub$taxon <- ifelse(is.na(tax_at(sub$subject)), "<unknown>", tax_at(sub$subject))
    sub$qtaxon <- qtax[match(sub$query, cand)]
    ntax <- tapply(sub$taxon, factor(sub$query, levels = cand), function(v) length(unique(v)))
    multi <- cand[!is.na(ntax) & ntax > 1L]
    if (length(multi)) {
      off <- sub[sub$query %in% multi & sub$taxon != sub$qtaxon, , drop = FALSE]
      st <- set_flag(st, multi, sprintf("matches multiple %s", level), "reference", off)
    }
    single <- cand[!is.na(ntax) & ntax == 1L]
    first <- sub[sub$query %in% single & !duplicated(sub$query), , drop = FALSE]
    wrong <- first$query[first$taxon != first$qtaxon]
    if (length(wrong)) {
      st <- set_flag(st, wrong, sprintf("matches different %s", level), "reference", first)
    }
    msg(verbose, "  ", length(multi), " sequence(s) match more than one ", level, "; ",
        length(wrong), " match a different ", level)
  }

  # report best hit for kept sequences too
  kept <- st$id[!st$flagged]
  j <- match(kept, rp$query)
  i <- match(kept, st$id)
  st$top_hit[i] <- rp$subject[j]
  st$identity[i] <- rp$identity[j]
  st$coverage[i] <- rp$coverage[j]

  list(status = st, lineages = lin)
}

order_hits <- function(h) {
  if (is.null(h) || nrow(h) == 0L) {
    h <- empty_hits()
    return(h)
  }
  need <- c("query", "subject", "identity", "coverage", "score", "pass")
  miss <- setdiff(need, names(h))
  if (length(miss)) {
    stop("Hit table is missing column(s): ", paste(miss, collapse = ", "), call. = FALSE)
  }
  h$query <- as.character(h$query)
  h$subject <- as.character(h$subject)
  h[order(h$query, !h$pass, -h$score, -h$identity, -h$coverage, h$subject), , drop = FALSE]
}

check_require <- function(require) {
  if (is.null(require)) return(NULL)
  if (is.character(require) && length(require) == 1L && is.null(names(require)) &&
      grepl(",", require, fixed = TRUE)) {
    # old command-line style "class,Magnoliopsida"
    p <- trimws(strsplit(require, ",", fixed = TRUE)[[1]])
    require <- stats::setNames(p[2], p[1])
  }
  if (!is.character(require) || length(require) != 1L || is.null(names(require))) {
    stop("`require` must look like c(class = \"Magnoliopsida\").", call. = FALSE)
  }
  list(rank = normalise_rank(names(require), "require"), value = unname(require))
}
