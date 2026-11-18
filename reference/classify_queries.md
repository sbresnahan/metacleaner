# Decide which query sequences to flag

Applies the metacleaner decision rules to search results. This is the
classification step used inside \[clean_reference()\]; call it directly
only if you produced the hit tables yourself (for example with
\[search_sequences()\]).

## Usage

``` r
classify_queries(
  query_ids,
  reference_hits,
  contaminant_hits = NULL,
  taxonomy,
  level = "genus",
  require = NULL,
  rescue = TRUE,
  verbose = TRUE
)
```

## Arguments

- query_ids:

  Character vector of all query IDs.

- reference_hits:

  Data frame of hits against the reference database, with at least
  columns \`query\`, \`subject\`, \`identity\`, \`coverage\`, \`score\`
  and \`pass\` (as returned by \[search_sequences()\]).

- contaminant_hits:

  Optional data frame of hits against the contaminant database (same
  columns), or \`NULL\`.

- taxonomy:

  Where to find taxonomy; see \[clean_reference()\].

- level:

  Rank at which query and hit taxonomy must agree. Default \`"genus"\`.

- require:

  Optional named character vector of length one, e.g. \`c(class =
  "Magnoliopsida")\`: every qualifying hit must belong to this taxon.

- rescue:

  Keep queries without a qualifying reference hit when their best hit
  agrees with their taxonomy at \`level\` (see rule 2). Default
  \`TRUE\`, as in metacleaner 1.x.

- verbose:

  Print progress messages.

## Value

A list with elements \`status\` (one row per query: \`id\`, \`flagged\`,
\`reason\`, \`database\`, \`top_hit\`, \`identity\`, \`coverage\`,
\`n_hits\`, \`hit_taxa\`) and \`lineages\` (the lineage table used).

## Decision rules

Each query is checked in this order; the first rule that applies gives
the reason it is flagged. 1. \*\*matches contaminant database\*\*: the
query has a qualifying hit (at \`identity\` and \`coverage\`) in the
contaminant database. 2. \*\*no match in reference database\*\*: the
query has no qualifying hit in the reference database. If \`rescue =
TRUE\` and the query's best non-qualifying hit belongs to the same taxon
as the query at \`level\`, the query is kept and that hit is used in the
checks below. 3. \*\*hit outside required \<rank\>\*\*: (only if
\`require\` is set) at least one qualifying hit is not in the required
taxon, or its taxonomy is unknown. 4. \*\*no taxonomy for query\*\*: the
query's own taxonomy at \`level\` is unknown. 5. \*\*matches multiple
\<level\>\*\*: the qualifying hits belong to more than one taxon at
\`level\` (a hit with unknown taxonomy counts as a separate taxon). 6.
\*\*matches different \<level\>\*\*: all qualifying hits belong to one
taxon that differs from the query's.

Queries not flagged by any rule are kept.

## Examples

``` r
hits <- data.frame(query = c("q1", "q2", "q2"), subject = c("r1", "r2", "r3"),
                   identity = 100, coverage = 100, score = 1, pass = TRUE)
tax <- data.frame(id = c("q1", "q2", "r1", "r2", "r3"), taxid = NA,
                  domain = "Eukaryota", phylum = NA, class = NA, order = NA,
                  family = NA, genus = c("Malus", "Malus", "Malus", "Malus", "Pyrus"),
                  species = NA)
res <- classify_queries(c("q1", "q2", "q3"), hits, taxonomy = tax, verbose = FALSE)
res$status[, c("id", "flagged", "reason")]
#>   id flagged                         reason
#> 1 q1   FALSE                           <NA>
#> 2 q2    TRUE         matches multiple genus
#> 3 q3    TRUE no match in reference database
```
