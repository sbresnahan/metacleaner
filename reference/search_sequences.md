# Search query sequences against a database of sequences

Finds, for every query sequence, the database sequences that it matches
with at least \`identity\` percent identity over at least \`coverage\`
percent of the query length (both strands are searched). This is the
search step used inside \[clean_reference()\]; most users do not need to
call it directly.

## Usage

``` r
search_sequences(
  query,
  database,
  identity = 100,
  coverage = 100,
  max_hits = 100,
  best_hit = TRUE,
  exclude_self = FALSE,
  engine = c("native", "blast"),
  threads = 1,
  k = 12,
  verbose = TRUE
)
```

## Arguments

- query:

  Query sequences: a FASTA file path or a named character vector.

- database:

  Database sequences: a FASTA file path or a named character vector.

- identity:

  Minimum percent identity (0-100). Default 100.

- coverage:

  Minimum percent query coverage (0-100). Default 100.

- max_hits:

  Maximum number of qualifying hits reported per query (best first).
  Default 100.

- best_hit:

  If \`TRUE\` (default), queries without any qualifying hit also report
  their best-scoring non-qualifying hit (with \`pass = FALSE\`), when
  one exists.

- exclude_self:

  If \`TRUE\`, a query is not allowed to match a database sequence with
  the same ID. Default \`FALSE\`.

- engine:

  \`"native"\` (default) or \`"blast"\`.

- threads:

  Number of CPU threads to use, or \`"auto"\` to use all available cores
  but one. Default 1. The search is split across threads on all
  operating systems (including macOS); results are identical whatever
  the number of threads.

- k:

  Word length for the native engine index (8-13). Default 12. Smaller
  values are more sensitive at low identity but slower.

- verbose:

  Print progress messages. Default \`TRUE\`.

## Value

A data frame with one row per reported hit and columns \`query\`,
\`subject\` (database ID), \`identity\`, \`coverage\`, \`score\`
(alignment score for the native engine, bit score for BLAST),
\`strand\`, \`qstart\`, \`qend\`, \`sstart\`, \`send\`, \`pass\`
(whether the hit meets both thresholds) and \`n_pass\` (total number of
qualifying hits for that query, which can exceed \`max_hits\`).

## Native engine

The default engine is written in C++ and needs no external software. The
database is indexed by k-mers (words of length \`k\`). For each query a
q-gram filter computes how many k-mers any qualifying match must share
with the query; only database sequences that pass this filter are
aligned (banded local alignment with megablast-like scores: match +1,
mismatch -2, gap 2.5). Because the filter is derived from the
thresholds, no qualifying match is missed when the thresholds are high
(identity of about 90 search is an exact containment test and is very
fast.

Identity is the number of identical positions divided by the alignment
length (as BLAST's \`pident\`). Upper/lower case and U/T are treated as
the same; an ambiguity code (e.g. \`N\`, \`R\`, \`Y\`) is identical only
to the same code. BLAST instead scores ambiguity codes as mismatches, so
a sequence containing them can be 100 engine but not with BLAST.
Coverage is the fraction of the query covered by the best alignment to
that database sequence. BLAST's \`qcovs\` adds together several separate
alignments to the same sequence, so coverage values close to the
threshold can differ slightly between the engines.

## BLAST engine

With \`engine = "blast"\`, BLAST+ (\`makeblastdb\` and \`blastn\`,
default megablast settings) must be installed and on the system
\`PATH\`; see \[blast_available()\]. Identity is \`pident\` of the
best-scoring HSP and coverage is \`qcovs\`, as in metacleaner 1.x.

## See also

\[clean_reference()\]

## Examples

``` r
q <- c(q1 = "ACGTACGTTAGCATGCATCGATCGATGCTAGCTAGG")
db <- c(s1 = "TTTTACGTACGTTAGCATGCATCGATCGATGCTAGCTAGGTTTT",
        s2 = "ACGTACGTTAGCATGCATCGTTCGATGCTAGCTAGG")
search_sequences(q, db, verbose = FALSE)
#>   query subject identity coverage score strand qstart qend sstart send pass
#> 1    q1      s1      100      100    72      +      1   36      5   40 TRUE
#>   n_pass
#> 1      1
search_sequences(q, db, identity = 95, coverage = 95, verbose = FALSE)
#>   query subject  identity coverage score strand qstart qend sstart send pass
#> 1    q1      s1 100.00000      100    72      +      1   36      5   40 TRUE
#> 2    q1      s2  97.22222      100    66      +      1   36      1   36 TRUE
#>   n_pass
#> 1      2
#> 2      2
```
