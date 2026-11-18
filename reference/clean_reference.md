# Clean a barcode reference database

The main metacleaner function. It removes sequences that are probably
mislabeled or contaminants from a DNA barcode reference database (for
example a plant ITS2 database built from NCBI GenBank records). For each
query sequence it:

## Usage

``` r
clean_reference(
  query,
  reference,
  contaminants = NULL,
  taxonomy,
  identity = 100,
  coverage = 100,
  level = "genus",
  require = NULL,
  rescue = TRUE,
  exclude_self = FALSE,
  engine = c("native", "blast"),
  threads = 1,
  max_hits = 100,
  outdir = NULL,
  prefix = NULL,
  verbose = TRUE
)
```

## Arguments

- query:

  The database to clean: a FASTA file path or a named character vector
  (see \[read_fasta()\]). Sequence names must be accessions that can be
  found in \`taxonomy\`.

- reference:

  Sequences of the wanted group (e.g. all plant ITS sequences from
  NCBI): FASTA file path or named character vector. It may contain the
  query sequences themselves.

- contaminants:

  Optional sequences of unwanted groups (e.g. fungal ITS): FASTA file
  path or named character vector. \`NULL\` (default) skips this step.

- taxonomy:

  Where to find the taxonomy of query and reference sequences. One of:
  \* a taxonomy file path (see \[read_taxonomy()\] for the formats), \*
  \`"ncbi"\` to look up accessions online (see
  \[taxonomy_from_ncbi()\]), \* the path to a taxonomizr database ending
  in \`.sql\` (see \[taxonomy_from_taxonomizr()\]), \* a lineage table
  (data frame) such as returned by \[read_taxonomy()\], \* or a list /
  vector of several of these, tried in order: e.g. \`c("query.tax",
  "ncbi")\` uses the file first and NCBI for anything not in the file.

- identity:

  Minimum percent identity for a hit (0-100). Default 100.

- coverage:

  Minimum percent of the query covered by a hit (0-100). Default 100.

- level:

  Taxonomic rank at which a query and its hits must agree: one of
  \`"domain"\`, \`"phylum"\`, \`"class"\`, \`"order"\`, \`"family"\`,
  \`"genus"\` (default) or \`"species"\`.

- require:

  Optional extra filter, e.g. \`c(class = "Magnoliopsida")\`: queries
  with any hit outside this taxon are flagged. (The old command-line
  form \`"class,Magnoliopsida"\` is also accepted.)

- rescue:

  If \`TRUE\` (default, as in metacleaner 1.x), a query with no hit at
  the thresholds is kept when its best lower-scoring hit belongs to the
  same taxon at \`level\`. Set to \`FALSE\` for stricter cleaning.

- exclude_self:

  If \`TRUE\`, a query may not match the reference sequence with the
  same accession. Default \`FALSE\`.

- engine:

  Search engine: \`"native"\` (default, built in) or \`"blast"\`
  (requires BLAST+; see \[blast_available()\]).

- threads:

  Number of CPU threads to use, or \`"auto"\` to use all available cores
  but one. Default 1. The search is split across threads on all
  operating systems (including macOS); results are identical whatever
  the number of threads.

- max_hits:

  Maximum hits per query used in the taxonomy checks (best-scoring
  first). Default 100.

- outdir:

  Optional output folder. If given, the results are also written there
  with \[write_results()\].

- prefix:

  File name prefix for the output files. Defaults to the query file name
  without extension, or \`"metacleaner"\`.

- verbose:

  Print progress messages. Default \`TRUE\`.

## Value

An object of class \`metacleaner_result\`: a list with \* \`clean\`: the
kept sequences (named character vector), \* \`flagged\`: a data frame
with one row per removed sequence and the reason, the hit that triggered
it and both taxonomies, \* \`status\`: one row per query sequence (kept
or not), \* \`taxonomy\`: lineage table of the kept sequences, \*
\`hits\`: the search results (\`reference\` and \`contaminants\`), \*
\`settings\`: the settings used.

Print it for a summary; use \[write_results()\] to save it.

## Details

1\. searches a database of unwanted sequences (\`contaminants\`, e.g.
fungal ITS) and flags queries that match one; 2. searches a database of
wanted sequences (\`reference\`, e.g. plant ITS) and flags queries
without a match; 3. compares the taxonomy of each query with the
taxonomy of its matches at the chosen \`level\` (genus by default) and
flags queries whose matches belong to a different or to several taxa.

A match ("hit") counts only if it reaches both \`identity\` and
\`coverage\`. The full decision rules are listed in
\[classify_queries()\].

## See also

\[write_results()\], \[read_fasta()\], \[read_taxonomy()\],
\[build_reference_db()\]

## Examples

``` r
ex <- function(f) system.file("extdata", f, package = "metacleaner")
res <- clean_reference(
  query        = ex("query.fasta"),
  reference    = ex("reference.fasta"),
  contaminants = ex("contaminants.fasta"),
  taxonomy     = ex("taxonomy.tax"),
  level        = "genus"
)
#> metacleaner: reading sequences
#>   52 query, 168 reference, 33 contaminant sequences
#> Step 1/3: searching the contaminant database
#> Step 2/3: searching the reference database
#> Step 3/3: comparing taxonomy at genus level
#>   3 sequence(s) match the contaminant database
#>   retrieving taxonomy for 171 sequence(s)
#>   0 sequence(s) have no qualifying match in the reference database (0 more rescued by their best hit)
#>   6 sequence(s) match more than one genus; 0 match a different genus
#> Done: kept 40 of 52 sequences, flagged 12 (0.3 s).
res
#> <metacleaner result>
#>   query sequences:   52
#>   kept:              40
#>   flagged:           12
#>     matches multiple genus                   6
#>     matches contaminant database             3
#>     no taxonomy for query at genus           3
#>   settings: identity >= 100%, coverage >= 100%, level = genus, engine = native
#> Use write_results() to save the cleaned database and report.
head(res$flagged)
#>           id
#> 1 KF984840.1
#> 2 MN228964.1
#> 3 MW545606.1
#> 4 MW545611.1
#> 5 MW545649.1
#> 6 MW545773.1
#>                                                                                    query_taxonomy
#> 1 Eukaryota|Ascomycota|Eurotiomycetes|Eurotiales|Trichocomaceae|Talaromyces|Talaromyces_rugulosus
#> 2             Eukaryota|Streptophyta|Magnoliopsida|Asterales|Asteraceae|NA|Dicomeae_sp._Funk13443
#> 3               Eukaryota|Streptophyta|Magnoliopsida|Asterales|Asteraceae|Ajania|Ajania_potaninii
#> 4           Eukaryota|Streptophyta|Magnoliopsida|Asterales|Asteraceae|Ajania|Ajania_scharnhorstii
#> 5 Eukaryota|Streptophyta|Magnoliopsida|Asterales|Asteraceae|Chrysanthemum|Chrysanthemum_zawadzkii
#> 6         Eukaryota|Streptophyta|Magnoliopsida|Asterales|Asteraceae|Artemisia|Artemisia_pubescens
#>      top_hit
#> 1 KF984840.1
#> 2 MN228964.1
#> 3 MW545600.1
#> 4 MW545600.1
#> 5 MW545607.1
#> 6 MW545787.1
#>                                                                                        hit_taxonomy
#> 1   Eukaryota|Ascomycota|Eurotiomycetes|Eurotiales|Trichocomaceae|Talaromyces|Talaromyces_rugulosus
#> 2               Eukaryota|Streptophyta|Magnoliopsida|Asterales|Asteraceae|NA|Dicomeae_sp._Funk13443
#> 3 Eukaryota|Streptophyta|Magnoliopsida|Asterales|Asteraceae|Chrysanthemum|Chrysanthemum_przewalskii
#> 4 Eukaryota|Streptophyta|Magnoliopsida|Asterales|Asteraceae|Chrysanthemum|Chrysanthemum_przewalskii
#> 5           Eukaryota|Streptophyta|Magnoliopsida|Asterales|Asteraceae|Ajania|Ajania_tripinnatisecta
#> 6         Eukaryota|Streptophyta|Magnoliopsida|Asterales|Asteraceae|Filifolium|Filifolium_sibiricum
#>                           reason     database identity coverage n_hits
#> 1   matches contaminant database contaminants      100      100      1
#> 2 no taxonomy for query at genus    reference      100      100      1
#> 3         matches multiple genus    reference      100      100      8
#> 4         matches multiple genus    reference      100      100      8
#> 5         matches multiple genus    reference      100      100     11
#> 6         matches multiple genus    reference      100      100      4
#>                hit_taxa
#> 1                  <NA>
#> 2             <unknown>
#> 3 Chrysanthemum; Ajania
#> 4 Chrysanthemum; Ajania
#> 5 Ajania; Chrysanthemum
#> 6 Artemisia; Filifolium
# \donttest{
# save the cleaned database and the report
write_results(res, outdir = tempdir(), prefix = "example")
#> Results written to /tmp/RtmpX8yzUw:
#>   example_clean.fasta
#>   example_clean.tax
#>   example_flagged.tsv
#>   example_summary.txt
# }
```
