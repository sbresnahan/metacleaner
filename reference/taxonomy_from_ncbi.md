# Look up taxonomy for accessions from NCBI

Looks up the lineage of each accession with the NCBI E-utilities web
service (no local database needed). Only the accessions you ask for are
downloaded, in batches. Requires an internet connection and the
suggested packages rentrez and xml2.

## Usage

``` r
taxonomy_from_ncbi(
  ids,
  api_key = Sys.getenv("ENTREZ_KEY"),
  cache = NULL,
  batch_size = 200,
  verbose = TRUE
)
```

## Arguments

- ids:

  Character vector of nucleotide accessions.

- api_key:

  Optional NCBI API key. Defaults to the \`ENTREZ_KEY\` environment
  variable if set.

- cache:

  Optional path to an \`.rds\` file used to store and reuse lookups
  between runs.

- batch_size:

  Number of IDs per web request. Default 200.

- verbose:

  Print progress messages. Default \`TRUE\`.

## Value

A lineage table (see \[read_taxonomy()\]) with one row per ID that could
be resolved.

## Details

NCBI limits anonymous users to 3 requests per second; with a free API
key (see \<https://account.ncbi.nlm.nih.gov/\>) the limit is 10 per
second. Use \`cache\` to store results in a file so that repeated runs
do not download the same information again.

## Examples

``` r
# \donttest{
if (requireNamespace("rentrez", quietly = TRUE) &&
    requireNamespace("xml2", quietly = TRUE)) {
  try(taxonomy_from_ncbi(c("MW545606.1", "MZ509176.1")))
}
#>   looking up 2 accession(s) at NCBI
#>     2 / 2
#>           id   taxid    domain       phylum         class     order     family
#> 1 MW545606.1  463586 Eukaryota Streptophyta Magnoliopsida Asterales Asteraceae
#> 2 MZ509176.1 2909194 Eukaryota Streptophyta Magnoliopsida   Fabales   Fabaceae
#>      genus           species
#> 1   Ajania  Ajania potaninii
#> 2 Pueraria Pueraria sp. KK8s
# }
```
