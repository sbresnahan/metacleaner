# Download sequences from NCBI GenBank

Downloads all NCBI Nucleotide records that match a search term (written
exactly as you would type it on the NCBI website) to a FASTA file, in
batches. Requires an internet connection and the suggested package
rentrez.

## Usage

``` r
download_ncbi_sequences(
  term,
  file,
  batch_size = 500,
  max_records = Inf,
  resume = FALSE,
  api_key = Sys.getenv("ENTREZ_KEY"),
  verbose = TRUE
)
```

## Arguments

- term:

  NCBI Nucleotide search term, e.g. \`'(ITS2) AND
  Embryophyta\[Organism\] AND ("0"\[SLEN\] : "10000"\[SLEN\])'\`.

- file:

  Output FASTA file.

- batch_size:

  Records per request (max 10000). Default 500.

- max_records:

  Maximum number of records to download. Default: all.

- resume:

  If \`TRUE\` and \`file\` exists, continue an interrupted download by
  appending to it. If \`FALSE\` (default) an existing file is an error.

- api_key:

  Optional NCBI API key. Defaults to the \`ENTREZ_KEY\` environment
  variable.

- verbose:

  Print progress. Default \`TRUE\`.

## Value

The number of records in \`file\`, invisibly.

## Details

Large downloads can take a long time. If a download is interrupted, run
the same command again with \`resume = TRUE\` to continue where it
stopped. An NCBI API key (free, see
\<https://account.ncbi.nlm.nih.gov/\>) makes downloads faster.

## See also

\[build_reference_db()\]

## Examples

``` r
# \donttest{
if (requireNamespace("rentrez", quietly = TRUE)) {
  out <- tempfile(fileext = ".fasta")
  try(download_ncbi_sequences("Malus[Organism] AND ITS2", out, max_records = 20))
}
#> NCBI search found 399 record(s); downloading 20
#>   20 / 20
# }
```
