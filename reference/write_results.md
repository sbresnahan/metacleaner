# Save metacleaner results to files

Writes the cleaned database and a report of the removed sequences:

## Usage

``` r
write_results(
  result,
  outdir,
  prefix = NULL,
  tax_format = c("taxonomizr", "lineage"),
  overwrite = TRUE,
  verbose = TRUE
)
```

## Arguments

- result:

  The output of \[clean_reference()\].

- outdir:

  Folder to write to (created if needed), e.g. \`"cleaned"\`.

- prefix:

  File name prefix. Default: the prefix stored in \`result\`.

- tax_format:

  Format of the \`\_clean.tax\` file: \`"taxonomizr"\` (default; the
  comma-separated format used by metacleaner 1.x and as MetaCurator
  input: \`taxid,accession,taxid,domain,...,species\`) or \`"lineage"\`
  (tab-separated \`accession\` and a \`k\_\_;p\_\_;...;s\_\_\` lineage
  string, as used by QIIME-style classifiers).

- overwrite:

  Overwrite existing files? Default \`TRUE\`.

- verbose:

  Print the names of the written files. Default \`TRUE\`.

## Value

The paths of the written files, invisibly.

## Details

\* \`\<prefix\>\_clean.fasta\`: the kept sequences (one line per
sequence); \* \`\<prefix\>\_clean.tax\`: taxonomy of the kept sequences;
\* \`\<prefix\>\_flagged.tsv\`: one row per removed sequence with its
taxonomy, the hit that caused it to be flagged, that hit's taxonomy and
the reason; \* \`\<prefix\>\_summary.txt\`: settings and counts per
reason.

## Examples

``` r
ex <- function(f) system.file("extdata", f, package = "metacleaner")
res <- clean_reference(ex("query.fasta"), ex("reference.fasta"),
                       ex("contaminants.fasta"), taxonomy = ex("taxonomy.tax"),
                       verbose = FALSE)
files <- write_results(res, outdir = tempdir(), prefix = "example")
#> Results written to /tmp/RtmpX8yzUw:
#>   example_clean.fasta
#>   example_clean.tax
#>   example_flagged.tsv
#>   example_summary.txt
files
#>                                 fasta                              taxonomy 
#> "/tmp/RtmpX8yzUw/example_clean.fasta"   "/tmp/RtmpX8yzUw/example_clean.tax" 
#>                               flagged                               summary 
#> "/tmp/RtmpX8yzUw/example_flagged.tsv" "/tmp/RtmpX8yzUw/example_summary.txt" 
```
