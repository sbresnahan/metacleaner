# Read sequences from a FASTA file

Reads a FASTA file (plain text or gzip-compressed, with sequences on one
line or wrapped over many lines) into a named character vector. By
default each sequence is named by its accession, i.e. the first word of
its header line, with old-style NCBI headers such as
\`gi\|12345\|gb\|AB123456.1\|\` reduced to \`AB123456.1\`.

## Usage

``` r
read_fasta(file, ids = c("accession", "header"), drop_version = FALSE)
```

## Arguments

- file:

  Path to a FASTA file. Files ending in \`.gz\` are decompressed
  automatically.

- ids:

  How sequences are named: \`"accession"\` (default) uses the accession
  (first word of the header); \`"header"\` keeps the full header line.

- drop_version:

  If \`TRUE\`, remove the version suffix from accessions (\`AB123456.1\`
  becomes \`AB123456\`). Default \`FALSE\`.

## Value

A named character vector: the names are sequence IDs and the values are
the sequences in upper case.

## See also

\[write_fasta()\], \[dereplicate()\]

## Examples

``` r
f <- system.file("extdata", "query.fasta", package = "metacleaner")
seqs <- read_fasta(f)
length(seqs)
#> [1] 52
head(names(seqs))
#> [1] "MZ509216.1" "MZ509176.1" "MW545606.1" "MW545649.1" "MW545773.1"
#> [6] "MW545611.1"
```
