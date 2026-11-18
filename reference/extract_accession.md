# Extract accessions from FASTA header lines

Returns the first word of each header. Old-style NCBI headers
(\`gi\|12345\|gb\|AB123456.1\|\`) and pipe-delimited headers
(\`gb\|AB123456.1\|\`) are reduced to the accession.

## Usage

``` r
extract_accession(header)
```

## Arguments

- header:

  Character vector of FASTA header lines (with or without the leading
  \`\>\`).

## Value

A character vector of accessions.

## Examples

``` r
extract_accession(c(">MN123456.1 Malus domestica ITS2",
                    "gi|12345|gb|AB123456.1|"))
#> [1] "MN123456.1" "AB123456.1"
```
