# Check whether BLAST+ is installed

metacleaner does not need BLAST+, but it can use it as an alternative
search engine (\`engine = "blast"\`). This function reports whether the
\`makeblastdb\` and \`blastn\` programs can be found on the system
\`PATH\`.

## Usage

``` r
blast_available(verbose = FALSE)
```

## Arguments

- verbose:

  If \`TRUE\`, print where the programs were found or how to install
  them.

## Value

\`TRUE\` or \`FALSE\`, invisibly when \`verbose = TRUE\`.

## Examples

``` r
blast_available()
#> [1] FALSE
```
