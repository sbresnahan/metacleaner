# Taxonomic ranks used by metacleaner

The seven ranks that metacleaner reports and that can be used for
\`level\` and \`require\`. \`"superkingdom"\` and \`"kingdom"\` are
accepted as synonyms of \`"domain"\` (NCBI renamed this rank in 2025).

## Usage

``` r
taxonomic_ranks()
```

## Value

A character vector of rank names, from broadest to narrowest.

## Examples

``` r
taxonomic_ranks()
#> [1] "domain"  "phylum"  "class"   "order"   "family"  "genus"   "species"
```
