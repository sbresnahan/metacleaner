# Reverse complement DNA sequences

Reverse complement DNA sequences

## Usage

``` r
reverse_complement(x)
```

## Arguments

- x:

  Character vector of DNA sequences (IUPAC ambiguity codes are
  complemented correctly).

## Value

A character vector of the same length, with names kept.

## Examples

``` r
reverse_complement(c(a = "AACGTN", b = "RYKM"))
#>        a        b 
#> "NACGTT"   "KMRY" 
```
