# Remove duplicate sequences

Keeps the first copy of every distinct sequence, like \`seqkit rmdup
-s\`.

## Usage

``` r
dereplicate(x, both_strands = TRUE)
```

## Arguments

- x:

  A named character vector of sequences, or a path to a FASTA file.

- both_strands:

  If \`TRUE\` (default), a sequence and its reverse complement count as
  duplicates.

## Value

A named character vector with duplicates removed.

## Details

Be careful when de-replicating a \*reference\* ("good") database that is
used by \[clean_reference()\]: if two genera share an identical
sequence, keeping only one copy hides exactly the ambiguity that the
"matches multiple genus" check is designed to detect. De-replicating the
contaminant database is harmless.

## Examples

``` r
seqs <- c(a = "ACGTT", b = "ACGTT", c = "AACGT", d = "GGGCC")
dereplicate(seqs)
#>       a       d 
#> "ACGTT" "GGGCC" 
dereplicate(seqs, both_strands = FALSE)
#>       a       c       d 
#> "ACGTT" "AACGT" "GGGCC" 
```
