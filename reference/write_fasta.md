# Write sequences to a FASTA file

Write sequences to a FASTA file

## Usage

``` r
write_fasta(x, file, width = 0)
```

## Arguments

- x:

  A named character vector of sequences (for example from
  \[read_fasta()\]). The names are written as header lines.

- file:

  Output path. If it ends in \`.gz\` the file is gzip-compressed.

- width:

  Maximum number of bases per line. The default, \`0\`, writes each
  sequence on a single line, which is what most metabarcoding tools
  expect.

## Value

The output path, invisibly.

## Examples

``` r
seqs <- c(seq1 = "ACGTACGT", seq2 = "TTGACCA")
out <- tempfile(fileext = ".fasta")
write_fasta(seqs, out)
readLines(out)
#> [1] ">seq1"    "ACGTACGT" ">seq2"    "TTGACCA" 
```
