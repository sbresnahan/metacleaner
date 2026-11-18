# Look up taxonomy for accessions in a local taxonomizr database

Uses a local SQLite database created by
\[taxonomizr::prepareDatabase()\] (usually called
\`accessionTaxa.sql\`). Building that database downloads tens of
gigabytes from NCBI, so for most users \[taxonomy_from_ncbi()\] or a
taxonomy file (\[read_taxonomy()\]) is easier.

## Usage

``` r
taxonomy_from_taxonomizr(ids, sqlfile)
```

## Arguments

- ids:

  Character vector of nucleotide accessions.

- sqlfile:

  Path to the taxonomizr SQLite database.

## Value

A lineage table (see \[read_taxonomy()\]).

## Examples

``` r
if (FALSE) { # \dontrun{
taxonomy_from_taxonomizr(c("MW545606.1"), "accessionTaxa.sql")
} # }
```
