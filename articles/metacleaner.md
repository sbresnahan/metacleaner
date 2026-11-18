# Getting started with metacleaner

## What metacleaner does

Reference databases for DNA metabarcoding (for example plant ITS2
databases built from NCBI GenBank) contain sequences with the wrong
name: fungal sequences submitted as plants, sequences labeled with the
wrong genus, and sequences that are identical in several genera and so
cannot identify any of them. Reads that match these entries are assigned
to the wrong taxon.

metacleaner checks every sequence in your database (the **query**) in
three steps:

1.  Does it match a sequence in a database of things that should *not*
    be there (the **contaminants**, e.g. fungal ITS)? Then it is
    removed.
2.  Does it match a sequence in a database of things that *should* be
    there (the **reference**, e.g. all plant ITS)? If not, it is
    removed.
3.  Do the reference sequences it matches have the same taxonomy as the
    query at the level you choose (genus by default)? If they belong to
    a different genus, or to several genera, it is removed.

A “match” means at least `identity` percent identical over at least
`coverage` percent of the query’s length. Both are 100 by default, so
only identical sequences count.

You need no programming experience and no other software: everything,
including the sequence search, is built into the package.

## Installing

Install from CRAN:

``` r

install.packages("metacleaner", dependencies = TRUE)
```

`dependencies = TRUE` also installs the optional packages used to
download data and taxonomy from NCBI.

If R reports that package ‘metacleaner’ is not available, it is still
under review at CRAN. Install it from GitHub instead:

``` r

install.packages("remotes")
remotes::install_github("sbresnahan/metacleaner", dependencies = TRUE)
```

Installing from GitHub compiles the package on your computer, so it
needs build tools:
[Rtools](https://cran.r-project.org/bin/windows/Rtools/) on Windows, or
the Xcode command line tools on macOS (run `xcode-select --install` in
Terminal). Most Linux systems already have them.

## Your files

You need:

- **query**: the FASTA file you want to clean. Each header must start
  with an NCBI accession, e.g. `>MW545606.1` (anything after the first
  space is ignored). Sequences may be on one line or wrapped over
  several lines.
- **reference**: a FASTA file of sequences from the group you want to
  keep.
- **contaminants** (optional but recommended): a FASTA file of sequences
  from groups you want to remove.
- **taxonomy**: where to find the taxonomy of each accession. This can
  be a taxonomy file (for example the `.tax` file made with MetaCurator
  or taxonomizr), or simply `"ncbi"` to look it up online.

[`build_reference_db()`](https://seantbresnahan.com/metacleaner/reference/build_reference_db.md)
can create the reference and contaminant files for you; see
[`vignette("building-databases")`](https://seantbresnahan.com/metacleaner/articles/building-databases.md).

## A worked example

metacleaner includes a small real example: plant ITS2 sequences from
NCBI (the query), a plant reference set, a fungal contaminant set and a
taxonomy file.
[`system.file()`](https://rdrr.io/r/base/system.file.html) just finds
these example files on your computer; for your own data, put your file
paths in quotes instead, e.g. `query = "C:/data/ITS2_database.fasta"`.

``` r

library(metacleaner)
ex <- function(f) system.file("extdata", f, package = "metacleaner")

result <- clean_reference(
  query        = ex("query.fasta"),
  reference    = ex("reference.fasta"),
  contaminants = ex("contaminants.fasta"),
  taxonomy     = ex("taxonomy.tax"),
  level        = "genus"
)
#> metacleaner: reading sequences
#>   52 query, 168 reference, 33 contaminant sequences
#> Step 1/3: searching the contaminant database
#> Step 2/3: searching the reference database
#> Step 3/3: comparing taxonomy at genus level
#>   3 sequence(s) match the contaminant database
#>   retrieving taxonomy for 171 sequence(s)
#>   0 sequence(s) have no qualifying match in the reference database (0 more rescued by their best hit)
#>   6 sequence(s) match more than one genus; 0 match a different genus
#> Done: kept 40 of 52 sequences, flagged 12 (0.2 s).
result
#> <metacleaner result>
#>   query sequences:   52
#>   kept:              40
#>   flagged:           12
#>     matches multiple genus                   6
#>     matches contaminant database             3
#>     no taxonomy for query at genus           3
#>   settings: identity >= 100%, coverage >= 100%, level = genus, engine = native
#> Use write_results() to save the cleaned database and report.
```

The flagged sequences, with the reason and the hit that caused it:

``` r

head(result$flagged[, c("id", "reason", "top_hit", "hit_taxa")])
#>           id                         reason    top_hit              hit_taxa
#> 1 KF984840.1   matches contaminant database KF984840.1                  <NA>
#> 2 MN228964.1 no taxonomy for query at genus MN228964.1             <unknown>
#> 3 MW545606.1         matches multiple genus MW545600.1 Chrysanthemum; Ajania
#> 4 MW545611.1         matches multiple genus MW545600.1 Chrysanthemum; Ajania
#> 5 MW545649.1         matches multiple genus MW545607.1 Ajania; Chrysanthemum
#> 6 MW545773.1         matches multiple genus MW545787.1 Artemisia; Filifolium
```

Save the cleaned database and the report:

``` r

write_results(result, outdir = "cleaned", prefix = "ITS2")
```

This writes four files into the folder `cleaned`:

| File | Contents |
|----|----|
| `ITS2_clean.fasta` | the sequences that were kept |
| `ITS2_clean.tax` | their taxonomy (same format as the old metacleaner / MetaCurator input) |
| `ITS2_flagged.tsv` | every removed sequence, its taxonomy, the hit that caused the removal, that hit’s taxonomy and the reason; open it in Excel |
| `ITS2_summary.txt` | the settings used and the number removed for each reason |

You can also do it all in one go by giving `outdir` to
[`clean_reference()`](https://seantbresnahan.com/metacleaner/reference/clean_reference.md).

## Choosing settings

- `identity` and `coverage` (default 100 and 100). Lowering them,
  e.g. to `identity = 99, coverage = 95`, also counts near-identical
  matches, so more sequences are checked against more relatives. The
  search is slower and stricter.
- `level` (default `"genus"`). The rank at which the query and its
  matches must agree. Use `"family"` for a more lenient clean or
  `"species"` for a stricter one.
- `require`, e.g. `require = c(class = "Magnoliopsida")`, removes any
  sequence that matches something outside flowering plants.
- `rescue` (default `TRUE`). A sequence without a full-length identical
  match is kept if its closest match is in the same genus. Set
  `rescue = FALSE` to remove all such sequences.
- `threads`. Use more CPU cores for large databases, e.g. `threads = 8`,
  or `threads = "auto"` to use all but one.

## The reasons a sequence is removed

| Reason | Meaning |
|----|----|
| matches contaminant database | identical (at your thresholds) to a sequence in the contaminant database |
| no match in reference database | not found in the reference database |
| hit outside required … | matches something outside the taxon given in `require` |
| no taxonomy for query at … | the sequence’s own taxonomy is unknown at `level` |
| matches multiple … | it matches reference sequences from more than one genus (or other `level`) |
| matches different … | all its matches belong to one other genus |

The checks run in this order and each sequence gets the first reason
that applies; see
[`?classify_queries`](https://seantbresnahan.com/metacleaner/reference/classify_queries.md)
for details.

## Coming from the command-line version

| metacleaner.py option | R argument |
|----|----|
| `-q/--query` | `query` |
| `-x/-f` (badblastdb) | `contaminants` (a FASTA file; no BLAST database needed) |
| `-y/-g` (goodblastdb) | `reference` |
| `-d/--taxadb` | `taxonomy = "path/to/accessionTaxa.sql"` |
| `-e/--pident` | `identity` |
| `-v/--qcovs` | `coverage` |
| `-l/--filterlevel` | `level` |
| `-w/--addfilter class,Magnoliopsida` | `require = c(class = "Magnoliopsida")` |
| `-t/--threads` | `threads` |
| `-o/--outdir` | `outdir` |
| `-c/--chunks`, `-s/--sortonly` | not needed |

To reproduce the BLAST-based search of the old version exactly, install
BLAST+ and use `engine = "blast"`.
