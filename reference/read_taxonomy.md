# Read a taxonomy file

Reads a file that gives the taxonomy (lineage) of each sequence and
returns it as a standard lineage table. Three layouts are recognized
automatically:

## Usage

``` r
read_taxonomy(file, format = c("auto", "taxonomizr", "lineage", "table"))
```

## Arguments

- file:

  Path to the taxonomy file.

- format:

  \`"auto"\` (default), \`"taxonomizr"\`, \`"lineage"\` or \`"table"\`.

## Value

A lineage table: a data frame with columns \`id\`, \`taxid\` and the
seven ranks returned by \[taxonomic_ranks()\]; missing values are
\`NA\`.

## Details

\* \*\*taxonomizr / MetaCurator input format\*\* (also written by
\[write_results()\]): comma-separated, no header, columns \`taxid\`,
\`accession\`, \`taxid\`, then domain, phylum, class, order, family,
genus and species. Example:
\`"49992","Z37472.1","49992","Eukaryota","Streptophyta",NA,"Lamiales",...\`
\* \*\*lineage strings\*\* (QIIME, Metaxa2 and similar): two
tab-separated columns, the ID and a \`;\`-separated lineage, optionally
with rank prefixes such as \`k\_\_\`, \`d\_\_\`, \`p\_\_\`, \`c\_\_\`,
\`o\_\_\`, \`f\_\_\`, \`g\_\_\`, \`s\_\_\`. \* \*\*a table with a header
row\*\* (tab- or comma-separated) with an ID column (\`id\`,
\`accession\`, \`seqid\` or \`qseqid\`) and columns named after ranks
(\`superkingdom\`/\`kingdom\`/\`domain\`, \`phylum\`, ..., \`species\`).

## Examples

``` r
f <- system.file("extdata", "taxonomy.tax", package = "metacleaner")
tax <- read_taxonomy(f)
head(tax)
#>           id   taxid    domain       phylum         class     order     family
#> 1 MZ509216.1  109224 Eukaryota Streptophyta Magnoliopsida   Fabales   Fabaceae
#> 2 MZ509176.1 2909194 Eukaryota Streptophyta Magnoliopsida   Fabales   Fabaceae
#> 3 MW545606.1  463586 Eukaryota Streptophyta Magnoliopsida Asterales Asteraceae
#> 4 MW545649.1  463591 Eukaryota Streptophyta Magnoliopsida Asterales Asteraceae
#> 5 MW545773.1  401928 Eukaryota Streptophyta Magnoliopsida Asterales Asteraceae
#> 6 MW545611.1 2068728 Eukaryota Streptophyta Magnoliopsida Asterales Asteraceae
#>           genus                  species
#> 1   Neustanthus Neustanthus phaseoloides
#> 2      Pueraria        Pueraria sp. KK8s
#> 3        Ajania         Ajania potaninii
#> 4 Chrysanthemum  Chrysanthemum zawadzkii
#> 5     Artemisia      Artemisia pubescens
#> 6        Ajania     Ajania scharnhorstii
```
