# Building reference and contaminant databases

metacleaner needs two sets of sequences besides the database you want to
clean:

- a **reference** set of the group you want to keep (e.g. all plant
  ITS1/ITS2 sequences in GenBank), and
- a **contaminant** set of groups that should not be in your database
  (e.g. all non-plant ITS sequences).

[`build_reference_db()`](https://sbresnahan.github.io/metacleaner/reference/build_reference_db.md)
downloads these from NCBI and tidies them in one step: it keeps only the
accession in each header, removes empty sequences and writes each
sequence on a single line. Previously this took a mix of R, `awk`,
`perl`, `sed` and `seqkit` commands.

These downloads need an internet connection and can take hours for large
searches. Get a free NCBI API key (<https://account.ncbi.nlm.nih.gov/>)
to make them about three times faster, and store it once per session:

``` r

Sys.setenv(ENTREZ_KEY = "your_key_here")
```

## Plant ITS1/ITS2 example (Crone et al. 2023)

``` r

library(metacleaner)

# Wanted: plant ITS1/5.8S/28S/ITS2 sequences
build_reference_db(
  term = '((((ITS1) OR 5.8S) OR 28S) OR ITS2) AND Embryophyta[Organism] AND ("0"[SLEN] : "10000"[SLEN])',
  file = "ITS_plants.fasta"
)

# Unwanted: the same markers from everything that is not a land plant
build_reference_db(
  term = '((((ITS1) OR 5.8S) OR 28S) OR ITS2) NOT Embryophyta[Organism] AND ("0"[SLEN] : "10000"[SLEN])',
  file = "ITS_nonplants.fasta",
  dereplicate = TRUE
)
```

The search `term` is written exactly as you would type it into the NCBI
Nucleotide website, so you can test it there first and see how many
records it returns.

If a download stops (for example because the connection dropped), run
the same command again: it continues where it left off.

## Should I remove duplicate sequences?

`dereplicate = TRUE` keeps only one copy of each distinct sequence. This
is fine, and saves time, for the **contaminant** set. Do **not**
dereplicate the **reference** set: if a sequence is shared by *Malus*
and *Pyrus*, keeping only the *Malus* copy hides the ambiguity, and
metacleaner can no longer flag a query as “matches multiple genus”.

## Taxonomy

The easiest option is `taxonomy = "ncbi"` in
[`clean_reference()`](https://sbresnahan.github.io/metacleaner/reference/clean_reference.md).
It looks up only the accessions that are needed. To avoid downloading
them again in later runs, look them up once and save them:

``` r

ids <- c(names(read_fasta("my_database.fasta")),
         names(read_fasta("ITS_plants.fasta")))
tax <- taxonomy_from_ncbi(ids, cache = "taxonomy_cache.rds")
write.table(tax, "taxonomy_table.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
```

You can then use `taxonomy = "taxonomy_table.tsv"`. If you already have
a local taxonomizr database, use
`taxonomy = "path/to/accessionTaxa.sql"`.

## Doing the steps yourself

``` r

download_ncbi_sequences("Malus[Organism] AND ITS2", "malus_raw.fasta")
seqs <- read_fasta("malus_raw.fasta")   # headers reduced to accessions
seqs <- dereplicate(seqs)               # optional
write_fasta(seqs, "malus.fasta")
```
