# Changelog

## metacleaner 2.0.0

First release as an R package. This is a complete rewrite of the
Python/R command-line scripts (metacleaner.py + taxafilter.R) released
as v1.0.0 on GitHub.

### New features

- One function,
  [`clean_reference()`](https://sbresnahan.github.io/metacleaner/reference/clean_reference.md),
  runs the whole workflow;
  [`write_results()`](https://sbresnahan.github.io/metacleaner/reference/write_results.md)
  saves the cleaned FASTA, taxonomy, a report of flagged sequences and a
  summary.
- Built-in C++ search engine: BLAST+, Python, pyfasta and shell tools
  are no longer required. Multithreaded on all platforms (`threads`, or
  `threads = "auto"`). BLAST+ can still be used with `engine = "blast"`.
- Taxonomy can come from a taxonomy file (taxonomizr/MetaCurator CSV,
  QIIME-style lineage strings or a table), from NCBI web services
  (`taxonomy = "ncbi"`, only the needed accessions are looked up), from
  a local taxonomizr database, or from several of these in turn.
- [`build_reference_db()`](https://sbresnahan.github.io/metacleaner/reference/build_reference_db.md)
  and
  [`download_ncbi_sequences()`](https://sbresnahan.github.io/metacleaner/reference/download_ncbi_sequences.md)
  replace the shell recipe for building reference and contaminant
  databases;
  [`read_fasta()`](https://sbresnahan.github.io/metacleaner/reference/read_fasta.md),
  [`write_fasta()`](https://sbresnahan.github.io/metacleaner/reference/write_fasta.md),
  [`dereplicate()`](https://sbresnahan.github.io/metacleaner/reference/dereplicate.md)
  and
  [`extract_accession()`](https://sbresnahan.github.io/metacleaner/reference/extract_accession.md)
  replace the awk/perl/sed/seqkit steps.
- FASTA files may be wrapped or gzip-compressed.

### Bug fixes relative to the command-line version

- Queries with no BLAST hit at all in the reference database were never
  flagged; they were written to the clean FASTA but were missing from
  the clean taxonomy file. They are now flagged as “no match in
  reference database”.
- BLAST output was processed in blocks of 10,000 lines, so a query whose
  hits were split across two blocks could be judged on only part of its
  hits.
- The automatic taxonomizr database build could never run.
- The taxonomy step crashed when no sequence was flagged.
- A dot in a directory name broke the search step.
- Flag reasons always said “genera”, whatever `--filterlevel` was used.
- Database names were used as regular expressions to split the flagged
  list, so similar database names (e.g. `ITS` and `ITS_plant`) were
  confused.
- The `superkingdom` rank no longer exists in current taxonomizr
  databases (it is now `domain`); this is handled transparently.
- pyfasta, required for chunking, does not run under Python 3.

### Changes in behavior

- The flagged-sequence report has one row per sequence (previously one
  row per hit) and new columns `database`, `identity`, `coverage`,
  `n_hits` and `hit_taxa`.
- Up to `max_hits = 100` qualifying hits per query are used in the
  taxonomy checks (previously BLAST’s `-max_target_seqs 10`).
