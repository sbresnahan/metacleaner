# Package index

## Clean a reference database

The two functions most users need: run the whole workflow, then save the
cleaned FASTA, taxonomy, list of removed sequences and summary.

- [`clean_reference()`](https://seantbresnahan.com/metacleaner/reference/clean_reference.md)
  : Clean a barcode reference database
- [`write_results()`](https://seantbresnahan.com/metacleaner/reference/write_results.md)
  : Save metacleaner results to files

## Build reference and contaminant databases

Download sequences from NCBI GenBank inside R.

- [`build_reference_db()`](https://seantbresnahan.com/metacleaner/reference/build_reference_db.md)
  : Build a reference (or contaminant) database from NCBI
- [`download_ncbi_sequences()`](https://seantbresnahan.com/metacleaner/reference/download_ncbi_sequences.md)
  : Download sequences from NCBI GenBank

## Taxonomy

Read taxonomy from a file, from NCBI web services, or from a local
taxonomizr database.

- [`read_taxonomy()`](https://seantbresnahan.com/metacleaner/reference/read_taxonomy.md)
  : Read a taxonomy file
- [`taxonomy_from_ncbi()`](https://seantbresnahan.com/metacleaner/reference/taxonomy_from_ncbi.md)
  : Look up taxonomy for accessions from NCBI
- [`taxonomy_from_taxonomizr()`](https://seantbresnahan.com/metacleaner/reference/taxonomy_from_taxonomizr.md)
  : Look up taxonomy for accessions in a local taxonomizr database
- [`taxonomic_ranks()`](https://seantbresnahan.com/metacleaner/reference/taxonomic_ranks.md)
  : Taxonomic ranks used by metacleaner

## Individual steps

The search and decision steps that clean_reference() runs, for advanced
use.

- [`search_sequences()`](https://seantbresnahan.com/metacleaner/reference/search_sequences.md)
  : Search query sequences against a database of sequences
- [`classify_queries()`](https://seantbresnahan.com/metacleaner/reference/classify_queries.md)
  : Decide which query sequences to flag
- [`blast_available()`](https://seantbresnahan.com/metacleaner/reference/blast_available.md)
  : Check whether BLAST+ is installed

## Working with FASTA files

- [`read_fasta()`](https://seantbresnahan.com/metacleaner/reference/read_fasta.md)
  : Read sequences from a FASTA file
- [`write_fasta()`](https://seantbresnahan.com/metacleaner/reference/write_fasta.md)
  : Write sequences to a FASTA file
- [`dereplicate()`](https://seantbresnahan.com/metacleaner/reference/dereplicate.md)
  : Remove duplicate sequences
- [`reverse_complement()`](https://seantbresnahan.com/metacleaner/reference/reverse_complement.md)
  : Reverse complement DNA sequences
- [`extract_accession()`](https://seantbresnahan.com/metacleaner/reference/extract_accession.md)
  : Extract accessions from FASTA header lines

## About the package

- [`metacleaner`](https://seantbresnahan.com/metacleaner/reference/metacleaner-package.md)
  [`metacleaner-package`](https://seantbresnahan.com/metacleaner/reference/metacleaner-package.md)
  : metacleaner: Automated Curation of DNA Barcode Reference Databases
