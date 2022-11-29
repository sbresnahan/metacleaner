# metacleaner
## Software for cleaning reference sequence databases generated by tools like MetaCurator for metabarcoding and metagenomics  

Reference sequence databases generated by tools like MetaCurator (https://github.com/RTRichar/MetaCurator) sometimes contain falsely labelled references (see: https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.13314). In particular, we found that pollen metabarcoding experiments using plant ITS1 and ITS2 region databases generated by MetaCurator would yield many Sanger reads corresponding to random sequences in fungi, or to ITS1/ITS2 sequences in plants of the wrong genera. 

metacleaner takes as input `.fasta` and `.tax` files generated by tools like MetaCurator, and searches for hits against "good" and "bad" sequence databases to filter undesired accessions before downstream use with tools like MetaBarcoder (https://github.com/RTRichar/MetabarcodeDBsV2).

### Workflow:
1) Query sequences are searched for hits against known undesired sequences using blastn; query sequences with hits above user-defined thresholds for percentage identity (pident) and query coverage (qcovs) are flagged as mislabeled.

2) Query sequences are searched for hits against known desired sequences using blastn. Query sequences with no hits above user-defined thresholds for percentage identity (pident) and query coverage (qcovs) are flagged as mislabeled, while hits above these thresholds are flagged as candidate clean sequences.

3) Taxonomy info for the top hits against desired sequences for the candidate clean sequences is retrieved using taxonomizr (https://github.com/sherrillmix/taxonomizr) and compared against the taxonomy info of the query sequence. If the taxonomy info of the query and subject are not similar at a user-defined level of taxonomy (one of superkingdom, phylum, class, order, family, genus, or species), the query sequence is flagged as mislabeled.

4) Flagged sequences are filtered from the query database.

### Requirements:
python v3.0 or greater, R v.4.0 or greater, pyfasta, and BLAST command line tools
#### Python packages:
pandas, numpy, and biopython
#### R packages:
taxonomizr

### Setup:
metacleaner requires as input:
- A properly formatted `.fasta` file (no line breaks) with NCBI accession IDs as sequence headers
- A `.tax` file containing taxonomy information for the query sequences
- A path to a directory containing either blastdb databases - corresponding to the known "good" and "bad" sequences - generated with `makeblastdb`, or `.fasta` files to create the databases from

Additionally, taxonomizr requires an SQLite database called `accessionTaxa.sql`. The directory containing this file can be specified using `-d` or `--taxadb`, or will otherwise be generated via `taxonomizr::prepareDatabase()` (warning: this requires up to 70Gb and will take several hours to complete). 

### Options:```
-q, --query <(required) path to query fasta file>
-o, --outdir <path to output directory>
-e, --pident <(default=100) percent identity to bad sequence database threshold for filtering>
-v, --qcovs <(default=100) query cover to bad sequence database threshold for filtering>
-r, --runmode <(default=1)>
   (1) = search for query sequences against good and bad blastdbs
   (2) = search for query sequences in only good blastdb
   (3) = search for query sequences in only bad blastdb
-s, --sortonly <(default=F) one of T/F. T = start at sorting step (requires previously generated blastn output files in output directory specified by -o)
-l, --filterlevel <(default=genus) one of superkingdom, phylum, class, order, family, genus, or species>
-b, --blastdbdir <(required) path to badblastdb directory>
-x, --badblastdb <(required) badblastdb>
-f, --badblastdbinput <path to fasta file for badblastdb (required if badblastdb does not already exist in directory specified by -b)>
-y, --goodblastdb <(required) goodblastdb>
-g, --goodblastdbinput <path to fasta file for goodblastdb (required if goodblastdb does not already exist in directory specified by -b)>
-d, --taxadb <(required) path to directory containing accessionTaxa.sql file for taxonomizr>
-c, --chunks <(default=100) number of chunks to split query file into (higher values may increase speed for larger query files)>
-t, --threads <(default=1) number of threads for blastn>
```
