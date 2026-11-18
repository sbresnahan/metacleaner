# Build the small example dataset in inst/extdata from the real NCBI data
# downloaded by data-raw/fetch_validation_data.R.
# Usage: Rscript data-raw/make_example_data.R <validation_data_dir> <extdata_dir>
library(metacleaner)
args <- commandArgs(TRUE)
vd <- if (length(args) >= 1) args[1] else "validation_data"
out <- if (length(args) >= 2) args[2] else "inst/extdata"
set.seed(2023)

ref <- read_fasta(file.path(vd, "reference.raw.fasta"))
con <- read_fasta(file.path(vd, "contaminants.raw.fasta"))
con <- con[!duplicated(names(con))]
tax <- read_taxonomy(file.path(vd, "taxonomy_table.tsv"))

# Clean the whole plant reference set against itself to find real examples
res <- clean_reference(ref, ref, con, taxonomy = tax, threads = 8, verbose = FALSE)
fl <- res$flagged
pick <- unlist(lapply(split(fl$id, fl$reason), function(x) utils::head(sample(x), 6)))
kept <- sample(names(res$clean), 40)
qids <- c(pick, kept)
# Three real fungal ITS2 records are deliberately added to the query set to
# illustrate the contaminant check (a fungal sequence that ended up in a plant
# database). They are also in the contaminant set, with their true taxonomy.
fungal <- sample(names(con), 3)
query <- c(ref[qids], con[fungal])

hits <- res$hits$reference
hits <- hits[hits$query %in% qids, ]
ref_ids <- unique(c(qids, hits$subject))
reference <- ref[ref_ids]
ch <- res$hits$contaminants
con_ids <- unique(c(fungal, ch$subject[ch$query %in% qids], sample(names(con), 30)))
contaminants <- con[con_ids]

dir.create(out, showWarnings = FALSE, recursive = TRUE)
write_fasta(query, file.path(out, "query.fasta"))
write_fasta(reference, file.path(out, "reference.fasta"))
write_fasta(contaminants, file.path(out, "contaminants.fasta"))
tx <- tax[match(unique(c(names(query), names(reference), names(contaminants))), tax$id), ]
tx <- tx[!is.na(tx$id), ]
d <- data.frame(tx$taxid, tx$id, tx$taxid, tx[, taxonomic_ranks()])
utils::write.table(d, file.path(out, "taxonomy.tax"), sep = ",", row.names = FALSE,
                   col.names = FALSE)
print(table(fl$reason[fl$id %in% pick]))
print(file.size(list.files(out, full.names = TRUE)))
