# Download the real NCBI data used to validate metacleaner 2.0.0 and to derive
# the small example dataset shipped in inst/extdata.
# Usage: Rscript data-raw/fetch_validation_data.R <outdir>
library(metacleaner)
args <- commandArgs(TRUE)
out <- if (length(args)) args[1] else "validation_data"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
fams <- "(Rosaceae[Organism] OR Ericaceae[Organism] OR Asteraceae[Organism] OR Fabaceae[Organism])"
terms <- list(
  query = paste("ITS2[All Fields] AND", fams, 'AND ("150"[SLEN] : "800"[SLEN])'),
  reference = paste("(ITS1[All Fields] OR ITS2[All Fields]) AND", fams,
                    'AND ("150"[SLEN] : "3000"[SLEN])'),
  contaminants = 'ITS2[All Fields] AND Fungi[Organism] AND ("150"[SLEN] : "800"[SLEN])'
)
n <- c(query = 3000, reference = 12000, contaminants = 5000)
for (nm in names(terms)) {
  f <- file.path(out, paste0(nm, ".raw.fasta"))
  download_ncbi_sequences(terms[[nm]], f, max_records = n[[nm]], resume = file.exists(f))
}
ids <- unique(unlist(lapply(names(terms), function(nm)
  names(read_fasta(file.path(out, paste0(nm, ".raw.fasta")))))))
tax <- taxonomy_from_ncbi(ids, cache = file.path(out, "ncbi_taxonomy_cache.rds"))
write.table(tax, file.path(out, "taxonomy_table.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE, na = "NA")
writeLines(c(format(Sys.time()), unlist(terms)), file.path(out, "download_info.txt"))
