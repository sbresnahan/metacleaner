ex <- function(f) system.file("extdata", f, package = "metacleaner")

test_that("clean_reference runs end to end on the example data", {
  res <- clean_reference(ex("query.fasta"), ex("reference.fasta"), ex("contaminants.fasta"),
                         taxonomy = ex("taxonomy.tax"), verbose = FALSE)
  expect_s3_class(res, "metacleaner_result")
  q <- read_fasta(ex("query.fasta"))
  expect_equal(length(res$clean) + nrow(res$flagged), length(q))
  expect_true(all(names(res$clean) %in% names(q)))
  expect_false(any(names(res$clean) %in% res$flagged$id))
  expect_equal(res$taxonomy$id, names(res$clean))
  expect_true(all(!is.na(res$flagged$reason)))
  expect_output(print(res), "metacleaner result")
  expect_output(print(summary(res)), "kept")
})

test_that("write_results writes consistent files", {
  res <- clean_reference(ex("query.fasta"), ex("reference.fasta"), ex("contaminants.fasta"),
                         taxonomy = ex("taxonomy.tax"), verbose = FALSE)
  d <- tempfile()
  f <- write_results(res, outdir = d, prefix = "t", verbose = FALSE)
  expect_true(all(file.exists(f)))
  fa <- read_fasta(f[["fasta"]])
  tx <- read_taxonomy(f[["taxonomy"]])
  expect_equal(names(fa), names(res$clean))
  expect_equal(tx$id, names(fa))
  fl <- utils::read.delim(f[["flagged"]])
  expect_equal(nrow(fl), nrow(res$flagged))
  f2 <- write_results(res, outdir = d, prefix = "t2", tax_format = "lineage", verbose = FALSE)
  tx2 <- read_taxonomy(f2[["taxonomy"]])
  expect_equal(tx2$genus, tx$genus)
  expect_error(write_results(res, outdir = d, prefix = "t", overwrite = FALSE), "already exist")
})

test_that("queries with no hit anywhere are flagged (bug in metacleaner 1.x)", {
  q <- read_fasta(ex("query.fasta"))
  set.seed(1)
  syn <- c(SYN1 = paste(sample(c("A", "C", "G", "T"), 300, TRUE), collapse = ""))
  res <- clean_reference(c(q[1:5], syn), ex("reference.fasta"), taxonomy = ex("taxonomy.tax"),
                         verbose = FALSE)
  expect_equal(res$flagged$reason[res$flagged$id == "SYN1"], "no match in reference database")
  expect_false("SYN1" %in% names(res$clean))
})

test_that("missing taxonomy argument gives a clear error", {
  expect_error(clean_reference(ex("query.fasta"), ex("reference.fasta")), "taxonomy")
})
