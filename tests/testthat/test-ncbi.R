test_that("existing files are not overwritten without resume", {
  skip_if_not_installed("rentrez")
  f <- tempfile(fileext = ".fasta")
  writeLines(c(">a", "ACGT"), f)
  expect_error(download_ncbi_sequences("ITS2", f, verbose = FALSE), "resume = TRUE")
})

test_that("NCBI download and database building work (online)", {
  skip_on_cran()
  skip_if_offline()
  skip_if_not_installed("rentrez")
  term <- 'Ajania[Organism] AND ITS2 AND ("200"[SLEN] : "800"[SLEN])'
  f <- tempfile(fileext = ".fasta")
  download_ncbi_sequences(term, f, batch_size = 2, max_records = 3, verbose = FALSE)
  s <- read_fasta(f)
  expect_length(s, 3)
  expect_true(all(nchar(s) >= 200))
  # resume does not duplicate records that are already present
  download_ncbi_sequences(term, f, batch_size = 2, max_records = 3, resume = TRUE,
                          verbose = FALSE)
  expect_length(read_fasta(f), 3)
  out <- tempfile(fileext = ".fasta")
  build_reference_db(term, out, dereplicate = TRUE, min_length = 200, max_records = 3,
                     verbose = FALSE)
  s2 <- read_fasta(out)
  expect_true(length(s2) >= 1 && length(s2) <= 3 && all(nchar(s2) >= 200))
  expect_false(file.exists(sub("\\.fasta$", ".raw.fasta", out)))
  out2 <- tempfile(fileext = ".fasta")
  expect_error(build_reference_db(term, out2, min_length = 1e6, max_records = 3,
                                  verbose = FALSE), "min_length")
  expect_false(file.exists(out2))
})
