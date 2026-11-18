test_that("read_fasta handles wrapped lines, CRLF, lowercase and header words", {
  f <- tempfile(fileext = ".fasta")
  writeLines(c(">AB123456.1 Malus domestica ITS2", "acgt", "ACGT", "",
               ">gi|12345|gb|XY999999.2| something", "TTTT\r",
               ";comment", ">empty", ">Z1", "NNAC GT"), f)
  expect_warning(s <- read_fasta(f), "empty")
  expect_equal(unname(s), c("ACGTACGT", "TTTT", "NNACGT"))
  expect_equal(names(s), c("AB123456.1", "XY999999.2", "Z1"))
  expect_equal(names(suppressWarnings(read_fasta(f, drop_version = TRUE)))[1:2],
               c("AB123456", "XY999999"))
  expect_equal(names(suppressWarnings(read_fasta(f, ids = "header")))[1],
               "AB123456.1 Malus domestica ITS2")
})

test_that("read_fasta reads gzip files and rejects non-FASTA input", {
  s <- c(a = "ACGT", b = "GGCC")
  f <- tempfile(fileext = ".fasta.gz")
  write_fasta(s, f)
  expect_equal(read_fasta(f), s)
  bad <- tempfile()
  writeLines(c("ACGT", ">x", "ACGT"), bad)
  expect_error(read_fasta(bad), "does not look like a FASTA")
  expect_error(read_fasta(tempfile()), "not found")
})

test_that("write_fasta round-trips and wraps", {
  s <- c(a = strrep("ACGT", 30), b = "GG")
  f <- tempfile(fileext = ".fa")
  write_fasta(s, f, width = 50)
  expect_equal(length(readLines(f)), 2 + 3 + 1)
  expect_equal(read_fasta(f), s)
  write_fasta(s, f)
  expect_equal(length(readLines(f)), 4)
})

test_that("reverse complement and dereplication", {
  expect_equal(reverse_complement(c(x = "AACGTN")), c(x = "NACGTT"))
  expect_equal(unname(reverse_complement("RYKMBVDH")), "DHBVKMRY")
  s <- c(a = "AACC", b = "AACC", c = "GGTT", d = "ACGA")
  expect_equal(names(dereplicate(s)), c("a", "d"))
  expect_equal(names(dereplicate(s, both_strands = FALSE)), c("a", "c", "d"))
})

test_that("extract_accession handles NCBI header styles", {
  expect_equal(extract_accession(c(">A1.1 desc", "gi|1|gb|B2.1|", "gb|C3.1|", "D4")),
               c("A1.1", "B2.1", "C3.1", "D4"))
})
