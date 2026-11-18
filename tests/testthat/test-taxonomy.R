test_that("taxonomizr-style CSV is read", {
  f <- tempfile(fileext = ".tax")
  writeLines(c('" 49992","Z37472.1","49992","Eukaryota","Streptophyta",NA,"Lamiales","Lamiaceae","Thymus","Thymus vulgaris"',
               '"3750","AB1.1","3750","Eukaryota","Streptophyta","Magnoliopsida","Rosales","Rosaceae","Malus","Malus domestica"'), f)
  t <- read_taxonomy(f)
  expect_equal(t$id, c("Z37472.1", "AB1.1"))
  expect_equal(t$taxid, c("49992", "3750"))
  expect_true(is.na(t$class[1]))
  expect_equal(t$species[2], "Malus domestica")
})

test_that("lineage strings with and without prefixes are read", {
  f <- tempfile(fileext = ".txt")
  writeLines(c("A1\tk__Eukaryota;p__Streptophyta;c__Magnoliopsida;o__Rosales;f__Rosaceae;g__Malus;s__Malus domestica",
               "A2\tEukaryota;Streptophyta;Magnoliopsida;Rosales;Rosaceae;Pyrus"), f)
  t <- read_taxonomy(f)
  expect_equal(t$genus, c("Malus", "Pyrus"))
  expect_equal(t$domain, c("Eukaryota", "Eukaryota"))
  expect_true(is.na(t$species[2]))
})

test_that("tables with a header are read, superkingdom becomes domain", {
  f <- tempfile(fileext = ".tsv")
  writeLines(c("accession\tsuperkingdom\tgenus\tspecies", "A1\tEukaryota\tMalus\tMalus x"), f)
  t <- read_taxonomy(f)
  expect_equal(t$domain, "Eukaryota")
  expect_equal(t$genus, "Malus")
  expect_true(is.na(t$family))
})

test_that("lookups match accessions with or without version", {
  tab <- data.frame(id = c("A1", "B2.3"), taxid = NA, domain = "E", phylum = NA,
                    class = NA, order = NA, family = NA, genus = c("Malus", "Pyrus"),
                    species = NA, stringsAsFactors = FALSE)
  l <- metacleaner:::get_lineages(c("A1.1", "B2.3", "C9"), list(tab), verbose = FALSE)
  expect_equal(l$id, c("A1.1", "B2.3"))
  expect_equal(l$genus, c("Malus", "Pyrus"))
})

test_that("rank names are validated", {
  expect_equal(metacleaner:::normalise_rank("Superkingdom"), "domain")
  expect_error(metacleaner:::normalise_rank("clade"), "not a recognized")
})

test_that("NCBI lookup works (online)", {
  skip_on_cran()
  skip_if_offline()
  skip_if_not_installed("rentrez")
  skip_if_not_installed("xml2")
  t <- taxonomy_from_ncbi(c("MW545606.1"), verbose = FALSE)
  expect_equal(nrow(t), 1L)
  expect_false(is.na(t$genus))
})
