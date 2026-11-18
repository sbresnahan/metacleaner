mk_tax <- function(id, genus, class = "Magnoliopsida") {
  data.frame(id = id, taxid = NA_character_, domain = "Eukaryota", phylum = "Streptophyta",
             class = class, order = NA, family = NA, genus = genus, species = NA,
             stringsAsFactors = FALSE)
}
hit <- function(q, s, pass = TRUE, score = 100) {
  data.frame(query = q, subject = s, identity = ifelse(pass, 100, 90),
             coverage = 100, score = score, pass = pass, stringsAsFactors = FALSE)
}

tax <- rbind(
  mk_tax(c("q_ok", "q_multi", "q_wrong", "q_none", "q_rescue", "q_norescue",
           "q_bad", "q_req"), c("Malus", "Malus", "Malus", "Malus", "Pyrus", "Pyrus",
                                "Malus", "Malus")),
  mk_tax(c("r_malus1", "r_malus2", "r_pyrus", "r_prunus"), c("Malus", "Malus", "Pyrus", "Prunus")),
  mk_tax("r_fern", "Malus", class = "Polypodiopsida"),
  mk_tax("q_notax", NA)
)
ref_hits <- rbind(
  hit("q_ok", c("r_malus1", "r_malus2")),
  hit("q_multi", c("r_malus1", "r_pyrus")),
  hit("q_wrong", "r_prunus"),
  hit("q_rescue", "r_pyrus", pass = FALSE),
  hit("q_norescue", "r_prunus", pass = FALSE),
  hit("q_req", c("r_malus1", "r_fern")),
  hit("q_notax", "r_malus1"),
  hit("q_unknownhit", "r_unknown")
)
con_hits <- hit("q_bad", "fungus1")
ids <- c("q_ok", "q_multi", "q_wrong", "q_none", "q_rescue", "q_norescue", "q_bad",
         "q_req", "q_notax")

test_that("each decision rule assigns the expected reason", {
  r <- classify_queries(ids, ref_hits, con_hits, taxonomy = tax, verbose = FALSE)$status
  reason <- stats::setNames(r$reason, r$id)
  expect_true(is.na(reason[["q_ok"]]))
  expect_equal(reason[["q_bad"]], "matches contaminant database")
  expect_equal(reason[["q_none"]], "no match in reference database")
  expect_true(is.na(reason[["q_rescue"]]))
  expect_equal(reason[["q_norescue"]], "no match in reference database")
  expect_equal(reason[["q_multi"]], "matches multiple genus")
  expect_equal(r$top_hit[r$id == "q_multi"], "r_pyrus")
  expect_equal(reason[["q_wrong"]], "matches different genus")
  expect_equal(reason[["q_notax"]], "no taxonomy for query at genus")
  expect_true(is.na(reason[["q_req"]]))
})

test_that("rescue = FALSE, require and level options", {
  r <- classify_queries(ids, ref_hits, con_hits, taxonomy = tax, rescue = FALSE,
                        verbose = FALSE)$status
  expect_equal(r$reason[r$id == "q_rescue"], "no match in reference database")
  r2 <- classify_queries(ids, ref_hits, con_hits, taxonomy = tax,
                         require = c(class = "Magnoliopsida"), verbose = FALSE)$status
  expect_equal(r2$reason[r2$id == "q_req"], "hit outside required class (Magnoliopsida)")
  r3 <- classify_queries(ids, ref_hits, con_hits, taxonomy = tax, require = "class,Magnoliopsida",
                         verbose = FALSE)$status
  expect_equal(r2, r3)
  r4 <- classify_queries(ids, ref_hits, con_hits, taxonomy = tax, level = "class",
                         verbose = FALSE)$status
  expect_true(is.na(r4$reason[r4$id == "q_multi"]))
  expect_equal(r4$reason[r4$id == "q_req"], "matches multiple class")
  expect_error(classify_queries(ids, ref_hits, taxonomy = tax, level = "tribe"), "not a recognized")
  r5 <- classify_queries(ids, ref_hits, taxonomy = tax, level = "superkingdom", verbose = FALSE)
  expect_s3_class(r5$status, "data.frame")
})

test_that("works when nothing is flagged or there are no hits", {
  r <- classify_queries("q_ok", ref_hits, taxonomy = tax, verbose = FALSE)$status
  expect_false(r$flagged)
  r2 <- classify_queries(c("q_ok", "q_none"), ref_hits[0, ], taxonomy = tax, verbose = FALSE)$status
  expect_true(all(r2$reason == "no match in reference database"))
})

test_that("hits with unknown taxonomy count as a separate taxon", {
  h <- rbind(hit("q_ok", "r_malus1"), hit("q_ok", "r_unknown"))
  r <- classify_queries("q_ok", h, taxonomy = tax, verbose = FALSE)$status
  expect_equal(r$reason, "matches multiple genus")
})
