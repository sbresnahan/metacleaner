rand_seq <- function(n, seed) {
  set.seed(seed)
  paste(sample(c("A", "C", "G", "T"), n, replace = TRUE), collapse = "")
}
mutate <- function(s, pos) {
  for (p in pos) {
    b <- substr(s, p, p)
    substr(s, p, p) <- if (b == "A") "C" else "A"
  }
  s
}

base <- rand_seq(300, 1)
db <- c(
  exact_in = paste0(rand_seq(50, 2), base, rand_seq(50, 3)),
  mm3 = mutate(base, c(50, 120, 200)),
  del2 = paste0(substr(base, 1, 150), substr(base, 153, 300)),
  rc = unname(reverse_complement(base)),
  partial = substr(base, 1, 200),
  unrelated = rand_seq(400, 4)
)
q <- c(q1 = base)

test_that("exact search finds containment on both strands", {
  h <- search_sequences(q, db, verbose = FALSE)
  expect_setequal(h$subject, c("exact_in", "rc"))
  expect_true(all(h$identity == 100 & h$coverage == 100 & h$pass))
  expect_equal(h$strand[h$subject == "rc"], "-")
  expect_equal(unique(h$n_pass), 2L)
})

test_that("approximate search reports correct identity and coverage", {
  h <- search_sequences(q, db, identity = 98, coverage = 60, verbose = FALSE)
  expect_setequal(h$subject, c("exact_in", "rc", "mm3", "del2", "partial"))
  expect_equal(h$identity[h$subject == "mm3"], 99)
  expect_equal(h$identity[h$subject == "del2"], 100 * 298 / 300)
  expect_equal(h$coverage[h$subject == "partial"], 100 * 200 / 300)
  h2 <- search_sequences(q, db, identity = 99.5, coverage = 95, verbose = FALSE)
  expect_setequal(h2$subject, c("exact_in", "rc"))
})

test_that("best_hit, max_hits and exclude_self work", {
  q2 <- c(q2 = mutate(base, c(10, 100, 150, 250)))
  h <- search_sequences(q2, db, verbose = FALSE)
  expect_equal(nrow(h), 1L)
  expect_false(h$pass)
  expect_equal(nrow(search_sequences(q2, db, best_hit = FALSE, verbose = FALSE)), 0L)
  h3 <- search_sequences(q, db, identity = 98, coverage = 60, max_hits = 2, verbose = FALSE)
  expect_equal(nrow(h3), 2L)
  expect_equal(unique(h3$n_pass), 5L)
  self <- c(exact_in = db[["exact_in"]])
  expect_equal(nrow(search_sequences(self, db, exclude_self = TRUE, best_hit = FALSE,
                                     verbose = FALSE)), 0L)
})

test_that("threads give identical results", {
  set.seed(5)
  qs <- vapply(1:50, function(i) substr(db[["exact_in"]], i, i + 199), "")
  names(qs) <- paste0("s", 1:50)
  a <- search_sequences(qs, db, identity = 97, coverage = 90, threads = 1, verbose = FALSE)
  b <- search_sequences(qs, db, identity = 97, coverage = 90, threads = 2, verbose = FALSE)
  expect_equal(a, b)
})

test_that("input validation gives helpful errors", {
  expect_error(search_sequences(q, db, identity = 120), "between")
  expect_error(search_sequences("not_a_file.fasta", db), "does not exist")
  expect_error(search_sequences(c("ACGT", "GGCC"), db), "named character vector")
})

test_that("native and BLAST engines agree on simple cases", {
  skip_on_cran()
  skip_if_not(blast_available())
  a <- search_sequences(q, db, identity = 98, coverage = 60, verbose = FALSE)
  b <- search_sequences(q, db, identity = 98, coverage = 60, engine = "blast", verbose = FALSE)
  expect_setequal(a$subject, b$subject)
  m <- match(a$subject, b$subject)
  expect_equal(a$identity, b$identity[m], tolerance = 1e-3)
})
