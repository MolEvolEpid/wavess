test_that("calc_div_metrics works", {
  aln1 <- calc_div_metrics(ape::as.DNAbin(list(
    seq1 = c("A", "C", "G", "T"),
    seq2 = c("A", "C", "G", "T"),
    seq3 = c("A", "C", "G", "T"),
    seq4 = c("A", "C", "G", "T")
  )), "seq1", c(seq1 = 0, seq2 = 1, seq3 = 1, seq4 = 1))
  expect_equal(aln1$diversity, c(NaN, 0))
  expect_equal(aln1$divergence, c(0, 0))
  aln2 <- calc_div_metrics(ape::as.DNAbin(list(
    seq1 = c("A", "C", "G", "T"),
    seq2 = c("C", "C", "G", "T"),
    seq3 = c("C", "C", "G", "T"),
    seq4 = c("C", "C", "G", "T")
  )), "seq1", c(seq1 = 0, seq2 = 1, seq3 = 1, seq4 = 1))
  expect_equal(aln2$diversity, c(NaN, 0))
  expect_equal(aln2$divergence, c(0, 0.25))
  aln3 <- calc_div_metrics(ape::as.DNAbin(list(
    seq1 = c("A", "C", "G", "T"),
    seq2 = c("A", "C", "G", "T"),
    seq3 = c("C", "C", "G", "T"),
    seq4 = c("C", "C", "G", "T")
  )), "seq1", c(seq1 = 0, seq2 = 1, seq3 = 1, seq4 = 1))
  expect_equal(aln3$diversity, c(NaN, 1 / 4 * 2 / 3))
  expect_equal(aln3$divergence, c(0, 0.25 * 2 / 3))
  aln4 <- calc_div_metrics(ape::as.DNAbin(list(
    seq1 = c("A", "C", "G", "T"),
    seq2 = c("A", "C", "G", "T"),
    seq3 = c("C", "T", "G", "T"),
    seq4 = c("C", "C", "G", "T")
  )), "seq1", c(seq1 = 0, seq2 = 1, seq3 = 1, seq4 = 1))
  expect_equal(aln4$diversity, c(NaN, 1 / 3))
  expect_equal(aln4$divergence, c(0, 0.25))
  expect_error(
    calc_div_metrics(ape::as.DNAbin(list(
      seq1 = c("A", "C", "G", "T"),
      seq2 = c("A", "C", "G", "T"),
      seq3 = c("A", "C", "G", "T"),
      seq4 = c("A", "C", "G", "T")
    )), "seq1", c(1, 2, 3)),
    "The length of gen must be the same as the number of sequences in the alignment"
  )
})
