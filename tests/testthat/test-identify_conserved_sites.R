test_that("get_seq_pos works", {
  expect_equal(
    get_seq_pos(tibble::tibble(seq = c("a", "-", "t")), "seq")$seq_pos,
    c(1, NA, 2)
  )
})

test_that("map_ref_founder works", {
  dnamat <- matrix(c("a", "t", "c", "-", "a", "-", "c", "g"), nrow = 2, byrow = TRUE)
  rownames(dnamat) <- c("seq1", "seq2")
  aln <- ape::as.DNAbin(dnamat)
  map <- map_ref_founder(aln, "seq1", "seq2")
  expect_equal(map$alignment_pos, 1:4)
  expect_equal(map$ref_pos, c(1:3, NA))
  expect_equal(map$founder_pos, c(1, NA, 2:3))
  expect_equal(map$ref_base, c("a", "t", "c", "-"))
  expect_equal(map$founder_base, c("a", "-", "c", "g"))
})

test_that("find_consensus works", {
  dnamat <- matrix(c("-", "a", "t", "c", "-", "a", "a", "-", "c", "g"), nrow = 2, byrow = TRUE)
  rownames(dnamat) <- c("seq1", "seq2")
  aln <- ape::as.DNAbin(dnamat)
  cons <- find_consensus(aln, "seq1")
  expect_equal(cons$founder_pos, 1:3)
  expect_equal(cons$founder_base, c("a", "t", "c"))
  expect_equal(cons$consensus_base, c("a", "-", "c"))
  expect_equal(cons$consensus_prop, c(1, 0.5, 1))
  dnamat2 <- matrix(c("a", "-", "t", "c", "-", "a", "a", "t", "c", "g"), nrow = 2, byrow = TRUE)
  rownames(dnamat2) <- c("seq1", "seq3")
  aln2 <- ape::as.DNAbin(dnamat2)
  cons <- find_consensus(aln, "seq3", "seq1", aln2)
  expect_equal(cons$founder_pos, 1:5)
  expect_equal(cons$founder_base, c("a", "a", "t", "c", "g"))
  expect_equal(cons$consensus_base, c("a", NA, "-", "c", NA))
  expect_equal(cons$consensus_prop, c(1, NA, 0.5, 1, NA))
})

test_that("identify_conserved_sites works", {
  dnamat <- matrix(c("-", "a", "-", "t", "a", "a", "a", "c", "-", "g"), nrow = 2, byrow = TRUE)
  rownames(dnamat) <- c("seq1", "seq2")
  aln <- ape::as.DNAbin(dnamat)
  expect_equal(identify_conserved_sites(aln, "seq1")$conserved, c("Yes", "No", "No"))
  expect_equal(identify_conserved_sites(aln, "seq1", thresh = 0.3)$conserved, c("Yes", NA, "Yes"))
})
