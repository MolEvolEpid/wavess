test_that("extract_seqs works", {
  gp120 <- slice_aln(hxb2_cons_founder, 6225, 7757)
  expect_equal(
    extract_seqs(gp120, "B.US.2011.DEMB11US006.KC473833", start = 1, end = 20),
    list(founder = "ATGAGAGCGATGGGGATCAT", ref = NULL)
  )
  expect_equal(
    extract_seqs(gp120, "B.US.2011.DEMB11US006.KC473833", start = 1500),
    list(founder = "AATTGAACCATTGGGAATAGCACCCACCAGGGCA", ref = NULL)
  )
  expect_equal(
    extract_seqs(gp120, "B.US.2011.DEMB11US006.KC473833",
      "B.FR.83.HXB2_LAI_IIIB_BRU.K03455",
      start = 1, end = 20
    ),
    list(founder = "ATGAGAGCGATGGGGATCAT", ref = "ATGAGAGTGAAGG-------")
  )
  expect_error(
    extract_seqs("not_aln"),
    "aln must be of the "
  )
  expect_error(
    extract_seqs(gp120, "not_founder"),
    "founder_name must be the name of a sequence in aln, but is not_founder"
  )
  expect_error(
    extract_seqs(gp120, "B.US.2011.DEMB11US006.KC473833", "not_name"),
    "ref_name must be the name of a sequence in aln, but is not_name"
  )
  expect_error(
    extract_seqs(gp120, "B.US.2011.DEMB11US006.KC473833", start = "not_number"),
    "start must be numeric, but is a character"
  )
  expect_error(
    extract_seqs(gp120, "B.US.2011.DEMB11US006.KC473833", end = 10000),
    "end must be <= the length of the alignment"
  )
})

test_that("slice_aln works", {
  expect_equal(
    dim(slice_aln(hxb2_cons_founder, start = 1, end = 20)),
    c(3, 20)
  )
  expect_equal(
    dim(slice_aln(hxb2_cons_founder,
      start = 1, end = 20,
      "B.US.2011.DEMB11US006.KC473833"
    )),
    c(1, 20)
  )
})
