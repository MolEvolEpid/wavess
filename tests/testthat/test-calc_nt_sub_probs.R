test_that("calc_nt_sub_probs works", {
  # this is a hack to get the data in the right format...
  hiv_env_flt_2021 <- ape::as.matrix.DNAbin(hiv_env_flt_2021)
  capture.output(probs <- calc_nt_sub_probs(hiv_env_flt_2021[1:3,]), file = nullfile())
  expect_equal(round(probs[1,2], 2), 0.24)
  expect_equal(round(probs[2,1], 2), 0.43)
  expect_error(calc_nt_sub_probs('not_an_aln'),
               'aln must be of the')
  expect_error(calc_nt_sub_probs(hiv_env_flt_2021[1:3,], 'tr'),
               'tr must be of the')
})
