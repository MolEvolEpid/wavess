test_that("calc_nt_sub_probs works", {
  # this is a hack to get the data in the right format...
  hiv_env_flt_2022 <- ape::as.matrix.DNAbin(hiv_env_flt_2022)
  capture.output(probs <- calc_nt_sub_probs(hiv_env_flt_2022[1:3, ]),
    file = nullfile()
  )
  expect_equal(unname(unlist(round(probs[1, 3], 2))), 0.23)
  expect_equal(unname(unlist(round(probs[2, 2], 2))), 0.43)
  expect_error(
    calc_nt_sub_probs("not_an_aln"),
    "aln must be of the"
  )
  expect_error(
    calc_nt_sub_probs(hiv_env_flt_2022[1:3, ], "tr"),
    "tr must be of the"
  )
})
