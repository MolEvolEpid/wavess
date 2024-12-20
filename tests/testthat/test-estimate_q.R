test_that("estimate_q works", {
  # this is a hack to get the data in the right format...
  hiv_env_flt_2022 <- ape::as.matrix.DNAbin(hiv_env_flt_2022)
  capture.output(q <- estimate_q(hiv_env_flt_2022[1:3, ]),
    file = nullfile()
  )
  expect_equal(round(q[1, 3], 2), 1.11)
  expect_equal(round(q[2, 2], 2), -1.83)
  expect_error(
    estimate_q("not_an_aln"),
    "aln must be of the"
  )
  expect_error(
    estimate_q(hiv_env_flt_2022[1:3, ], "tr"),
    "tr must be of the"
  )
})

test_that("calc_nt_sub_probs_from_q works", {
  hiv_env_flt_2022 <- ape::as.matrix.DNAbin(hiv_env_flt_2022)
  capture.output(q <- estimate_q(hiv_env_flt_2022[1:3, ]),
                 file = nullfile()
  )
  probs <- calc_nt_sub_probs_from_q(q, 3.5e-5)
  expect_equal(round(probs[1, 3], 2), 0.65)
  expect_equal(round(probs[2, 2], 2), 0)
})
