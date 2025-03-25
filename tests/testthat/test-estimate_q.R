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

test_that("rate_to_probability works", {
  expect_equal(round(rate_to_probability(2), 2), 0.86)
  expect_equal(round(rate_to_probability(2, 2), 2), 0.98)
  expect_equal(signif(rate_to_probability(3.5e-5), 2), 3.5e-5)
  expect_error(
    rate_to_probability("wrong"),
    "rate must be numeric, but is a character"
  )
  expect_error(
    rate_to_probability(1, -1),
    "time must be a positive number"
  )
})

test_that("calc_q_from_rates works", {
  expect_equal(
    calc_q_from_rates(hiv_mut_rates, 2.4e-5, 1.2),
    structure(c(
      -1.39356019923586, 0.916781309351891, 2.93388807340711,
      0.549959814571701, 0.165021502856394, -3.20893522566672, 0.0182889027942412,
      1.83380054585274, 1.10021446664849, 0.0915971984019694, -3.31892358055296,
      0.550113237779834, 0.128324229611011, 2.20055671790487, 0.366746604293706,
      -2.93387359824134
    ), dim = c(4L, 4L), dimnames = list(c(
      "A", "C",
      "G", "T"
    ), c("A", "C", "G", "T")))
  )
})
