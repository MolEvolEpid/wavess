test_that("rate_to_probability works", {
  expect_equal(rate_to_probability(2.1), 0.2571585)
  expect_equal(round(rate_to_probability(2.1, 2), 2), 0.27)
  expect_equal(signif(rate_to_probability(3.5e-5), 2), 3.5e-5)
  expect_error(rate_to_probability('wrong'),
               'rate must be numeric, but is a character')
  expect_error(rate_to_probability(1, -1),
               'k must be a positive number')
})
