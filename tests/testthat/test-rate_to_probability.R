test_that("rate_to_probability works", {
  expect_equal(round(rate_to_probability(2),2), 0.86)
  expect_equal(round(rate_to_probability(2,2),2), 0.98)
  expect_equal(signif(rate_to_probability(3.5e-5), 2), 3.5e-5)
  expect_error(rate_to_probability('wrong'),
               'rate must be numeric, but is a character')
  expect_error(rate_to_probability(1, -1),
               'time must be a positive number')
})
