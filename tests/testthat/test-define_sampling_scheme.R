test_that("define_sampling_scheme() works", {
  ss <- define_sampling_scheme(define_growth_curve())
  expect_equal(ss$n_sample_active[1],1)
  expect_equal(ss$n_sample_active[301], 20)
  expect_error(define_sampling_scheme(tibble::tibble(n=1)),
               "`growth_curve` must contain the columns `generation` and `active_cell_count`, but contains instead: n")
})
