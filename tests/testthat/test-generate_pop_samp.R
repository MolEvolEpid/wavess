test_that("define_growth_curve works", {
  ps <- define_growth_curve()
  expect_equal(dim(ps), c(5001, 2))
  expect_equal(colnames(ps), c("generation", "active_cell_count"))
  expect_equal(ps$generation, 0:5000)
  expect_equal(ps$active_cell_count[1], 10)
  expect_equal(ps$active_cell_count[30], 1915)
  expect_equal(define_growth_curve(n_gen = 300, n0 = 2)$active_cell_count[1], 2)
})


test_that("define_sampling_scheme() works", {
  ss <- define_sampling_scheme()
  expect_equal(ss$day[1], 0)
  expect_equal(ss$day[11], 3650)
  expect_equal(ss$n_sample_active[1], 20)
  expect_equal(ss$n_sample_active[11], 20)
})
