test_that("generate_pop_samp works", {
  ps <- generate_pop_samp()
  expect_equal(dim(ps), c(3001, 3))
  expect_equal(colnames(ps), c(
    "generation", "active_cell_count",
    "n_sample_active"
  ))
  expect_equal(ps$generation, 0:3000)
  expect_equal(ps$active_cell_count[1], 10)
  expect_equal(ps$active_cell_count[30], 1915)
  expect_equal(generate_pop_samp(n_gen = 300, n0 = 2)$active_cell_count[1], 2)
})

test_that("define_growth_curve works", {
  gc <- define_growth_curve()
  expect_equal(nrow(gc), 3001)
  expect_equal(colnames(gc), c("generation", "active_cell_count"))
  expect_equal(gc$generation, 0:3000)
  expect_equal(gc$active_cell_count[1], 1)
  expect_equal(gc$active_cell_count[10], 452)
})

test_that("define_sampling_scheme() works", {
  ss <- define_sampling_scheme(define_growth_curve())
  expect_equal(ss$n_sample_active[1], 1)
  expect_equal(ss$n_sample_active[301], 20)
})
