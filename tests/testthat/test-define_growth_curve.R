test_that("define_growth_curve works", {
  gc <- define_growth_curve()
  expect_equal(nrow(gc), 3001)
  expect_equal(colnames(gc), c('generation', 'active_cell_count'))
  expect_equal(gc$generation, 0:3000)
  expect_equal(gc$active_cell_count[1], 1)
  expect_equal(round(gc$active_cell_count[30]/2000, 1), 0.9)
  expect_equal(define_growth_curve(curve_type = 'constant')$active_cell_count,
               rep(2000, 3001))
  expect_equal(define_growth_curve(curve_type = 'linear')$active_cell_count[3001],
               2000)
})
