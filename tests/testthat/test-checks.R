test_that("check_define_growth_curve_inputs works", {
  expect_error(define_growth_curve('not_a_curve'),
               "`curve_type` must be logistic, linear, or constant. You provided: not_a_curve")
  expect_error(define_growth_curve(gN = 'a'), 'gN must be numeric, but is a character')
  expect_error(define_growth_curve(gN = -1), 'gN must be a positive number, but is -1')
  expect_error(define_growth_curve(K = 'a'), 'K must be numeric, but is a character')
  expect_error(define_growth_curve(K = -1), 'K must be a positive number, but is -1')
  expect_error(define_growth_curve(n0 = 'a'), 'n0 must be numeric, but is a character')
  expect_error(define_growth_curve(n0 = 10000), 'n0 must be a number \u2264K, but is 10000')
  expect_error(define_growth_curve(gS = 'a'), 'gS must be numeric, but is a character')
  expect_error(define_growth_curve(gS = 10000), 'gS must be a number \u2264gN, but is 10000')
  expect_error(define_growth_curve(pK = 'a'), 'pK must be numeric, but is a character')
  expect_error(define_growth_curve(pK = 2), 'pK must be in the range')
})

test_that("check_define_sampling_scheme_inputs works", {
  gc <- define_growth_curve()
  expect_no_error(define_sampling_scheme(gc))
  expect_error(define_sampling_scheme('char'), 'growth_curve must be a data frame or tibble, but is a character')
  expect_error(define_sampling_scheme(dplyr::rename(gc, g = generation)),
               '`growth_curve` must contain the columns `generation` and `active_cell_count`, but contains instead: g, active_cell_count')
  expect_error(define_sampling_scheme(gc, sampling_frequency = 'a'), 'sampling_frequency must be numeric, but is a character')
  expect_error(define_sampling_scheme(gc, sampling_frequency = -1), 'sampling_frequency must be a positive number, but is -1')
  expect_error(define_sampling_scheme(gc, sampling_frequency = 'a'), 'sampling_frequency must be numeric, but is a character')
  expect_error(define_sampling_scheme(gc, sampling_frequency = 10000), 'sampling_frequency must be ')
  expect_error(define_sampling_scheme(gc, max_samp = 'a'), 'max_samp must be numeric, but is a character')
  expect_error(define_sampling_scheme(gc, max_samp = -1), 'max_samp must be a positive number, but is -1')
})

test_that("check_is_numeric works", {
  expect_no_error(check_is_numeric(0, 'test'))
  expect_error(check_is_numeric('a', 'test'), 'test must be numeric, but is a character')
})

test_that("check_is_pos works", {
  expect_no_error(check_is_pos(1, 'test'))
  expect_error(check_is_pos(0, 'test'), 'test must be a positive number, but is 0')
})

test_that("check_is_df works", {
  expect_no_error(check_is_df(define_growth_curve(), 'test'))
  expect_error(check_is_df(0, 'test'), 'test must be a data frame or tibble, but is a numeric')
})
