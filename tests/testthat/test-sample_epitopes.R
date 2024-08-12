test_that("sample_epitopes works", {
  set.seed(20240812)
  expect_message(epitopes <- sample_epitopes(get_epitope_frequencies(env_features$position)))
  expect_equal(colnames(epitopes),
               c('epi_start_nt', 'epi_end_nt', 'max_fitness_cost'))
  expect_equal(sum(is.na(epitopes$epi_start_nt)), 0)
  expect_equal(sum(is.na(epitopes$epi_end_nt)), 0)
  expect_equal(epitopes$max_fitness_cost,
               c(0.1, 0.133333333333333, 0.166666666666667, 0.2, 0.233333333333333,
                 0.266666666666667, 0.3, 0.333333333333333, 0.366666666666667,
                 0.4))
  expect_message(sample_epitopes(get_epitope_frequencies(env_features$position), cost_type = 'random'))
  expect_message(sample_epitopes(get_epitope_frequencies(env_features$position), end_aa_pos = 1000))
  expect_error(sample_epitopes(get_epitope_frequencies(env_features$position), max_resamples = 1),
               "Too many resamples required.")
  expect_message(sample_epitopes(get_epitope_frequencies(env_features$position),
                  end_aa_pos = (7785-6225)/3,
                  ref_founder_map = map_ref_founder(hxb2_founder, names(hxb2_founder)[1], names(hxb2_founder)[2])))
  expect_message(sample_epitopes(get_epitope_frequencies(env_features$position),
                                 aa_epitope_length = 11))
})

test_that("get_epitope_frequencies works", {
  expect_no_error(epi_freqs <- get_epitope_frequencies(c(1,3,3)))
  expect_equal(epi_freqs$aa_position, 1:3)
  expect_equal(epi_freqs$n_features, c(1,0,2))
  expect_equal(epi_freqs$epitope_probability, c(1/3, 0, 2/3))
})

test_that("convert_ref_to_founder_epitopes works", {
  expect_equal(convert_ref_to_founder_epitopes(tibble::tibble(epi_start_nt = 1,
                                                 epi_end_nt = 10,
                                                 max_fitness_cost = 0.1),
                                  tibble::tibble(ref_pos = 0:10,
                                                 founder_pos = 1:11)),
               tibble::tibble(epi_start_nt=2,
                              epi_end_nt=11,
                              max_fitness_cost=0.1))
})
