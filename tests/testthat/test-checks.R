test_that("check_define_growth_curve_inputs works", {
  expect_error(define_growth_curve('not_a_curve'),
               "`curve_type` must be logistic or constant. You provided: not_a_curve")
  expect_error(define_growth_curve(gN = 'a'), 'gN must be numeric, but is a character')
  expect_error(define_growth_curve(gN = -1), 'gN must be a positive')
  expect_error(define_growth_curve(K = 'a'), 'K must be numeric, but is a character')
  expect_error(define_growth_curve(K = -1), 'K must be a positive')
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
  expect_error(define_sampling_scheme(gc, sampling_frequency = -1), 'sampling_frequency must be a positive number')
  expect_error(define_sampling_scheme(gc, sampling_frequency = 'a'), 'sampling_frequency must be numeric, but is a character')
  expect_error(define_sampling_scheme(gc, sampling_frequency = 10000), 'sampling_frequency must be ')
  expect_error(define_sampling_scheme(gc, max_samp = 'a'), 'max_samp must be numeric, but is a character')
  expect_error(define_sampling_scheme(gc, max_samp = -1), 'max_samp must be a positive number')
})

test_that("check_is_numeric works", {
  expect_no_error(check_is_numeric(0, 'test'))
  expect_error(check_is_numeric('a', 'test'), 'test must be numeric, but is a character')
})

test_that("check_is_pos works", {
  expect_no_error(check_is_pos(1, 'test'))
  expect_error(check_is_pos(0, 'test'), 'test must be a positive number')
})

test_that("check_is_df works", {
  expect_no_error(check_is_df(define_growth_curve(), 'test'))
  expect_error(check_is_df(0, 'test'), 'test must be a data frame or tibble, but is a numeric')
})

test_that("check_is_string works", {
  expect_no_error(check_is_string('str', 'test'))
  expect_error(check_is_string(0, 'test'), 'test must be a string, but is a numeric')
})

test_that("check_is_0to1 works", {
  expect_no_error(check_is_0to1(0, 'test'))
  expect_no_error(check_is_0to1(1, 'test'))
  expect_no_error(check_is_0to1(0.5, 'test'))
  expect_error(check_is_0to1('wrong', 'test'), 'test must be numeric, but is a character')
  expect_error(check_is_0to1(-1, 'test'), 'test must be in the range')
  expect_error(check_is_0to1(10, 'test'), 'test must be in the range')
})

test_that("check_is_dnabin works", {
  expect_no_error(check_is_dnabin(hxb2_founder, 'hxb2_founder'))
  expect_error(check_is_dnabin(0, 'test'), 'test must be of the ')
})

test_that("check_name_in_alignment works", {
  expect_no_error(check_name_in_alignment(hxb2_founder, 'B.FR.83.HXB2_LAI_IIIB_BRU.K03455'))
  expect_error(check_name_in_alignment(hxb2_founder, 'test', 'hxb2_founder', 'test_name'), 'test_name must be the name of a sequence in hxb2_founder, but is test')
})

test_that("check_identify_conserved_sites_inputs works", {
  expect_no_error(check_identify_conserved_sites_inputs(hiv_env_flt_2021, 'B.US.2011.DEMB11US006.KC473833', thresh = 0.99,
                                                        ref = 'B.FR.83.HXB2_LAI_IIIB_BRU.K03455', founder_aln = hxb2_founder,
                                                        founder_start_pos = 1))
  expect_error(check_identify_conserved_sites_inputs(hiv_env_flt_2021, 'B.US.2011.DEMB11US006.KC473833', thresh = 0.99,
                                                     ref = 'B.FR.83.HXB2_LAI_IIIB_BRU.K03455', founder_aln = NULL,
                                                     founder_start_pos = 1), 'When `ref` is specified, `founder_aln` must also be specified')
  expect_error(check_identify_conserved_sites_inputs(hiv_env_flt_2021, 'B.US.2011.DEMB11US006.KC473833', thresh = 0.99,
                                                     ref = NULL, founder_aln = hxb2_founder,
                                                     founder_start_pos = 1), 'When `founder_aln` is specified, `ref` must also be specified')
  expect_error(check_identify_conserved_sites_inputs(hiv_env_flt_2021, 'B.FR.83.HXB2_LAI_IIIB_BRU.K03455', thresh = 2,
                                                     ref = NULL, founder_aln = NULL,
                                                     founder_start_pos = 1), '`thresh` must be a number between 0 and 1 inclusive')
})

test_that("check_get_seq_pos_inputs works", {
  expect_no_error(check_get_seq_pos_inputs(tibble::tibble(hxb2 = unlist(as.character(hxb2_founder[1]))), 'hxb2'))
})

test_that("check_map_ref_founder_inputs works", {
  expect_no_error(check_map_ref_founder_inputs(hxb2_founder, 'B.FR.83.HXB2_LAI_IIIB_BRU.K03455', 'B.US.2011.DEMB11US006.KC473833'))
})

test_that("check_map_ref_founder_inputs works", {
  expect_no_error(check_sample_epitopes_inputs(get_epitope_frequencies(env_features$position),
                                               start_aa_pos = 1,
                                               end_aa_pos = NULL,
                                               num_epitopes = 10,
                                               aa_epitope_length = 10,
                                               min_max_fit_cost = 0.1,
                                               max_max_fit_cost = 0.4,
                                               cost_type = 'linear',
                                               max_resamples = 100,
                                               ref_founder_map = NULL))
  expect_error(check_sample_epitopes_inputs(tibble::tibble(test=0:1),
                                               start_aa_pos = 1,
                                               end_aa_pos = NULL,
                                               num_epitopes = 10,
                                               aa_epitope_length = 10,
                                               min_max_fit_cost = 0.1,
                                               max_max_fit_cost = 0.4,
                                               cost_type = 'linear',
                                               max_resamples = 100,
                                               ref_founder_map = NULL),
                  "The following columns must be included in ")
  expect_error(check_sample_epitopes_inputs(get_epitope_frequencies(env_features$position),
                                               start_aa_pos = 1,
                                               end_aa_pos = NULL,
                                               num_epitopes = 10,
                                               aa_epitope_length = 10,
                                               min_max_fit_cost = 0.1,
                                               max_max_fit_cost = 0.4,
                                               cost_type = 'test',
                                               max_resamples = 100,
                                               ref_founder_map = NULL),
               'Cost type must be either "linear" or "random", but you supplied: test')
  expect_error(check_sample_epitopes_inputs(get_epitope_frequencies(env_features$position),
                                               start_aa_pos = 1,
                                               end_aa_pos = -1,
                                               num_epitopes = 10,
                                               aa_epitope_length = 10,
                                               min_max_fit_cost = 0.1,
                                               max_max_fit_cost = 0.4,
                                               cost_type = 'linear',
                                               max_resamples = 100,
                                               ref_founder_map = NULL),
               "end_aa_pos must be a positive number")
  expect_error(check_sample_epitopes_inputs(get_epitope_frequencies(env_features$position),
                                               start_aa_pos = 1,
                                               end_aa_pos = NULL,
                                               num_epitopes = 10,
                                               aa_epitope_length = 10,
                                               min_max_fit_cost = 0.1,
                                               max_max_fit_cost = 0.4,
                                               cost_type = 'linear',
                                               max_resamples = 100,
                                               ref_founder_map = tibble::tibble(test=0:1)),
               "The following columns must be included in ")
  expect_no_error(check_sample_epitopes_inputs(get_epitope_frequencies(env_features$position),
                                               start_aa_pos = 1,
                                               end_aa_pos = NULL,
                                               num_epitopes = 10,
                                               aa_epitope_length = 10,
                                               min_max_fit_cost = 0.1,
                                               max_max_fit_cost = 0.4,
                                               cost_type = 'linear',
                                               max_resamples = 100,
                                               ref_founder_map = tibble::tibble(ref_pos=0:1,
                                                                                founder_pos=0:1)))
})
