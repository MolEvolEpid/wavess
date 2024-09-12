test_that("check_generate_pop_samp_inputs works", {
  expect_error(generate_pop_samp('not_a_curve'),
               "`curve_type` must be logistic or constant. You provided: not_a_curve")
  expect_error(generate_pop_samp(gN = 'a'), 'gN must be numeric, but is a character')
  expect_error(generate_pop_samp(gN = -1), 'gN must be a positive')
  expect_error(generate_pop_samp(K = 'a'), 'K must be numeric, but is a character')
  expect_error(generate_pop_samp(K = -1), 'K must be a positive')
  expect_error(generate_pop_samp(n0 = 'a'), 'n0 must be numeric, but is a character')
  expect_error(generate_pop_samp(n0 = 10000), 'n0 must be a number \u2264K, but is 10000')
  expect_error(generate_pop_samp(g50 = 'a'), 'g50 must be numeric, but is a character')
  expect_error(generate_pop_samp(g50 = 10000), 'g50 must be a number \u2264gN, but is 10000')
  expect_error(generate_pop_samp(sampling_frequency = 'a'), 'sampling_frequency must be numeric, but is a character')
  expect_error(generate_pop_samp(sampling_frequency = -1), 'sampling_frequency must be a positive number')
  expect_error(generate_pop_samp(sampling_frequency = 'a'), 'sampling_frequency must be numeric, but is a character')
  expect_error(generate_pop_samp(sampling_frequency = 10000), 'sampling_frequency must be ')
  expect_error(generate_pop_samp(max_samp = 'a'), 'max_samp must be numeric, but is a character')
  expect_error(generate_pop_samp(max_samp = -1), 'max_samp must be a positive number')
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

test_that("check_sample_epitopes_inputs works", {
  expect_no_error(check_sample_epitopes_inputs(get_epitope_frequencies(env_features$position),
                                               start_aa_pos = 1,
                                               end_aa_pos = NULL,
                                               num_epitopes = 10,
                                               aa_epitope_length = 10,
                                               max_fit_cost = 0.4,
                                               cost_type = 'linear',
                                               max_resamples = 100,
                                               ref_founder_map = NULL))
  expect_error(check_sample_epitopes_inputs(tibble::tibble(test=0:1),
                                               start_aa_pos = 1,
                                               end_aa_pos = NULL,
                                               num_epitopes = 10,
                                               aa_epitope_length = 10,
                                               max_fit_cost = 0.4,
                                               cost_type = 'linear',
                                               max_resamples = 100,
                                               ref_founder_map = NULL),
                  "The following columns must be included in ")
  expect_error(check_sample_epitopes_inputs(get_epitope_frequencies(env_features$position),
                                               start_aa_pos = 1,
                                               end_aa_pos = NULL,
                                               num_epitopes = 10,
                                               aa_epitope_length = 10,
                                               max_fit_cost = 0.4,
                                               cost_type = 'test',
                                               max_resamples = 100,
                                               ref_founder_map = NULL),
               'Cost type must be either "linear" or "random", but you supplied: test')
  expect_error(check_sample_epitopes_inputs(get_epitope_frequencies(env_features$position),
                                               start_aa_pos = 1,
                                               end_aa_pos = -1,
                                               num_epitopes = 10,
                                               aa_epitope_length = 10,
                                               max_fit_cost = 0.4,
                                               cost_type = 'linear',
                                               max_resamples = 100,
                                               ref_founder_map = NULL),
               "end_aa_pos must be a positive number")
  expect_error(check_sample_epitopes_inputs(get_epitope_frequencies(env_features$position),
                                               start_aa_pos = 1,
                                               end_aa_pos = NULL,
                                               num_epitopes = 10,
                                               aa_epitope_length = 10,
                                               max_fit_cost = 0.4,
                                               cost_type = 'linear',
                                               max_resamples = 100,
                                               ref_founder_map = tibble::tibble(test=0:1)),
               "The following columns must be included in ")
  expect_no_error(check_sample_epitopes_inputs(get_epitope_frequencies(env_features$position),
                                               start_aa_pos = 1,
                                               end_aa_pos = NULL,
                                               num_epitopes = 10,
                                               aa_epitope_length = 10,
                                               max_fit_cost = 0.4,
                                               cost_type = 'linear',
                                               max_resamples = 100,
                                               ref_founder_map = tibble::tibble(ref_pos=0:1,
                                                                                founder_pos=0:1)))
})

test_that('check_calc_nt_subst_probs_inputs works', {
  expect_no_error(check_calc_nt_subst_probs_inputs(hiv_env_flt_2021, ape::rtree(3), 'GTR+I+R(4)', 'none'))
  expect_error(check_calc_nt_subst_probs_inputs('hiv_env_flt_2021', ape::rtree(3), 'GTR+I+R(4)', 'none'),
               'aln must be of the ')
  expect_error(check_calc_nt_subst_probs_inputs(hiv_env_flt_2021, 'ape::rtree(3)', 'GTR+I+R(4)', 'none'),
               'tr must be of the ')
  expect_error(check_calc_nt_subst_probs_inputs(hiv_env_flt_2021, ape::rtree(3), 0, 'none'),
               'model must be a string, but is a numeric')
  expect_error(check_calc_nt_subst_probs_inputs(hiv_env_flt_2021, ape::rtree(3), 'GTR+I+R(4)', 'wrong'),
               'Rearrangement must be one of ')
})

test_that('check_run_wavess_inputs works', {
  hiv_env_flt_2021 <- ape::as.matrix.DNAbin(hiv_env_flt_2021)
  ps <- define_sampling_scheme(define_growth_curve(gN = 100), sampling_frequency = 50)
  fs <- c('ACGT','ATTT')
  suppressMessages(el <- sample_epitopes(get_epitope_frequencies(env_features$position)))
  capture.output(ntsp <- calc_nt_subst_probs(hiv_env_flt_2021[1:3,]), file = nullfile())
  expect_no_error(check_run_wavess_inputs(ps, fs, ntsp, NULL, 0.99, NULL, 1,
                                   NULL, 30, 0.01, 90,
                                   3.5e-5, 1.4e-5,
                                   0.001, 0.01, 0.01, 0.01, NULL))
  expect_error(check_run_wavess_inputs('ps', fs, ntsp, NULL, 0.99, NULL, 1,
                   NULL, 30, 0.01, 90,
                   3.5e-5, 1.4e-5,
                   0.001, 0.01, 0.01, 0.01, NULL),
        'pop_samp must be a data frame or tibble, but is a character')
  expect_error(check_run_wavess_inputs(ps |> dplyr::rename(gen = generation), fs, ntsp, NULL, 0.99, NULL, 1,
                                NULL, 30, 0.01, 90,
                                3.5e-5, 1.4e-5,
                                0.001, 0.01, 0.01, 0.01, NULL),
               'pop_samp must contain the columns generation, active_cell_count, n_sample_active')
  expect_error(check_run_wavess_inputs(ps |> dplyr::mutate(generation = sample(generation)), fs, ntsp, NULL, 0.99, NULL, 1,
                                NULL, 30, 0.01, 90,
                                3.5e-5, 1.4e-5,
                                0.001, 0.01, 0.01, 0.01, NULL),
               'pop_samp')
  expect_error(check_run_wavess_inputs(ps, 'ADAA', ntsp, NULL, 0.99, NULL, 1,
                                NULL, 30, 0.01, 90,
                                3.5e-5, 1.4e-5,
                                0.001, 0.01, 0.01, 0.01, NULL),
               'founder_seqs must only contain the characters ACGT')
  expect_error(check_run_wavess_inputs(ps, c('AG', 'ATAA'), ntsp, NULL, 0.99, NULL, 1,
                         NULL, 30, 0.01, 90,
                         3.5e-5, 1.4e-5,
                         0.001, 0.01, 0.01, 0.01, NULL),
               'All founder sequences must be the same length')
  expect_error(check_run_wavess_inputs(ps, 'ATAA', ntsp |> dplyr::select('A'), NULL, 0.99, NULL, 1,
                                NULL, 30, 0.01, 90,
                                3.5e-5, 1.4e-5,
                                0.001, 0.01, 0.01, 0.01, NULL),
        'nt_sub_probs must have colnames A,C,G,T')
  expect_error(check_run_wavess_inputs(ps, 'ATAA', ntsp[1,], NULL, 0.99, NULL, 1,
                                NULL, 30, 0.01, 90,
                                3.5e-5, 1.4e-5,
                                0.001, 0.01, 0.01, 0.01, NULL),
               'nt_sub_probs must have rownames A,C,G,T')
  expect_no_error(check_run_wavess_inputs(ps, 'ATAA', ntsp, 1, 0.99, NULL, 1,
                                NULL, 30, 0.01, 90,
                                3.5e-5, 1.4e-5,
                                0.001, 0.01, 0.01, 0.01, NULL))
  expect_no_error(check_run_wavess_inputs(ps, 'ATAA', ntsp, c(1,2), 0.99, NULL, 1,
                                   NULL, 30, 0.01, 90,
                                   3.5e-5, 1.4e-5,
                                   0.001, 0.01, 0.01, 0.01, NULL))
  expect_error(check_run_wavess_inputs(ps, 'ATAA', ntsp, 'a', 0.99, NULL, 1,
                                   NULL, 30, 0.01, 90,
                                   3.5e-5, 1.4e-5,
                                   0.001, 0.01, 0.01, 0.01, NULL),
        'conserved_sites must be numeric, but is a character')
  expect_error(check_run_wavess_inputs(ps, 'ATAA', ntsp, 1, 10, NULL, 1,
                                NULL, 30, 0.01, 90,
                                3.5e-5, 1.4e-5,
                                0.001, 0.01, 0.01, 0.01, NULL),
               'conserved_cost must be in the range')
  expect_no_error(check_run_wavess_inputs(ps, 'ATAA', ntsp, NULL, 0.99, 'ATTT', 1,
                                   NULL, 30, 0.01, 90,
                                   3.5e-5, 1.4e-5,
                                   0.001, 0.01, 0.01, 0.01, NULL))
  expect_error(check_run_wavess_inputs(ps, 'ATAA', ntsp, NULL, 0.99, 'ATT', 1,
                                   NULL, 30, 0.01, 90,
                                   3.5e-5, 1.4e-5,
                                   0.001, 0.01, 0.01, 0.01, NULL),
               'ref_seq must be the same length as the founder sequence')
  expect_no_error(check_run_wavess_inputs(ps, 'ATAA', ntsp, NULL, 0.99, NULL, 1,
                                   el, 30, 0.01, 90,
                                   3.5e-5, 1.4e-5,
                                   0.001, 0.01, 0.01, 0.01, NULL))
  expect_error(check_run_wavess_inputs(ps, 'ATAA', ntsp, NULL, 0.99, NULL, 1,
                                   el |> dplyr::select(epi_start_nt), 30, 0.01, 90,
                                   3.5e-5, 1.4e-5,
                                   0.001, 0.01, 0.01, 0.01, NULL),
                  'epitope_locations must contain the columns epi_start_nt, epi_end_nt, max_fitness_cost')
  expect_error(check_run_wavess_inputs(ps, 'ATAA', ntsp, NULL, 0.99, NULL, 1,
                                   NULL, 30, 0.01, 90,
                                   -1, 1.4e-5,
                                   0.001, 0.01, 0.01, 0.01, NULL),
               'prob_mut must be in the range')
  expect_error(check_run_wavess_inputs(ps, 'ATAA', ntsp, NULL, 0.99, NULL, 1,
                                el, 30, 0.01, -1,
                                3.5e-5, 1.4e-5,
                                0.001, 0.01, 0.01, 0.01, NULL),
        'gen_full_potency must be a number')
  expect_no_error(check_run_wavess_inputs(ps, 'ATAA', ntsp, NULL, 0.99, NULL, 1,
                                   NULL, 30, 0.01, 90,
                                   3.5e-5, 1.4e-5,
                                   0.001, 0.01, 0.01, 0.01, 1))
  expect_error(check_run_wavess_inputs(ps, 'ATAA', ntsp, NULL, 0.99, NULL, 1,
                                   NULL, 30, 0.01, 90,
                                   3.5e-5, 1.4e-5,
                                   0.001, 0.01, 0.01, 0.01, 'a'),
               'seed must be numeric, but is a character')
})
