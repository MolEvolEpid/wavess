# test_that("check_generate_pop_samp_inputs works", {
#   expect_error(
#     generate_pop_samp(n_gen = "a"),
#     "n_gen must be numeric, but is a character"
#   )
#   expect_error(generate_pop_samp(n_gen = -1), "n_gen must be a positive")
#   expect_error(
#     generate_pop_samp(carry_cap = "a"),
#     "carry_cap must be numeric, but is a character"
#   )
#   expect_error(
#     generate_pop_samp(carry_cap = -1),
#     "carry_cap must be a positive"
#   )
#   expect_error(
#     generate_pop_samp(n0 = "a"),
#     "n0 must be numeric, but is a character"
#   )
#   expect_error(
#     generate_pop_samp(n0 = 10000),
#     "n0 must be a number "
#   )
#   expect_error(
#     generate_pop_samp(max_growth_rate = "a"),
#     "max_growth_rate must be numeric, but is a character"
#   )
#   expect_error(
#     generate_pop_samp(sampling_frequency = "a"),
#     "sampling_frequency must be numeric, but is a character"
#   )
#   expect_error(
#     generate_pop_samp(sampling_frequency = -1),
#     "sampling_frequency must be a positive number"
#   )
#   expect_error(
#     generate_pop_samp(sampling_frequency = "a"),
#     "sampling_frequency must be numeric, but is a character"
#   )
#   expect_error(
#     generate_pop_samp(sampling_frequency = 10000),
#     "sampling_frequency must be "
#   )
#   expect_error(
#     generate_pop_samp(max_samp = "a"),
#     "max_samp must be numeric, but is a character"
#   )
#   expect_error(
#     generate_pop_samp(max_samp = -1),
#     "max_samp must be a positive number"
#   )
# })

test_that("check_is_numeric works", {
  expect_no_error(check_is_numeric(0, "test"))
  expect_error(
    check_is_numeric("a", "test"),
    "test must be numeric, but is a character"
  )
})

test_that("check_is_pos works", {
  expect_no_error(check_is_pos(1, "test"))
  expect_error(check_is_pos(0, "test"), "test must be a positive number")
})

test_that("check_is_df works", {
  expect_no_error(check_is_df(define_growth_curve(), "test"))
  expect_error(
    check_is_df(0, "test"),
    "test must be a data frame or tibble, but is a numeric"
  )
})

test_that("check_is_string works", {
  expect_no_error(check_is_string("str", "test"))
  expect_error(
    check_is_string(0, "test"),
    "test must be a string, but is a numeric"
  )
})

test_that("check_is_0to1 works", {
  expect_no_error(check_is_0to1(0, "test"))
  expect_no_error(check_is_0to1(1, "test"))
  expect_no_error(check_is_0to1(0.5, "test"))
  expect_error(
    check_is_0to1("wrong", "test"),
    "test must be numeric, but is a character"
  )
  expect_error(check_is_0to1(-1, "test"), "test must be in the range")
  expect_error(check_is_0to1(10, "test"), "test must be in the range")
})

test_that("check_is_dnabin works", {
  expect_no_error(check_is_dnabin(hxb2_cons_founder, "hxb2_cons_founder"))
  expect_error(check_is_dnabin(0, "test"), "test must be of the ")
})

test_that("check_name_in_alignment works", {
  expect_no_error(check_name_in_alignment(
    hxb2_cons_founder,
    "B.FR.83.HXB2_LAI_IIIB_BRU.K03455"
  ))
  expect_error(
    check_name_in_alignment(
      hxb2_cons_founder, "test",
      "hxb2_cons_founder", "test_name"
    ),
    "test_name must be the name of a sequence in hxb2_cons_founder, but is test"
  )
})

test_that("check_conserved_sites_inputs works", {
  expect_no_error(check_conserved_sites_inputs(hiv_env_flt_2022,
    "B.US.2011.DEMB11US006.KC473833",
    thresh = 0.99,
    ref = "B.FR.83.HXB2_LAI_IIIB_BRU.K03455", founder_aln = hxb2_cons_founder
  ))
  expect_error(
    check_conserved_sites_inputs(hiv_env_flt_2022,
      "B.US.2011.DEMB11US006.KC473833",
      thresh = 0.99,
      ref = "B.FR.83.HXB2_LAI_IIIB_BRU.K03455", founder_aln = NULL
    ),
    "When `ref` is specified, `founder_aln` must also be specified"
  )
  expect_error(
    check_conserved_sites_inputs(hiv_env_flt_2022,
      "B.US.2011.DEMB11US006.KC473833",
      thresh = 0.99,
      ref = NULL, founder_aln = hxb2_cons_founder
    ),
    "When `founder_aln` is specified, `ref` must also be specified"
  )
  expect_error(
    check_conserved_sites_inputs(hiv_env_flt_2022,
      "B.FR.83.HXB2_LAI_IIIB_BRU.K03455",
      thresh = 2,
      ref = NULL, founder_aln = NULL
    ),
    "`thresh` must be a number between 0 and 1 inclusive"
  )
})

test_that("check_get_seq_pos_inputs works", {
  expect_no_error(check_get_seq_pos_inputs(
    tibble::tibble(hxb2 = unlist(as.character(hxb2_cons_founder[1]))), "hxb2"
  ))
})

test_that("check_map_ref_founder_inputs works", {
  expect_no_error(check_map_ref_founder_inputs(
    hxb2_cons_founder,
    "B.FR.83.HXB2_LAI_IIIB_BRU.K03455", "B.US.2011.DEMB11US006.KC473833"
  ))
})

test_that("check_sample_epitopes_inputs works", {
  expect_no_error(check_sample_epitopes_inputs(
    get_epitope_frequencies(env_features$Position - 1),
    start_aa_pos = 1,
    end_aa_pos = NULL,
    num_epitopes = 10,
    aa_epitope_length = 10,
    max_fit_cost = 0.4,
    max_resamples = 100,
    ref_founder_map = NULL
  ))
  expect_error(
    check_sample_epitopes_inputs(tibble::tibble(test = 0:1),
      start_aa_pos = 1,
      end_aa_pos = NULL,
      num_epitopes = 10,
      aa_epitope_length = 10,
      max_fit_cost = 0.4,
      max_resamples = 100,
      ref_founder_map = NULL
    ),
    "The following columns must be included in "
  )
  expect_error(
    check_sample_epitopes_inputs(get_epitope_frequencies(env_features$Position - 1),
      start_aa_pos = 1,
      end_aa_pos = -1,
      num_epitopes = 10,
      aa_epitope_length = 10,
      max_fit_cost = 0.4,
      max_resamples = 100,
      ref_founder_map = NULL
    ),
    "end_aa_pos must be a positive number"
  )
  expect_error(
    check_sample_epitopes_inputs(get_epitope_frequencies(env_features$Position - 1),
      start_aa_pos = 1,
      end_aa_pos = NULL,
      num_epitopes = 10,
      aa_epitope_length = 10,
      max_fit_cost = 0.4,
      max_resamples = 100,
      ref_founder_map = tibble::tibble(test = 0:1)
    ),
    "The following columns must be included in "
  )
  expect_no_error(check_sample_epitopes_inputs(
    get_epitope_frequencies(env_features$Position - 1),
    start_aa_pos = 1,
    end_aa_pos = NULL,
    num_epitopes = 10,
    aa_epitope_length = 10,
    max_fit_cost = 0.4,
    max_resamples = 100,
    ref_founder_map = tibble::tibble(
      ref_pos = 0:1,
      founder_pos = 0:1
    )
  ))
})

test_that("check_estimate_q_inputs works", {
  expect_no_error(check_estimate_q_inputs(
    hiv_env_flt_2022, ape::rtree(3), "GTR+I+R(4)", "none"
  ))
  expect_error(
    check_estimate_q_inputs(
      "hiv_env_flt_2022", ape::rtree(3), "GTR+I+R(4)", "none"
    ),
    "aln must be of the "
  )
  expect_error(
    check_estimate_q_inputs(
      hiv_env_flt_2022, "ape::rtree(3)", "GTR+I+R(4)", "none"
    ),
    "tr must be of the "
  )
  expect_error(
    check_estimate_q_inputs(hiv_env_flt_2022, ape::rtree(3), 0, "none"),
    "model must be a string, but is a numeric"
  )
  expect_error(
    check_estimate_q_inputs(
      hiv_env_flt_2022, ape::rtree(3), "GTR+I+R(4)", "wrong"
    ),
    "Rearrangement must be one of "
  )
})

test_that("check_run_wavess_inputs works", {
  # hiv_env_flt_2022 <- ape::as.matrix.DNAbin(hiv_env_flt_2022)
  inf_pop_size <- define_growth_curve(n_gen = 100)
  samp_scheme <- define_sampling_scheme(sampling_frequency = 50, n_days = 100)
  fs <- rep("ACGT", 10)
  suppressMessages(el <- sample_epitopes(
    get_epitope_frequencies(env_features$Position - 1)
  ))
  hiv_q_mat <- calc_q_from_rates(hiv_mut_rates, 2.4e-5, 1.2)
  expect_no_error(check_run_wavess_inputs(
    inf_pop_size, samp_scheme, fs, 1.2, hiv_q_mat,
    3.5e-5, 1.4e-5,
    NULL, 0.99, NULL, 1,
    NULL, 30, 0.01, 90,
    0.001, 0.01, 0.01, 0.01, NULL
  ))
  expect_error(
    check_run_wavess_inputs(
      inf_pop_size, samp_scheme |> dplyr::mutate(n_sample_active = 0), fs, 1.2, hiv_q_mat,
      3.5e-5, 1.4e-5,
      NULL, 0.99, NULL, 1,
      NULL, 30, 0.01, 90,
      0.001, 0.01, 0.01, 0.01, NULL
    ),
    "you must sample at least one day"
  )
  expect_error(check_run_wavess_inputs(
    inf_pop_size, samp_scheme, fs, 1.2, hiv_q_mat,
    3.5e-5, 1.4e-5,
    10000, 0.99, NULL, 1,
    NULL, 30, 0.01, 90,
    0.001, 0.01, 0.01, 0.01, NULL
  ))
  expect_error(
    check_run_wavess_inputs(
      "ps", fs, 1.2, hiv_q_mat,
      3.5e-5, 1.4e-5,
      NULL, 0.99, NULL, 1,
      NULL, 30, 0.01, 90,
      0.001, 0.01, 0.01, 0.01, NULL
    ),
    "inf_pop_size must be a data frame or tibble, but is a character"
  )
  expect_error(
    check_run_wavess_inputs(
      inf_pop_size |> dplyr::rename(gen = generation), fs, 1.2, hiv_q_mat,
      3.5e-5, 1.4e-5,
      NULL, 0.99, NULL, 1,
      NULL, 30, 0.01, 90,
      0.001, 0.01, 0.01, 0.01, NULL
    ),
    "inf_pop_size must contain the columns generation, active_cell_count"
  )
  expect_error(
    check_run_wavess_inputs(
      inf_pop_samp |> dplyr::mutate(generation = sample(generation)), fs, 1.2, hiv_q_mat,
      3.5e-5, 1.4e-5,
      NULL, 0.99, NULL, 1,
      NULL, 30, 0.01, 90,
      0.001, 0.01, 0.01, 0.01, NULL
    ),
    "pop_samp"
  )
  expect_error(
    check_run_wavess_inputs(
      inf_pop_size, samp_scheme, rep("ADAA", 10), 1.2, hiv_q_mat,
      3.5e-5, 1.4e-5,
      NULL, 0.99, NULL, 1,
      NULL, 30, 0.01, 90,
      0.001, 0.01, 0.01, 0.01, NULL
    ),
    "founder_seqs must only contain the characters ACGT"
  )
  expect_error(check_run_wavess_inputs(
    inf_pop_size, samp_scheme, fs, 'notanum', hiv_q_mat,
    3.5e-5, 1.4e-5,
    NULL, 0.99, NULL, 1,
    NULL, 30, 0.01, 90,
    0.001, 0.01, 0.01, 0.01, NULL
  ),
  "generation_time must be numeric, but is")
  expect_error(
    check_run_wavess_inputs(
      inf_pop_size, samp_scheme, rep("ATAA", 10), 1.2, hiv_q_mat[1:2, ],
      3.5e-5, 1.4e-5,
      NULL, 0.99, NULL, 1,
      NULL, 30, 0.01, 90,
      0.001, 0.01, 0.01, 0.01, NULL
    ),
    "q must have rownames A,C,G,T"
  )
  hiv_q_mat_tmp <- hiv_q_mat
  colnames(hiv_q_mat_tmp) <- rev(colnames(hiv_q_mat_tmp))
  expect_error(
    check_run_wavess_inputs(
      inf_pop_size, samp_scheme, rep("ATAA", 10), 1.2, hiv_q_mat_tmp,
      3.5e-5, 1.4e-5,
      NULL, 0.99, NULL, 1,
      NULL, 30, 0.01, 90,
      0.001, 0.01, 0.01, 0.01, NULL
    ),
    "the row names and column names of q must be in the same order"
  )
  expect_error(
    check_run_wavess_inputs(
      inf_pop_size, samp_scheme, rep("ATAA", 10), 1.2, hiv_q_mat[1, ],
      3.5e-5, 1.4e-5,
      NULL, 0.99, NULL, 1,
      NULL, 30, 0.01, 90,
      0.001, 0.01, 0.01, 0.01, NULL
    ),
    "q must be a matrix, but is a numeric"
  )
  expect_no_error(check_run_wavess_inputs(
    inf_pop_size, samp_scheme, rep("ATAA", 10), 1.2, hiv_q_mat,
    3.5e-5, 1.4e-5,
    c("1" = "a"), 0.99, NULL, 1,
    NULL, 30, 0.01, 90,
    0.001, 0.01, 0.01, 0.01, NULL
  ))
  expect_no_error(check_run_wavess_inputs(
    inf_pop_size, samp_scheme, rep("ATAA", 10), 1.2, hiv_q_mat,
    3.5e-5, 1.4e-5,
    c("1" = "A", "2" = "C"), 0.99, NULL, 1,
    NULL, 30, 0.01, 90,
    0.001, 0.01, 0.01, 0.01, NULL
  ))
  expect_error(
    check_run_wavess_inputs(
      inf_pop_size, samp_scheme, rep("ATAA", 10), 1.2, hiv_q_mat,
      3.5e-5, 1.4e-5,
      "a", 0.99, NULL, 1,
      NULL, 30, 0.01, 90,
      0.001, 0.01, 0.01, 0.01, NULL
    ),
    "conserved_sites must be a named vector"
  )
  expect_error(
    check_run_wavess_inputs(
      inf_pop_size, samp_scheme, rep("ATAA", 10), 1.2, hiv_q_mat,
      3.5e-5, 1.4e-5,
      c("1" = "A", "2" = "C"), 10, NULL, 1,
      NULL, 30, 0.01, 90,
      0.001, 0.01, 0.01, 0.01, NULL
    ),
    "conserved_cost must be in the range"
  )
  expect_no_error(check_run_wavess_inputs(
    inf_pop_size, samp_scheme, rep("ATAA", 10), 1.2, hiv_q_mat,
    3.5e-5, 1.4e-5,
    NULL, 0.99, "ATTT", 0.99,
    NULL, 30, 0.01, 90,
    0.001, 0.01, 0.01, 0.01, NULL
  ))

  expect_error(check_run_wavess_inputs(
    inf_pop_size, samp_scheme, rep("ATAA", 10), 1.2, hiv_q_mat,
    3.5e-5, 1.4e-5,
    NULL, 0.99, "ATTT", 1,
    NULL, 30, 0.01, 90,
    0.001, 0.01, 0.01, 0.01, NULL
  ), "replicative_cost must be in the range")

  expect_error(
    check_run_wavess_inputs(
      inf_pop_size, samp_scheme, rep("ATAA", 10), 1.2, hiv_q_mat,
      3.5e-5, 1.4e-5,
      NULL, 0.99, "ATT", 0.99,
      NULL, 30, 0.01, 90,
      0.001, 0.01, 0.01, 0.01, NULL
    ),
    "ref_seq must be the same length as the founder sequence"
  )
  expect_no_error(check_run_wavess_inputs(
    inf_pop_size, samp_scheme, rep("ATAA", 10), 1.2, hiv_q_mat,
    3.5e-5, 1.4e-5,
    NULL, 0.99, NULL, 1,
    el, 30, 0.01, 90,
    0.001, 0.01, 0.01, 0.01, NULL
  ))
  expect_error(
    check_run_wavess_inputs(
      inf_pop_size, samp_scheme, rep("ATAA", 10), 1.2, hiv_q_mat,
      3.5e-5, 1.4e-5,
      NULL, 0.99, NULL, 1,
      el |> dplyr::select(epi_start_nt), 30, 0.01, 90,
      0.001, 0.01, 0.01, 0.01, NULL
    ),
    "epitope_locations must contain the columns epi_start_nt, epi_end_nt, "
  )
  expect_error(
    check_run_wavess_inputs(
      inf_pop_size, samp_scheme, rep("ATAA", 10), 1.2, hiv_q_mat[, 1:2],
      3.5e-5, 1.4e-5,
      NULL, 0.99, NULL, 1,
      NULL, 30, 0.01, 90,
      0.001, 0.01, 0.01, 0.01, NULL
    ),
    "q must have colnames A,C,G,T"
  )
  expect_error(
    check_run_wavess_inputs(
      inf_pop_size, samp_scheme, rep("ATAA", 10), 1.2, hiv_q_mat,
      3.5e-5, 1.4e-5,
      NULL, 0.99, NULL, 1,
      el, 30, 0.01, -1,
      0.001, 0.01, 0.01, 0.01, NULL
    ),
    "gen_full_potency must be a positive number"
  )
  expect_no_error(check_run_wavess_inputs(
    inf_pop_size, samp_scheme, rep("ATAA", 10), 1.2, hiv_q_mat,
    3.5e-5, 1.4e-5,
    NULL, 0.99, NULL, 1,
    NULL, 30, 0.01, 90,
    0.001, 0.01, 0.01, 0.01, 1
  ))
  expect_error(
    check_run_wavess_inputs(
      inf_pop_size, samp_scheme, rep("ATAA", 10), 1.2, hiv_q_mat,
      3.5e-5, 1.4e-5,
      NULL, 0.99, NULL, 1,
      NULL, 30, 0.01, 90,
      0.001, 0.01, 0.01, 0.01, "a"
    ),
    "seed must be numeric, but is a character"
  )
  # ps$active_cell_count[1] <- 2
  expect_error(
    check_run_wavess_inputs(
      inf_pop_size, samp_scheme, c("AG", rep("ATAA", 9)), 1.2, hiv_q_mat,
      3.5e-5, 1.4e-5,
      NULL, 0.99, NULL, 1,
      NULL, 30, 0.01, 90,
      0.001, 0.01, 0.01, 0.01, NULL
    ),
    "All founder sequences must be the same length"
  )
})

test_that("check_extract_seqs_inputs works", {
  expect_no_error(extract_seqs(
    hxb2_cons_founder,
    "B.US.2011.DEMB11US006.KC473833"
  ))
  expect_no_error(extract_seqs(
    hxb2_cons_founder,
    "B.US.2011.DEMB11US006.KC473833", NULL, 3, 10
  ))
  expect_no_error(extract_seqs(
    hxb2_cons_founder,
    "B.US.2011.DEMB11US006.KC473833",
    "B.US.2011.DEMB11US006.KC473833", 3, 10
  ))
  expect_error(
    extract_seqs("not_aln", "B.US.2011.DEMB11US006.KC473833"),
    "aln must be of the"
  )
  expect_error(
    extract_seqs(hxb2_cons_founder, "wrong"),
    "founder_name must be the name of a sequence in aln, but is wrong"
  )
  expect_error(
    extract_seqs(hxb2_cons_founder, "B.US.2011.DEMB11US006.KC473833", "wrong"),
    "ref_name must be the name of a sequence in aln, but is wrong"
  )
  expect_error(
    extract_seqs(hxb2_cons_founder, "B.US.2011.DEMB11US006.KC473833",
      start = "wrong", end = 10
    ),
    "start must be numeric, but is a character"
  )
  expect_error(
    extract_seqs(hxb2_cons_founder, "B.US.2011.DEMB11US006.KC473833",
      start = 1, end = 10000
    ),
    "end must be <= the length of the alignment"
  )
})

test_that("check_slice_aln_inputs works", {
  expect_no_error(slice_aln(hxb2_cons_founder, 3, 10))
  expect_no_error(slice_aln(
    hxb2_cons_founder, 3, 10,
    "B.US.2011.DEMB11US006.KC473833"
  ))
  expect_error(
    slice_aln("not_aln", 3, 10),
    "aln must be of the"
  )
  expect_error(
    slice_aln(hxb2_cons_founder, 3, 10, "wrong"),
    "seqs must be the name of a sequence in aln, but is wrong"
  )
  expect_error(
    slice_aln(hxb2_cons_founder, start = "wrong", end = 10),
    "start must be numeric, but is a character"
  )
  expect_error(
    slice_aln(hxb2_cons_founder, start = 1, end = 10000),
    "end must be <= the length of the alignment"
  )
})
