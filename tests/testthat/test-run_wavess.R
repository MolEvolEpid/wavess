test_that("run_wavess works", {
  if (!reticulate::virtualenv_exists(Sys.getenv("WAVESS_PYTHON", unset = "r-wavess"))){
    skip('wavess python module not available for testing')
  }
  hiv_env_flt_2021 <- ape::as.matrix.DNAbin(hiv_env_flt_2021)
  samp_scheme <- define_sampling_scheme(define_growth_curve(gN = 100), sampling_frequency = 50)
  capture.output(probs <- calc_nt_sub_probs(hiv_env_flt_2021[1:3,]), file = nullfile())
  out <- run_wavess(samp_scheme, c('ATCG', 'ATTT'), probs)
  expect_equal(out$counts$generation, c(0,50,100))
  expect_equal(out$counts$active_cell_count, c(1,1998,1999))
  expect_equal(dim(out$counts), c(3,13))
  expect_equal(length(out$seqs), 43)
  no_lat <- run_wavess(samp_scheme, 'ATCG', probs, prob_act_to_lat = 0)
  expect_equal(unique(no_lat$counts$latent_cell_count), 0)
  expect_equal(unique(no_lat$counts$active_turned_latent), 0)
  expect_equal(unique(no_lat$counts$latent_turned_active), 0)
  expect_equal(unique(no_lat$counts$latent_died), 0)
  expect_equal(unique(no_lat$counts$latent_proliferated), 0)
  expect_no_error(run_wavess(samp_scheme, 'ATCG', probs))
  expect_no_error(run_wavess(samp_scheme, c('ATCG', 'ATTT'), probs, ref_seq = 'AAAA'))
  expect_no_error(run_wavess(samp_scheme, c('ATCG', 'ATTT'), probs, ref_seq = 'AAAA'))
  expect_error(run_wavess(samp_scheme, c('ATCG', 'ATTT'), probs, conserved_sites = 0),
               'conserved_sites must be a positive number')
  expect_no_error(run_wavess(samp_scheme, c('ATCG', 'ATTT'), probs, ref_seq = 'AAAA', conserved_sites = 1))
  expect_no_error(run_wavess(samp_scheme, c('ATCG', 'ATTT'), probs, conserved_sites = c(1,2)))
  expect_no_error(run_wavess(samp_scheme, c('ATCG', 'ATTT'), probs,
                             epitope_locations = tibble::tibble(epi_start_nt=1, epi_end_nt=4, max_fitness_cost=0.4)))
})
