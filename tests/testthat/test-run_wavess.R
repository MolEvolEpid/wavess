test_that("run_wavess works", {
  if (!reticulate::virtualenv_exists(Sys.getenv("WAVESS_PYTHON",
    unset = "r-wavess"
  ))) {
    skip("wavess python module not available for testing")
  }
  hiv_env_flt_2022 <- ape::as.matrix.DNAbin(hiv_env_flt_2022)
  samp_scheme <- define_sampling_scheme(define_growth_curve(n_gen = 100),
    sampling_frequency = 50
  )
  expect_error(
    run_wavess(samp_scheme, c("ATCG", "ATTT")),
    "Initial population size must equal the number of founder sequences"
  )
  out <- run_wavess(samp_scheme, "ATCG")
  expect_equal(out$counts$generation, c(0, 50, 100))
  expect_equal(out$counts$active_cell_count, c(1, 2000, 2000))
  expect_equal(dim(out$counts), c(3, 13))
  expect_equal(dim(out$seqs), c(42, 4))
  no_lat <- run_wavess(samp_scheme, "ATCG", prob_act_to_lat = 0)
  expect_equal(unique(no_lat$counts$latent_cell_count), 0)
  expect_equal(unique(no_lat$counts$active_turned_latent), 0)
  expect_equal(unique(no_lat$counts$latent_turned_active), 0)
  expect_equal(unique(no_lat$counts$latent_died), 0)
  expect_equal(unique(no_lat$counts$latent_proliferated), 0)
  expect_no_error(run_wavess(samp_scheme, "ATCG"))
  expect_no_error(run_wavess(samp_scheme, "atcg"))
  expect_no_error(run_wavess(samp_scheme, "atcg", ref_seq = "aaaa"))
  expect_no_error(run_wavess(samp_scheme, "ATCG", ref_seq = "AAAA"))
  expect_no_error(run_wavess(samp_scheme, "ATCG", ref_seq = "AAAA"))
  expect_error(
    run_wavess(samp_scheme, "ATCG", conserved_sites = c("-1" = "a")),
    "conserved_sites must be a number >= 0, but is -1"
  )
  expect_no_error(run_wavess(samp_scheme, "ATCG",
    ref_seq = "AAAA", conserved_sites = c("1" = "a")
  ))
  expect_no_error(run_wavess(samp_scheme, "ATCG",
    conserved_sites = c("1" = "a", "2" = "c")
  ))
  expect_no_error(run_wavess(samp_scheme, "ATCG",
    epitope_locations = tibble::tibble(
      epi_start_nt = 0, epi_end_nt = 3,
      max_fitness_cost = 0.4
    )
  ))
  samp_scheme$active_cell_count[1] <- 2
  expect_no_error(run_wavess(samp_scheme, c("ATCG", "AAAA")))
  expect_error(
    run_wavess(samp_scheme, c("ATCG", "A")),
    "All founder sequences must be the same length"
  )
})
