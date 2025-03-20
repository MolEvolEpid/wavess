test_that("run_wavess works", {
  if (!reticulate::virtualenv_exists(Sys.getenv("WAVESS_PYTHON",
    unset = "r-wavess"
  ))) {
    skip("wavess python module not available for testing")
  }
  hiv_env_flt_2022 <- ape::as.matrix.DNAbin(hiv_env_flt_2022)
  inf_pop_size <- define_growth_curve(n_gen = 200)
  samp_scheme <- define_sampling_scheme(sampling_frequency = 50, n_days = 100)
  expect_error(
    run_wavess(inf_pop_size, samp_scheme, c("ATCG", "ATTT")),
    "Initial population size must equal the number of founder sequences"
  )
  expect_error(
    run_wavess(
      define_growth_curve(n_gen = 50),
      define_sampling_scheme(sampling_frequency = 10),
      rep("ATCG", 10)
    ),
    "you requested to sample at a time after max"
  )
  out <- run_wavess(inf_pop_size, samp_scheme, rep("ATCG", 10), generation_time = 1)
  expect_equal(out$counts$generation, c(0, 50, 100))
  expect_equal(out$counts$active_cell_count, c(10, 2000, 2000))
  expect_equal(dim(out$counts), c(3, 13))
  expect_equal(dim(out$seqs), c(60, 4))
  no_lat <- run_wavess(inf_pop_size, samp_scheme, rep("ATCG", 10), act_to_lat = 0)
  expect_equal(unique(no_lat$counts$latent_cell_count), 0)
  expect_equal(unique(no_lat$counts$active_turned_latent), 0)
  expect_equal(unique(no_lat$counts$latent_turned_active), 0)
  expect_equal(unique(no_lat$counts$latent_died), 0)
  expect_equal(unique(no_lat$counts$latent_proliferated), 0)
  expect_no_error(run_wavess(inf_pop_size, samp_scheme, rep("ATCG", 10)))
  expect_no_error(run_wavess(inf_pop_size, samp_scheme, rep("atcg", 10)))
  expect_no_error(run_wavess(inf_pop_size, samp_scheme, rep("ATCG", 10), ref_seq = "aaaa"))
  expect_no_error(run_wavess(inf_pop_size, samp_scheme, rep("ATCG", 10), ref_seq = "AAAA"))
  expect_no_error(run_wavess(inf_pop_size, samp_scheme, rep("ATCG", 10), ref_seq = "AAAA"))
  expect_error(
    run_wavess(inf_pop_size, samp_scheme, rep("ATCG", 10), conserved_sites = c("-1" = "a")),
    "conserved_sites must be a number >= 0, but is -1"
  )
  expect_no_error(run_wavess(inf_pop_size, samp_scheme, rep("ATCG", 10),
    ref_seq = "AAAA", conserved_sites = c("1" = "a")
  ))
  expect_no_error(run_wavess(inf_pop_size, samp_scheme, rep("ATCG", 10),
    conserved_sites = c("1" = "a", "2" = "c")
  ))
  expect_no_error(run_wavess(inf_pop_size, samp_scheme, rep("ATCG", 10),
    epitope_locations = tibble::tibble(
      epi_start_nt = 0, epi_end_nt = 3,
      max_fitness_cost = 0.4
    )
  ))
  # samp_scheme$active_cell_count[1] <- 2
  expect_no_error(run_wavess(inf_pop_size, samp_scheme, c(rep("ATCG", 9), "AAAA")))
  expect_error(
    run_wavess(inf_pop_size, samp_scheme, c(rep("ATCG", 9), "A")),
    "All founder sequences must be the same length"
  )
  set.seed(1234)
  out <- run_wavess(inf_pop_size, samp_scheme, rep("ATCGAT", 10), generation_time = 1,
                    conserved_sites = c('1' = 'A'), ref_seq = "GGGGGG", epitope_locations = tibble(epi_start_nt = 0, epi_end_nt = 3, max_fitness_cost = 0.3))
  expect_equal(all(out$counts$mean_fitness_active != 1), TRUE)
  expect_equal(all(out$counts$mean_conserved_active == 1), TRUE)
  expect_equal(any(out$counts$mean_immune_active != 1), TRUE)
  expect_equal(all(out$counts$mean_replicative_active != 1), TRUE)
  expect_equal(all(out$fitness$overall != 1), TRUE)
  expect_equal(all(out$fitness$conserved == 1), TRUE)
  expect_equal(any(out$fitness$immune != 1), TRUE)
  expect_equal(all(out$fitness$replicative != 1), TRUE)
  })
