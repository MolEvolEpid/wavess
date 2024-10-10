test_that("sample_epitopes works", {
  set.seed(20240812)
  expect_message(epitopes <- sample_epitopes(
    get_epitope_frequencies(env_features$position)
  ))
  expect_equal(
    colnames(epitopes),
    c("epi_start_nt", "epi_end_nt", "max_fitness_cost")
  )
  expect_equal(sum(is.na(epitopes$epi_start_nt)), 0)
  expect_equal(sum(is.na(epitopes$epi_end_nt)), 0)
  expect_equal(epitopes$epi_start_nt,
               c(829, 238, 958, 1357, 394, 922, 502, 172, 541, 340))
  expect_equal(epitopes$epi_end_nt,
               c(829, 238, 958, 1357, 394, 922, 502, 172, 541, 340) + 30-1)
  expect_equal(
    epitopes$max_fitness_cost,
    c(0.040, 0.080, 0.120, 0.160, 0.200, 0.240, 0.280, 0.320, 0.360, 0.400)
  )
  expect_message(sample_epitopes(get_epitope_frequencies(env_features$position),
    cost_type = "random"
  ))
  expect_message(sample_epitopes(get_epitope_frequencies(env_features$position),
    end_aa_pos = 1000
  ))
  expect_error(
    sample_epitopes(get_epitope_frequencies(env_features$position),
      max_resamples = 1
    ),
    "Too many resamples required."
  )
  expect_message(sample_epitopes(get_epitope_frequencies(env_features$position),
    end_aa_pos = (7758 - 6225) / 3,
    ref_founder_map = map_ref_founder(
      slice_aln(hxb2_cons_founder, 6225, 7787),
      labels(hxb2_cons_founder)[1],
      labels(hxb2_cons_founder)[2]
    )
  ))
  expect_message(sample_epitopes(get_epitope_frequencies(env_features$position),
    aa_epitope_length = 11
  ))
})

test_that("get_epitope_frequencies works", {
  expect_no_error(epi_freqs <- get_epitope_frequencies(c(1, 3, 3)))
  expect_equal(epi_freqs$aa_position, 1:3)
  expect_equal(epi_freqs$n_features, c(1, 0, 2))
  expect_equal(epi_freqs$epitope_probability, c(1 / 3, 0, 2 / 3))
})

test_that("reindex_epitopes works", {
  expect_equal(
    reindex_epitopes(1, 2, 0.1,
      tibble::tibble(
        ref_pos = 1:10,
        founder_pos = 1:10
      )
    ),
    tibble::tibble(
      epi_start_nt = 1,
      epi_end_nt = 6,
      max_fitness_cost = 0.1
    )
  )
  expect_equal(
    reindex_epitopes(2, 3, 0.1,
                     tibble::tibble(
                       ref_pos = 1:10,
                       founder_pos = c(1:3,NA,NA,4:8)
                     )
    ),
    tibble::tibble(
      epi_start_nt = 1,
      epi_end_nt = 9,
      max_fitness_cost = 0.1
    )
  )
  expect_equal(
    reindex_epitopes(3, 1, 0.1,
                     tibble::tibble(
                       ref_pos = 1:10,
                       founder_pos = c(1:3,NA,NA,4:8)
                     )
    ),
    tibble::tibble(
      epi_start_nt = 4,
      epi_end_nt = 6,
      max_fitness_cost = 0.1
    )
  )
  expect_error(reindex_epitopes(5, 1, 0.1,
                   tibble::tibble(
                     ref_pos = 1:10,
                     founder_pos = c(1:3,NA,NA,4:8)
                   )),
               "Not all reference epitope start and end positions are in ")
})
