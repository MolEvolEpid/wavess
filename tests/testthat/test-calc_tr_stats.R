test_that("calc_tr_summary_stats works", {
  expect_warning(tss <- calc_tr_stats(
    ape::read.tree(text = "((t2:1,t1:1):1,t3:1);"),
    c("t1" = 1, "t2" = 2, "t3" = 2)
  ), "Generation 1 has only one tip, cannot calculate diversity.")
  expect_equal(tss$stat_name, c(
    "mean_leaf_depth",
    "mean_bl", "mean_int_bl", "mean_ext_bl",
    "mean_divergence", "mean_diversity",
    "divergence_slope", "diversity_slope",
    "transition_score"
  ))
  expect_equal(tss$stat_value, c(5 / 3, 1, 1, 1, 1.75, 3, -0.500000000000001, NA, 1 / 2))

  tss <- calc_tr_stats(
    ape::read.tree(text = "((t2:1,t1:1):1,(t3:2,t4:2):1);"),
    c("t1" = 1, "t2" = 1, "t3" = 2, "t4" = 2)
  )
  expect_equal(tss$stat_value, c(2, 4 / 3, 1, 1.5, 2.5, 3, 1, 2, 1 / 2))

  expect_error(
    calc_tr_stats(
      ape::read.tree(text = "((t2:1,t1:1):1,t3:1);"),
      c(1, "t2" = 2, "t3" = 2, "t4" = 3, "t5" = 3)
    ),
    "timepoints must be a vector named by tr tip labels, and must "
  )
})

test_that("calc_tr_dists works", {
  expect_equal(
    calc_tr_dists(ape::read.tree(text = "((t2:1,t1:1):1,t3:1);")),
    structure(list(mean_rtt = (2 * 2 + 1) / 3, mean_ttt = (2 + 3 + 3) / 3), class = c(
      "tbl_df",
      "tbl", "data.frame"
    ), row.names = c(NA, -1L))
  )
})
