test_that("calc_tr_summary_stats works", {
  tss <- calc_tr_summary_stats(
    ape::read.tree(text = "((t2:1,t1:1):1,t3:1);"),
    c("t1" = 1, "t2" = 2, "t3" = 2)
  )
  expect_equal(tss$stat_name, c(
    "sackin", "int_bl", "ext_bl", "int_over_ext",
    "parsimony_score"
  ))
  expect_equal(tss$stat_value, c(5, 1, 1, 1, 1))
  expect_equal(
    calc_int_over_ext(ape::read.tree(text = "((t2:1,t1:1):1,t3:4);")),
    c(1, 2, 0.5)
  )
  expect_equal(
    calc_parsimony(
      ape::read.tree(text = "(t1:1,((t3:1,(t4:1,t2:1):1):1,t5:1):1);"),
      c("t1" = 1, "t2" = 2, "t3" = 2, "t4" = 3, "t5" = 3)
    ),
    3
  )
  expect_error(
    calc_parsimony(
      ape::read.tree(text = "((t2:1,t1:1):1,t3:1);"),
      c(1, "t2" = 2, "t3" = 2, "t4" = 3, "t5" = 3)
    ),
    "timepoints must be a vector named by tr tip labels, and must "
  )
})
