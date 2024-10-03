test_that("calc_tr_summary_stats works", {
  tss <- calc_tr_summary_stats(ape::read.tree(text = "((t2:1,t1:1):1,t3:1);"), c("t1" = 1, "t2" = 2, "t3" = 2))
  expect_equal(tss$stat_name, c("sackin_norm", "int_over_ext", "mean_prop_survived"))
  ### expect_equal(tss$stat_value, c(1.67, 1, 1))
  expect_equal(calc_int_over_ext(ape::read.tree(text = "((t2:1,t1:1):1,t3:4);")), 0.5)
  tr <- ape::read.tree(text = "(t1:1,((t3:1,(t4:1,t2:1):1):1,t5:1):1);")
  expect_equal(
    calc_prop_survived(
      tr,
      c("t1" = 1, "t2" = 2, "t3" = 2, "t4" = 3, "t5" = 3)
    )$prop_survived,
    c(NA, 1, 1)
  )
  expect_equal(
    calc_prop_survived(
      tr,
      c("t1" = 1, "t2" = 3, "t3" = 2, "t4" = 3, "t5" = 2)
    )$prop_survived,
    c(NA, 1, 0.5)
  )
  ### calc_prop_survived(ape::read.tree(text="(t1:1,((t3:1,(t4:1,t2:1):1):1,t5:1):1);"), c('t1'=1,'t2'=1,'t3'=2,'t4'=3,'t5'=3))
  expect_equal(get_clusters(tr, c("t1" = 1, "t2" = 3, "t3" = 2, "t4" = 3, "t5" = 2))$pure_subtree_info$subtr_size, c(1, 1, 2, 1))
  expect_equal(get_clusters(tr, c("t1" = 1, "t2" = 3, "t3" = 2, "t4" = 3, "t5" = 2),
    grps = c("t1" = 1, "t2" = 3, "t3" = 2, "t4" = 3, "t5" = 2)
  )$pure_subtree_info$subtr_size, c(1, 1, 1))
  expect_equal(get_clusters(tr, c("t1" = 1, "t2" = 3, "t3" = 2, "t4" = 3, "t5" = 2),
    grps = c("t1" = 1, "t2" = 3, "t3" = 2, "t4" = 3, "t5" = 4)
  )$pure_subtree_info$subtr_size, c(1, 1, 1, 1))
  tr$node.label <- rep(0.8, 4)
  # note that this isn't necessarily what we want when computing proportion of lineages that survived, but it's the way the function works
  expect_equal(
    get_clusters(tr, c("t1" = 1, "t2" = 3, "t3" = 2, "t4" = 3, "t5" = 2),
      bootstrap = 0.9
    )$pure_subtree_info$subtr_size,
    c(1, 1, 1, 1, 1)
  )
  expect_equal(get_clusters(tr, c("t1" = 1, "t2" = 3, "t3" = 2, "t4" = 3, "t5" = 2),
    pureness = 0.4
  )$pure_subtree_info$subtr_size, c(1, 2, 2))
})
