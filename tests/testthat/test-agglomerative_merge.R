test_that("Can collapse one block", {
  g1 <- make_full_graph(3)
  p1 <- 1:3
  expect_error(cb1$new_partition, NA)

  p2 <- rep(1,3)
  expect_warning(collapse_step(g1, p2),  regexp = "^Cannot merge*.")
})
