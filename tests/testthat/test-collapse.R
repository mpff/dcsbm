test_that("Can collapse one block", {
  g1 <- make_full_graph(3)
  p1 <- 1:3
  invisible(capture.output(
    expect_error(collapse_step(g1, p1), NA)
  ))

  g2 <- make_full_graph(3)
  p2 <- rep(1,3)
  invisible(capture.output(
    expect_warning(collapse_step(g2, p2), regexp = "^Cannot merge*.")
  ))

  g3 <- make_full_graph(4)
  p3 <- 1:4
  invisible(capture.output(
    cb3 <- collapse_step(g3, p3, n.merges = 2)
  ))
  expect_equal(max(cb3$new_partition), max(p3) - 2)
})
