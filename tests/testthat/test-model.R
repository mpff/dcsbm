test_that("Agglomerative merging works", {
  g1 <- make_full_graph(4)
  invisible(capture.output(
    expect_error(sbm(g1), NA)
  ))

  g2 <- sample_ppm2(60, 0.9, 20, 5)
  g2 <- simplify(g2)
  invisible(capture.output(
    m2 <- sbm(g2, n.moves = 1, n.sweeps = 0)
  ))
  n_merges2 <- length(m2$partition)
  expect_equal(max(m2$partition[[n_merges2]]), 1)
})
