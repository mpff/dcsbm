test_that("Agglomerative merging works", {
  g1 <- make_full_graph(4)
  invisible(capture.output(
    expect_error(sbm(g1), NA)
  ))

  g2 <- make_full_graph(4, directed = TRUE)
  invisible(capture.output(
    expect_error(sbm(g2), NA)
  ))

  g3 <- sample_ppm2(60, 0.9, 20, 5)
  g3 <- simplify(g3)
  invisible(capture.output(
    m3 <- sbm(g3, n.moves = 1, n.sweeps = 0)
  ))
  n_merges3 <- length(m3$partition)
  expect_equal(max(m3$partition[[n_merges3]]), 1)
})


test_that("Can create block sequences", {
  bs1 <- block_sequence(5)
  expect_equal(bs1, c(5,3,2,1))

  bs2 <- block_sequence(5, sigma = 4)
  expect_equal(bs2, c(5,2,1))

  bs3 <- block_sequence(10, 6, sigma = 4)
  expect_equal(bs3, c(10,7,6))

  bs4 <- block_sequence(10, 6)
  expect_equal(bs4, c(10,8,7,6))
})
