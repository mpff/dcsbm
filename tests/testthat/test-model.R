test_that("Agglomerative merging works", {
  g1 <- make_full_graph(4)
  invisible(capture.output(
    expect_error(sbm(g1), NA),
    expect_error(sbm(g1, n.sweeps = 1), NA)
  ))

  g2 <- make_full_graph(4, directed = TRUE)
  invisible(capture.output(
    expect_error(sbm(g2), NA),
    expect_error(sbm(g2, n.sweeps = 1), NA)
  ))

  g3 <- sample_ppm(60, 0.9, 20, 5)
  invisible(capture.output(
    m3 <- sbm(g3, n.moves = 1, n.sweeps = 1)
  ))
  n_merges3 <- length(m3$partition)
  expect_equal(max(m3$partition[[n_merges3]]), 1)

  invisible(capture.output(
    expect_error(sbm(g1, degree_correction = TRUE), NA)
  ))

  g4 <- sample_ppm(60,0.9, 10, 3, directed = TRUE, loops = TRUE)
  invisible(capture.output(
    m4 <- sbm(g4, n.moves = 1, n.sweeps = 1)
  ))
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
