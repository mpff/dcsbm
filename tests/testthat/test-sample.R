test_that("Generating parted partition models works", {

  # See https://github.com/igraph/rigraph/blob/dev/tests/testthat/test_sbm.game.R
  bs <- c(4,6)
  g1 <- sample_ppm(10, p=1, q=1, block.sizes=bs,
                   directed=FALSE, loops=FALSE)
  expect_true(graph.isomorphic(g1, make_full_graph(10, directed=FALSE, loops=FALSE)))

  g2 <- sample_ppm(10, p=1, q=1, block.sizes=bs,
                   directed=FALSE, loops=TRUE)
  g2x <- make_full_graph(10, directed=FALSE, loops=TRUE)
  expect_equal(g2[sparse=FALSE], g2x[sparse=FALSE])

  g3 <- sample_ppm(10, p=1, q=1, block.sizes=bs,
                   directed=TRUE, loops=FALSE)
  g3x <- make_full_graph(10, directed=TRUE, loops=FALSE)
  expect_equal(g3[sparse=FALSE], g3x[sparse=FALSE])

  g4 <- sample_ppm(10, p=1, q=1, block.sizes=bs,
                   directed=TRUE, loops=TRUE)
  g4x <- make_full_graph(10, directed=TRUE, loops=TRUE)
  expect_equal(g4[sparse=FALSE], g4x[sparse=FALSE])
})


test_that("Generating weighted parted partition models works", {

  bs <- c(4,6)
  g1 <- sample_ppm(10, p=1, q=1, block.sizes=bs,
                   directed=FALSE, loops=FALSE, n.trials = 3)
  g1x <- make_full_graph(10, directed=FALSE, loops=FALSE)
  E(g1x)$weight <- 3
  expect_equal(g1[sparse=FALSE], g1x[sparse=FALSE])

  g2 <- sample_ppm(10, p=1, q=1, block.sizes=bs,
                   directed=FALSE, loops=TRUE, n.trials = 3)
  g2x <- make_full_graph(10, directed=FALSE, loops=TRUE)
  E(g2x)$weight <- 3
  expect_equal(g2[sparse=FALSE], g2x[sparse=FALSE])

  g3 <- sample_ppm(10, p=1, q=1, block.sizes=bs,
                   directed=TRUE, loops=FALSE, n.trials = 3)
  g3x <- make_full_graph(10, directed=TRUE, loops=FALSE)
  E(g3x)$weight <- 3
  expect_equal(g3[sparse=FALSE], g3x[sparse=FALSE])

  g4 <- sample_ppm(10, p=1, q=1, block.sizes=bs,
                   directed=TRUE, loops=TRUE, n.trials = 3)
  g4x <- make_full_graph(10, directed=TRUE, loops=TRUE)
  E(g4x)$weight <- 3
  expect_equal(g4[sparse=FALSE], g4x[sparse=FALSE])
})
