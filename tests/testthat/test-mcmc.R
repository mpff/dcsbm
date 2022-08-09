test_that("Making an MCMC sweep works", {
  skip("Added mcmcm_alternate.R")

  g1 <- make_full_graph(4, directed=FALSE, loops=FALSE)
  p1 <- c(rep(1,2), rep(2,2))
  expect_error(mcmc_sweep(g1, p1), NA)

  g2 <- make_ring(4, directed=TRUE)
  p2 <- c(rep(1,2), rep(2,2))
  expect_error(mcmc_sweep(g2, p2), NA)

  g3 <- make_ring(6, directed=TRUE)
  p3 <- c(rep(1,2), rep(2,2), rep(3,2))
  expect_error(mcmc_sweep(g3, p3), NA)

  g4 <- make_empty_graph(4)
  p4 <- c(rep(1,2), rep(2,2))
  expect_equal(mcmc_sweep(g4, p4), c(1,1,1,1))
})

test_that("MCMC sweep PPM example works", {
  skip("Added mcmcm_alternate.R")

  Nsize <- 10
  G <- sample_ppm(Nsize, 0.99, 0.01, c(0.3*Nsize, 0.5*Nsize, 0.2*Nsize))
  p <- c(rep(1, 0.3*Nsize), rep(2, 0.5*Nsize), rep(3, 0.2*Nsize))
  expect_error(mcmc_sweep(G, p), NA)
})

test_that("Random choice probability calculation works", {
  skip("Added mcmcm_alternate.R")

  g1 <- make_full_graph(4, directed=FALSE, loops=FALSE)
  p1 <- c(rep(1,2), rep(2,2))
  E1 <- block_edge_counts(g1, p1)
  pr1 <- random_probs(E1)
  expect_equal(pr1, c(1/3, 1/3))
})

test_that("Transmission probability calculation works", {
  skip("Added mcmcm_alternate.R")

  g1 <- make_full_graph(4, directed=FALSE, loops=FALSE)
  p1 <- c(rep(1,2), rep(2,2))
  E1 <- block_edge_counts(g1, p1)
  Pe1 <- edge_probs(E1)
  expect_equal(Pe1[,1], c(1/3, 2/3))
  expect_equal(Pe1[,2], c(2/3, 1/3))
})
