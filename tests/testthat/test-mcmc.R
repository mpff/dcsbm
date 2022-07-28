test_that("Making an MCMC sweep works", {
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
  expect_equal(mcmc_sweep(g4, p4), p4)
})

test_that("MCMC sweep PPM example works", {
  Nsize <- 10
  G <- sample_ppm(Nsize, 0.1, 0.01, c(0.3*Nsize, 0.5*Nsize, 0.2*Nsize))
  p <- c(rep(1, 0.3*Nsize), rep(2, 0.5*Nsize), rep(3, 0.2*Nsize))
  expect_error(mcmc_sweep(G, p), NA)
})
