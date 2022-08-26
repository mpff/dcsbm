test_that("Making an multiple MCMC sweeps works", {
  g1 <- make_full_graph(4, directed=FALSE, loops=FALSE)
  p1 <- c(rep(1,2), rep(2,2))
  expect_error(mcmc_sweep(g1, p1, 2, n.sweeps = 2), NA)

  g4 <- make_empty_graph(4, directed=FALSE)
  p4 <- c(rep(1,2), rep(2,2))
  expect_error(mcmc_sweep(g4, p4, 2, n.sweeps = 2))
})


test_that("MCMC sweeps PPM2 example works", {
  Nsize <- 10
  G <- sample_ppm(Nsize, c = 0.9, k = 10, B = 3)
  p <- c(rep(1, 0.3*Nsize), rep(2, 0.5*Nsize), rep(3, 0.2*Nsize))

  expect_error(mcmc_sweep(G, p, 3, n.sweeps = 2), class = "simpleError")

  G <- simplify(G)
  expect_error(mcmc_sweep(G, p, 3, n.sweeps = 2), NA)
})


test_that("Making an single MCMC sweep works", {
  g1 <- make_full_graph(4, directed=FALSE, loops=FALSE)
  p1 <- c(rep(1,2), rep(2,2))
  expect_error(mcmc_single_sweep(g1, p1, 2), NA)

  g4 <- make_empty_graph(4, directed=FALSE)
  p4 <- c(rep(1,2), rep(2,2))
  expect_error(mcmc_single_sweep(g4, p4, 2))
})


test_that("Single MCMC sweep PPM2 example works", {
  Nsize <- 10
  G <- sample_ppm(Nsize, c = 0.9, k = 10, B = 3)
  p <- c(rep(1, 0.3*Nsize), rep(2, 0.5*Nsize), rep(3, 0.2*Nsize))

  expect_error(mcmc_single_sweep(G, p, 3), class = "simpleError")

  G <- simplify(G)
  expect_error(mcmc_single_sweep(G, p, 3), NA)
})


test_that("Can calculate proposal results", {
  g1 <- make_full_graph(3)
  p1 <- c(1,1,2)
  S1 <- get_entropy(g1, p1)
  el1 <- block_edge_list(g1, p1, 2)

  # Change to equivalent partition.
  res1 <- get_proposal_results(1, 2, S1, g1, p1, el1, eps = 1, beta = 1)
  exp1 <- list("entropy_delta" = 0, "prob_to_accept" = 1, "new_entropy" = S1)
  expect_equal(res1, exp1)

  # Actual change to worse partition.
  res2 <- get_proposal_results(3, 1, S1, g1, p1, el1, eps = 1, beta = 1)
  exp2.dS <- 4.5 * H_binary(6/9) - 2 * H_binary(0.5)
  exp2.pta <- min(exp(-exp2.dS) * (1/8) / (1/2), 1)
  exp2.newS <-  4.5 * H_binary(6/9)
  exp2 <- list(
    "entropy_delta" = exp2.dS,
    "prob_to_accept" = exp2.pta,
    "new_entropy" = exp2.newS)
  expect_equal(res2, exp2)

  # No change.
  res3 <- get_proposal_results(1, 1, S1, g1, p1, el1, eps = 1, beta = 1)
  exp3 <- list("entropy_delta" = 0, "prob_to_accept" = 1, "new_entropy" = S1)
  expect_equal(res3, exp3)
})


test_that("creating block edge list works", {
  g1 <- make_ring(3)
  p1 <- c(1,2,3)
  bel1 <- block_edge_list(g1, p1)
  check <- sapply(seq_along(bel1), function(b){
    bel1[[b]] == as_ids(incident(g1, b))
  })
  expect_true(any(check))

  g2 <- make_full_graph(3)
  p2 <- c(1,1,2)
  bel2 <- block_edge_list(g2, p2)
  expect_equal(bel2[[1]], c(1, 2, 3))
  expect_equal(bel2[[2]], c(2, 3))
})


test_that("updating block edge list works", {
  g1 <- make_ring(3)
  p1 <- c(1,2,3)
  bel1 <- block_edge_list(g1, p1)
  bel2 <- update_block_edge_list(bel1, 1, g1, 2, 1)

  check <- sapply(seq_along(bel1), function(b){
    bel1[[b]] == as_ids(incident(g1, b))
  })
  expect_true(any(check))

})


test_that("swap block works", {
  p1 <- c(1, 1, 2, 2, 3)
  v1 <- 3
  pnb1 <- 1
  p1new <- swap_blocks(p1, v1, pnb1)
  expect_equal(p1new, c(1, 1, 1, 2, 3))
})

