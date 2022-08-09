test_that("Calculating undirected traditional entropy works", {
  g1 <- make_full_graph(4, directed=FALSE, loops=FALSE)
  p1 <- c(rep(1,2), rep(2,2))
  S1 <- get_entropy(g1, p1)
  expect_equal(S1, 2 * sum(H_binary(c(0.5, 0.5))))
})

test_that("Calculating entropy does not fail, when groups missing in partition", {
  g1 <- make_ring(6, directed=FALSE)
  p1 <- c(rep(1,3), rep(3,3))
  expect_error(get_entropy(g1, p1), NA)

  g2 <- make_ring(6, directed=FALSE)
  p2 <- c(rep(1,3), rep(2,3))
  expect_error(get_entropy(g2, p2), NA)

  S1 <- get_entropy(g1, p1)
  S2 <- get_entropy(g2, p2)
  expect_equal(S1, S2)

  g3 <- make_empty_graph(6, directed=FALSE)
  g3 <- add.edges(g3, c(2,3))
  p3 <- c(rep(1,3), rep(2,3))
  expect_error(get_entropy(g3, p3), NA)
})
