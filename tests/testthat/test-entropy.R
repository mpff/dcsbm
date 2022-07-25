test_that("Calculating undirected traditional entropy works", {
  g1 <- make_full_graph(4, directed=FALSE, loops=FALSE)
  p1 <- c(rep(1,2), rep(2,2))
  S1 <- get_entropy(g1, p1)
  expect_equal(S1, 2 * sum(H_binary(c(0.5, 0.5))))
})

test_that("Calculating entropy does not fail, when groups missing in partition", {
  g1 <- make_ring(6, directed=TRUE)
  p1 <- c(rep(1,3), rep(3,3))
  expect_error(get_entropy(g1, p1), NULL)
})
