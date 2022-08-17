test_that("Calculating undirected traditional entropy works", {
  g1 <- make_full_graph(4, directed=FALSE, loops=FALSE)
  p1 <- c(rep(1,2), rep(2,2))
  S1 <- get_entropy(g1, p1)
  expect_equal(S1, 2 * sum(H_binary(c(0.5, 0.5))))

  E2 <- matrix(c(6,0,0,0), nrow = 2)
  n2 <- c(3, 0)
  S2 <- entropy_trad(E2, n2)
  expect_equal(S2, 0.5 * 9 * H_binary(6/9))
})


test_that("Calculating directed traditional entropy works", {
  g1 <- make_full_graph(4, directed=TRUE, loops=FALSE)
  p1 <- c(rep(1,2), rep(2,2))
  S1 <- get_entropy(g1, p1)
  expect_equal(S1, 4 * sum(H_binary(c(0.5, 0.5))))
})


test_that("Calculating oneway degree corrected entropy works", {
  g1 <- make_ring(4)
  p1 <- c(rep(1,2), rep(2,2))
  S1 <- get_entropy(g1, p1, degree_correction = "oneway")
  S1_expected <- - 4 - 4 * log(2) - 0.5 * 4 * 2 * log(0.125)
  expect_equal(S1, S1_expected)
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
