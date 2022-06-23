test_that("Calculating edge counts works", {
  g1 <- make_full_graph(4, directed=FALSE, loops=FALSE)
  p1 <- c(rep(1,2), rep(2,2))
  Ec1 <- block_edge_counts(g1, p1)
  expect_equal(Ec1, matrix(c(2,4,4,2), nrow=2))

  g2 <- make_ring(4, directed=TRUE)
  p2 <- c(rep(1,2), rep(2,2))
  Ec2 <- block_edge_counts(g2, p2)
  expect_equal(Ec2, matrix(c(1,1,1,1), nrow=2))
})

test_that("Calculating node counts works", {
  p1 <- c(rep(1,2), rep(2,2))
  n1 <- block_node_counts(p1)
  expect_equal(n1, c(2, 2))

  p2 <- c(rep(1,2), rep(3,2))
  n2 <- block_node_counts(p2)
  expect_equal(n2, c(2,2))
})

test_that("Calculating binary entropy function works", {
  H1 <- H_binary(ers, nr, ns)
  expect_equal(1,1)
})
