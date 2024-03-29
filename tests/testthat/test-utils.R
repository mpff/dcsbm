test_that("Calculating edge counts works", {
  g1 <- make_full_graph(4, directed=FALSE, loops=FALSE)
  p1 <- c(rep(1,2), rep(2,2))
  Ec1 <- block_edge_counts(g1, p1)
  expect_equal(Ec1, matrix(c(2,4,4,2), nrow=2))

  g2 <- make_ring(4, directed=TRUE)
  p2 <- c(rep(1,2), rep(2,2))
  Ec2 <- block_edge_counts(g2, p2)
  expect_equal(Ec2, matrix(c(1,1,1,1), nrow=2))

  g3 <- make_ring(6, directed=TRUE)
  p3 <- c(rep(1,2), rep(2,2), rep(3,2))
  Ec3 <- block_edge_counts(g3, p3)
  expect_equal(Ec3, matrix(c(1,0,1,1,1,0,0,1,1), nrow=3))

  g4 <- make_empty_graph(4)
  p4 <- c(rep(1,2), rep(2,2))
  Ec4 <- block_edge_counts(g4, p4)
  expect_equal(Ec4, diag(0, nrow=2))

  g5 <- make_empty_graph(4)
  p5 <- c(rep(1,2), rep(3,2))
  Ec5 <- block_edge_counts(g5, p5)
  expect_equal(Ec5, diag(0, nrow=3))

  g6 <- make_empty_graph(4)
  p6 <- c(rep(1,2), rep(2,2))
  Ec6 <- block_edge_counts(g6, p6, n.blocks = 3)
  expect_equal(Ec5, diag(0, nrow=3))

  # Manual example
  g7 <- make_full_graph(4)
  p7 <- c(3, 3, 3, 1)
  expect_error(block_edge_counts(g7, p7, n.blocks = 4), NA)
})

test_that("Calculating node counts works", {
  p1 <- c(rep(1,2), rep(2,2))
  n1 <- block_node_counts(p1)
  expect_equal(n1, c(2, 2))

  p2 <- c(rep(1,2), rep(3,2))
  n2 <- block_node_counts(p2)
  expect_equal(n2, c(2,0,2))

  p3 <- c(rep(1,2), rep(2,2))
  n3 <- block_node_counts(p3, n.blocks = 3)
  expect_equal(n3, c(2,2,0))
})

test_that("Calculating effective degree sequence works", {
  g1 <- sample_ppm(30, 0.5, 10, 3)
  p1 <- sample_at_least_once(1:3, 30)
  ds1 <- block_degree_sequence(g1, p1)
  expect_equal(sum(ds1$total), 3)

  g2 <- make_full_graph(4)
  p2 <- c(1,3,4,4)
  ds2 <- block_degree_sequence(g2, p2)
  expect_equal(ds2$total, c(1, 1, 0.5, 0.5))

  g3 <- sample_ppm(30, 0.5, 10, 3, directed = TRUE)
  p3 <- sample_at_least_once(1:3, 30)
  ds3 <- block_degree_sequence(g3, p3)
  expect_equal(sum(ds3$outdegree), 3)
  expect_equal(sum(ds3$indegree), 3)

  g4 <- make_ring(4, directed = TRUE)
  p4 <- c(1,1,2,2)
  ds4 <- block_degree_sequence(g4, p4)
  expect_equal(ds4$outdegree, c(0.5, 0.5, 0.5, 0.5))
  expect_equal(ds4$indegree, c(0.5, 0.5, 0.5, 0.5))
})


test_that("Calculating binary entropy works", {
  H1 <- H_binary(c(0, 0.5, 1))
  expect_equal(H1, c(0, -log(0.5), 0))
})

test_that("Resample has no surprises", {
  x <- 1:10
  expect_equal(length(resample(x[x >  8])), 2)
  expect_equal(length(resample(x[x >  9])), 1)
  expect_equal(length(resample(x[x >  10])), 0)
})

test_that("Sample at least once has no surprises", {
  x <- 1:10

  smpl1 <- sample_at_least_once(x, 10)
  expect_equal(x, sort(smpl1))

  smpl2 <- sample_at_least_once(x, 20, prob = rep(0.1, 10))
  expect_equal(length(smpl2), 20)

  expect_error(sample_at_least_once(x, 5), regexp = "^Not enough*.")
})

test_that("Check partition works", {
  p1 <- c(1, 1, 2, 2, 3, 3)
  p1c <- check_partition(p1)
  expect_equal(p1c, p1)

  p2 <- c(1, 1, 1, 3, 3, 3)
  p2c <- check_partition(p2)
  expect_equal(p2c, c(1, 1, 1, 2, 2, 2))

  p3 <- c(1, 1, 4, 4, 6, 6, 3, 3)
  p3c <- check_partition(p3)
  expect_equal(p3c, c(1, 1, 3, 3, 4, 4, 2, 2))
})
