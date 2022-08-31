test_that("Calculating traditional entropy works", {

  # Undirected
  g1 <- make_full_graph(4, directed=FALSE, loops=FALSE)
  p1 <- c(rep(1,2), rep(2,2))
  S1 <- get_entropy(g1, p1)
  expect_equal(S1, 2 * sum(H_binary(c(0.5, 0.5))))

  E3 <- matrix(c(6,0,0,0), nrow = 2)
  n3 <- c(3, 0)
  S3 <- entropy_trad(E3, n3)
  expect_equal(S3, 0.5 * 9 * H_binary(6/9))

  # Directed
  g2 <- make_full_graph(4, directed=TRUE, loops=FALSE)
  p2 <- c(rep(1,2), rep(2,2))
  S2 <- get_entropy(g2, p2)
  expect_equal(S2, 4 * sum(H_binary(c(0.5, 0.5))))

  E4 <- matrix(c(6,0,0,0), nrow = 2)
  n4 <- c(3, 0)
  S4 <- entropy_trad(E4, n4, directed = TRUE)
  expect_equal(S4,  9 * H_binary(6/9))

  # Weighted
  g5 <- make_full_graph(2, directed=FALSE)
  g5 <- add.edges(g5, c(1,2))
  p5 <- c(1,2)
  S5 <- get_entropy(g5, p5)
  expect_equal(S5, 3 * H_binary(1/3))

  E6 <- matrix(c(6,0,0,0), nrow = 2)
  n6 <- c(3, 0)
  S6 <- entropy_trad(E6, n6, directed = FALSE, simple = FALSE)
  expect_equal(S6,  .5 * (9 + 6) * H_binary(9/(9+6)))

  # Weighted + Directed
  S7 <- entropy_trad(E6, n6, directed = TRUE, simple = FALSE)
  expect_equal(S7,  (9 + 6) * H_binary(9/(9+6)))
})


test_that("Calculating degree corrected entropy works", {
  g1 <- make_ring(4, directed=FALSE)
  p1 <- c(rep(1,2), rep(2,2))
  S1 <- get_entropy(g1, p1, degree_correction = TRUE)
  expect_equal(S1,- 4 - 4 * log(2) - 0.5 * 4 * 2 * log(0.125))

  g2 <- make_ring(3, directed=TRUE)
  g2 <- add.edges(g2, c(1,3))
  p2 <- c(1,1,2)
  S2 <- get_entropy(g2, p2, degree_correction = TRUE)
  expect_equal(S2,- 4 - 2 * log(2) - (1 * log(1/2/3) + 2 * log(2/2/1) + 1 * log(1/2/3)) )

  E3 <- matrix(c(6,0,0,0), nrow = 2)
  d3 <- list("total" = c(2, 2, 2))
  S3 <- entropy_corrected(E3, d3)
  expect_equal(S3, - 3 - 3 * log(2) - 0.5 * 6 * log(1/6))

  E4 <- matrix(c(6,0,0,0), nrow = 2)
  d4 <- list("in" = c(1, 1, 1), "out" = c(1,1,1))
  S4 <- entropy_corrected(E4, d4, directed = TRUE)
  expect_equal(S4, - 6 - 6 * log(1/6))

  # Weighted
  g5 <- make_full_graph(2)
  g5 <- add.edges(g5, c(1,2))
  p5 <- c(1,2)
  S5 <- get_entropy(g5, p5, degree_correction = TRUE)
  expect_equal(S5, - 2 - 2 * log(2) - .5 * 2 * 2 * log(1/2))

  # Weighted + Directed
  g6 <- make_full_graph(2, directed = TRUE)
  g6 <- add.edges(g6, c(1,2))
  p6 <- c(1,2)
  S6 <- get_entropy(g6, p6, degree_correction = TRUE)
  expect_equal(S6, - 3 - 2 * log(2) - 2 * log(2) - log(0.25))
})


test_that("Calculating entropy does not fail, when groups missing in partition", {
  g <- make_ring(6)

  p1 <- c(rep(1,3), rep(3,3))
  expect_error(get_entropy(g, p1), NA)

  p2 <- c(rep(1,3), rep(2,3))
  expect_error(get_entropy(g, p2), NA)

  S1 <- get_entropy(g, p1)
  S2 <- get_entropy(g, p2)
  expect_equal(S1, S2)

  g3 <- add.edges(g, c(2,3))
  p3 <- c(rep(1,3), rep(2,3))
  expect_error(get_entropy(g3, p3), NA)
})

