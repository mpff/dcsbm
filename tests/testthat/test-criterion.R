test_that("Can calculate MDL criterion", {
  # Undirected
  g1 <- make_full_graph(4, directed=FALSE, loops=FALSE)
  p1 <- c(rep(1,2), rep(2,2))
  c1 <- calculate_mdl(g1, p1)
  expect_equal(c1, 2 * sum(H_binary(c(0.5, 0.5))) + 6 * h_func(2 * 3 / (2 * 6)) + 4 * log(2))

  # Directed
  g2 <- make_full_graph(4, directed=TRUE, loops=FALSE)
  p2 <- c(rep(1,2), rep(2,2))
  c2 <- calculate_mdl(g2, p2)
  expect_equal(c2, 4 * sum(H_binary(c(0.5, 0.5))) + 12 * h_func(2 * 2 / 12) + 4 * log(2))

  # Weighted
  g3 <- make_full_graph(2, directed=FALSE)
  g3 <- add.edges(g3, c(1,2))
  p3 <- c(1,2)
  c3 <- calculate_mdl(g3, p3)
  expect_equal(c3, 3 * H_binary(1/3) + 2 * h_func(2 * 3 / 2 / 2) + 2 * log(2))
})


test_that("Can calculate degree corrected MDL", {
  # Undirected
  g1 <- make_ring(4, directed=FALSE)
  p1 <- c(rep(1,2), rep(2,2))
  c1 <- calculate_mdl(g1, p1, degree_correction = TRUE)
  expect_equal(c1,- 4 - 4 * log(2) - 0.5 * 4 * 2 * log(0.125) +
                 4 * h_func(2 * 3 / (2 * 4)) + 4 * log(2) -
                 4 * 1 * log(1))

  # Weighted
  g2 <- make_full_graph(2, directed=FALSE)
  g2 <- add.edges(g2, c(1,2))
  p2 <- c(1,2)
  c2 <- calculate_mdl(g2, p2, degree_correction = TRUE)
  expect_equal(c2, - 2 - 2 * log(2) - .5 * 2 * 2 * log(1/2) +
                 2 * h_func(2 * 3 / (2 * 2)) + 2 * log(2) -
                 2 * 1 * log(1))

  # Weighted + Directed
  g3 <- make_full_graph(2, directed = TRUE)
  g3 <- add.edges(g3, c(1,2))
  p3 <- c(1,2)
  c3 <- calculate_mdl(g3, p3, degree_correction = TRUE)
  expect_equal(c3, - 3 - 2 * log(2) - 2 * log(2) - log(0.25) +
                 3 * h_func(2 * 2 / 3) + 2 * log(2) -
                 2 * 4 * 0.5 * log(0.5))
})
