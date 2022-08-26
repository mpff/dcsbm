test_that("Generating parted partition models (Piexoto) works", {

  g1 <- sample_ppm(10, c = 0.5, k = 10, B = 2)
  cg1 <- contract(g1, c(1,1,1,1,1,2,2,2,2,2))
  adj1 <- as_adj(cg1, sparse = FALSE)
  diag(adj1) <- 2 * diag(adj1)

  expect_equal(degree(cg1), c(100, 100))
  expect_equal(adj1, matrix(50, nrow=2, ncol=2))

  g2 <- sample_ppm(10, c=0.5, k=10, B=2, directed=TRUE)
  expect_true(is.directed(g2))
})

test_that("Generating parted partition models with degree variability works", {

  g2 <- sample_dcppm(10, c=0.5, k=10, B=2, k_coef=1, directed=TRUE)
  expect_true(is.directed(g2))

  g3 <- sample_dcppm(10, c=0.5, k=20, B=2, k_coef=1, loops=TRUE)
  expect_true(any(is.loop(g3)))

  g4 <- sample_dcppm(10, c = 0.5, k = 10, B = 2, k_coef=1, directed=TRUE, loops=TRUE)
  expect_true(is.directed(g4) & any(is.loop(g4)))
})
