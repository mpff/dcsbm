test_that("Generating graph based on parted partition model works", {

  pt <- sample(1:2, 10, replace = T)
  expect_error(sample_ppa(pt, 0.3, 0.01), NULL)

})
