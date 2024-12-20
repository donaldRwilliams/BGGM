# test-bma_posterior.R

library(testthat)
library(BGGM)

test_that("bma_posterior returns expected output", {
  # Load example data
  Y <- ptsd[,1:4]
  
  # Fit the model using ggm_search
  fit <- ggm_search(x = Y, prior_prob = 0.9, iter = 5000)
  
  # Generate posterior samples using bma_posterior
  bma_result <- bma_posterior(fit, iter = 100, progress = FALSE)
 
  # Check that the result is a list
  expect_type(bma_result, "list")

  # Check that the list contains 'bma_mean' and 'samples'
  expect_named(bma_result, c("bma_mean", "samples"))

  # Check that 'bma_mean' is a matrix with appropriate dimensions
  expect_true(is.matrix(bma_result$bma_mean))
  expect_equal(dim(bma_result$bma_mean), c(ncol(Y), ncol(Y)))

  # Check that 'samples' is a 3D array with appropriate dimensions
  expect_true(length(dim(bma_result$samples)) == 3)
  expect_equal(dim(bma_result$samples), c(ncol(Y), ncol(Y), 100))
})
