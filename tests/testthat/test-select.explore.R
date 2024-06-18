library(testthat)
library(BGGM)

test_that("select.explore returns expected structure for two-sided alternative", {
  Y <- matrix(rnorm(100), ncol = 4)
  fit <- explore(Y, progress = FALSE)
  result <- select.explore(fit, alternative = "two.sided")
  
  expect_s3_class(result, "select.explore")
  expect_named(result, c("pcor_mat_zero", "pcor_mat", "pcor_sd_fisher","Adj_10", "Adj_01", "BF_10", "BF_01", "BF_cut", "alternative", "call", "type", "formula", "analytic", "object"))
  expect_true(is.matrix(result$pcor_mat_zero))
  expect_true(is.matrix(result$pcor_mat))
  expect_true(is.matrix(result$Adj_10))
  expect_true(is.matrix(result$Adj_01))
})

test_that("select.explore returns expected structure for greater alternative", {
  Y <- matrix(rnorm(100), ncol = 4)
  fit <- explore(Y, progress = FALSE)
  result <- select.explore(fit, alternative = "greater")
  
  expect_s3_class(result, "select.explore")
  expect_named(result, c("pcor_mat_zero", "pcor_mat", "pcor_sd_fisher", "Adj_20", "Adj_02", "BF_20", "BF_02", "BF_cut", "alternative", "call", "type", "formula", "analytic", "object"))
  expect_true(is.matrix(result$pcor_mat_zero))
  expect_true(is.matrix(result$pcor_mat))
  expect_true(is.matrix(result$Adj_20))
  expect_true(is.matrix(result$Adj_02))
})

test_that("select.explore returns expected structure for exhaustive alternative", {
  Y <- matrix(rnorm(100), ncol = 4)
  fit <- explore(Y, progress = FALSE)
  result <- select.explore(fit, alternative = "exhaustive")
  
  expect_s3_class(result, "select.explore")
  expect_named(result, c("post_prob", "neg_mat", "pos_mat", "null_mat", "alternative", "pcor_mat", "pcor_sd_fisher", "call", "prob", "type", "formula", "analytic", "object"))
  expect_true(is.data.frame(result$post_prob))
  expect_true(is.matrix(result$neg_mat))
  expect_true(is.matrix(result$pos_mat))
  expect_true(is.matrix(result$null_mat))
})
