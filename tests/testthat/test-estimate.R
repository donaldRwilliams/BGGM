library(testthat)
library(BGGM)

# Set up a test matrix or data frame
test_data <- ptsd[1:20, 1:5]


# Test default call
test_that("estimate with default parameters works", {
  result <- estimate(test_data)
  expect_s3_class(result, c("BGGM", "estimate", "default"))
})

# Test invalid type handling
test_that("estimate with incorrect type returns error", {
  expect_error( estimate(test_data, type = "invalid") )
})

test_that("estimate returns an object of correct class", {
  Y <- matrix(rnorm(100), ncol=5)
  result <- estimate(Y, type = "continuous", iter = 10)
  expect_true("estimate" %in% class(result))
})

test_that("estimate handles different types correctly", {
    Y <- matrix(rnorm(100), ncol=5)

    result_cont <- estimate(Y, type = "continuous", iter = 10)
    expect_true("estimate" %in% class(result_cont))

    result_bin <- estimate(Y > 0, type = "binary", iter = 10)
    expect_true("estimate" %in% class(result_bin))

    ## Expect error for unsupported types
    expect_error(estimate(Y, type = "unsupported", iter = 10))
})

test_that("estimate handles ordinal data without controls", {
    set.seed(123)
    Y <- matrix(sample(1:4, size = 200, replace = TRUE), ncol = 4)
    colnames(Y) <- paste0("item", seq_len(ncol(Y)))

    fit <- estimate(Y, type = "ordinal", iter = 10, progress = FALSE)

    expect_s3_class(fit, c("BGGM", "estimate", "default"))
    expect_identical(fit$type, "ordinal")
    expect_equal(dim(fit$pcor_mat), c(4, 4))
    expect_true(is.array(fit$post_samp$pcors))
})

## Ordinal data check
test_that("estimate handles ordinal data with control variables", {
    set.seed(321)
    Y <- matrix(sample(1:3, size = 180, replace = TRUE), ncol = 3)
    colnames(Y) <- paste0("node", seq_len(ncol(Y)))
    control <- sample(0:1, size = nrow(Y), replace = TRUE)
    dat <- data.frame(Y, control = control)

    fit <- estimate(dat,
        type = "ordinal",
        formula = ~control,
        iter = 10,
        progress = FALSE
    )

    expect_s3_class(fit, c("BGGM", "estimate", "default"))
    expect_identical(fit$type, "ordinal")
    expect_equal(fit$p, ncol(Y))
    expect_equal(dim(fit$pcor_mat), c(ncol(Y), ncol(Y)))
    expect_true(all(c("(Intercept)", "control") %in% colnames(fit$X)))
})


test_that("estimate output structure is as expected", {
    Y <- matrix(rnorm(100), ncol = 5)
    result <- estimate(Y, type = "continuous", iter = 10)
    expect_true(is.list(result))
    expect_true("pcor_mat" %in% names(result))
    expect_true("post_samp" %in% names(result))
})

## Error Handling
test_that("estimate handles errors and missing data", {
  Y <- matrix(c(rnorm(98), NA, NA), ncol=5)
  
  ## Check handling with NA without imputation
  ## No warning coded 
  ## expect_warning(estimate(Y, type = "continuous", impute = FALSE, iter = 10))
  
  # Check handling with NA with imputation
  result_impute <- estimate(Y, type = "continuous", impute = TRUE, iter = 10)
  expect_true("estimate" %in% class(result_impute))

  ## Check handling with NA without imputation
  expect_warning(
    estimate(Y, type = "mixed", impute = FALSE, iter = 10)
  )
})

## Testing Summary and Plot Functions
test_that("summary.estimate returns expected output", {
  Y <- matrix(rnorm(100), ncol=5)
  fit <- estimate(Y, iter = 10)
  summ <- summary(fit)
  
  expect_true("summary_estimate" %in% class(summ))
  expect_true(is.data.frame(summ$dat_results))
})

test_that("plot.summary.estimate works without error", {
  Y <- matrix(rnorm(100), ncol=5)
  fit <- estimate(Y, iter = 10)
  summ <- summary(fit)
  
  expect_silent(plot(summ))
})

