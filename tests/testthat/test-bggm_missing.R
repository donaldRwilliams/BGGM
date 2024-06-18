library(testthat)
library(BGGM)

test_that("bggm_missing handles missing data with method 'estimate'", {
    set.seed(123)
    
    # Create a sample dataset with missing values
    Y <- matrix(rnorm(100), ncol = 5)
    Y[sample(length(Y), 20)] <- NA
    
    # Impute missing data
    imp <- mice::mice(Y, m = 5, print = FALSE)
    
    # Test bggm_missing with method 'estimate'
    result <- bggm_missing(imp, method = "estimate", iter = 250, progress = FALSE)
    
    expect_s3_class(result, "estimate")
    expect_true("post_samp" %in% names(result))
    expect_true("pcors" %in% names(result$post_samp))
})

test_that("bggm_missing handles missing data with method 'explore'", {
    set.seed(123)
    
    # Create a sample dataset with missing values
    Y <- matrix(rnorm(100), ncol = 5)
    Y[sample(length(Y), 20)] <- NA
    
    # Impute missing data
    imp <- mice::mice(Y, m = 5, print = FALSE)
    
    # Test bggm_missing with method 'explore'
    result <- bggm_missing(imp, method = "explore", iter = 250, progress = FALSE)
    
    expect_s3_class(result, "explore")
    expect_true("post_samp" %in% names(result))
    expect_true("pcors" %in% names(result$post_samp))
})
