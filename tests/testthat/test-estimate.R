library(testthat)
library(BGGM)

# Set up a test matrix or data frame
test_data <- ptsd[, 1:5]


# Test default call
test_that("estimate with default parameters works", {
  result <- estimate(test_data)
  expect_s3_class(result, c("BGGM", "estimate", "default"))
})

# Test invalid type handling
test_that("estimate with incorrect type returns error", {
  expect_error( estimate(test_data, type = "invalid") )
})

# Add more tests...
