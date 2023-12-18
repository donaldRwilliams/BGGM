library(testthat)
library(BGGM)

test_data <- ptsd[, 1:5]
fit <- estimate(test_data)

# Test summary call
test_that("summary.estimate basic functionality", {
  summary_result <- summary.estimate(fit)
  expect_s3_class(summary_result, c("BGGM", "estimate", "summary_estimate", "summary.estimate"))
  expect_true("Relation" %in% colnames(summary_result$dat_results))
})

# Test col_names = FALSE
test_that("summary.estimate handles col_names properly", {
  summary_result <- summary.estimate(fit, col_names = FALSE)
  expected_relation_format <- grepl("--", summary_result$dat_results$Relation[1])
  expect_true(expected_relation_format)
})

# Add more tests...
