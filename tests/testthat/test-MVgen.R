test_that("MVgen amputation roughly respects requested missingness", {
  set.seed(42)
  X <- matrix(rnorm(100), nrow = 10, ncol = 10)
  X_amp <- mi4p::MVgen(X, prop_NA = 0.2)
  prop <- mean(is.na(X_amp))
  # Because selection has replacement, realized missingness <= requested
  expect_gt(prop, 0)
  expect_lte(prop, 0.2 + 1e-8)
  expect_equal(dim(X_amp), dim(X))
})
