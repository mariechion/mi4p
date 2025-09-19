test_that("proj_matrix returns numeric vector matching input length", {
  # Build minimal metadata (sTab) and design with two conditions
  sTab <- data.frame(
    Sample.name = paste0("X", 1:6),
    Condition   = factor(rep(c("A","B"), each = 3), levels = c("A","B")),
    Bio.Rep     = 1:6,
    stringsAsFactors = FALSE
  )
  # Two simple 2x2 positive diagonal matrices
  V <- list(diag(c(1,2)), diag(c(0.5, 3)))
  out <- mi4p::proj_matrix(V, sTab)
  expect_type(out, "double")
  expect_length(out, length(V))
  expect_true(all(is.finite(out)))
})
