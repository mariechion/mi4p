# Dummy mean and variance functions to avoid heavy dependencies in tests
dummy_mean <- function(ind, peptide, tabdata, metacond) {
  # Return group means for the given peptide and imputation
  tapply(tabdata[peptide, , ind], metacond, mean)
}

dummy_var <- function(ind, peptide, data, metacond) {
  # Return diagonal matrix of within-group variances for the given peptide/imp
  v <- tapply(data[peptide, , ind], metacond, stats::var)
  diag(as.numeric(v), nrow = length(v), ncol = length(v))
}

make_toy_array <- function(P = 3L, S = 6L, M = 4L) {
  # 2 groups of equal size
  metacond <- factor(rep(c("A","B"), each = S/2L), levels = c("A","B"))
  set.seed(123)
  base <- matrix(rnorm(P * S, mean = rep(c(0, 1), each = S/2)), nrow = P, ncol = S)
  # Create M simple imputations by adding small noise
  arr <- array(NA_real_, dim = c(P, S, M))
  for (m in seq_len(M)) arr[,,m] <- base + rnorm(P * S, sd = 0.1)
  list(arr = arr, metacond = metacond)
}

test_that("rubin1.one/all with dummy mean produce expected shapes", {
  toy <- make_toy_array(P = 5, S = 6, M = 3)
  mu1 <- mi4p::rubin1.one(1, toy$arr, funcmean = dummy_mean, metacond = toy$metacond)
  expect_length(mu1, nlevels(toy$metacond))

  mu_all <- mi4p::rubin1.all(toy$arr, toy$metacond, funcmean = dummy_mean, is.parallel = FALSE)
  expect_equal(nrow(mu_all), dim(toy$arr)[1])
  expect_equal(ncol(mu_all), nlevels(toy$metacond))
})

test_that("rubin2wt.one/all with dummy variance produce matrices of right size", {
  toy <- make_toy_array(P = 4, S = 6, M = 5)
  W1 <- mi4p::rubin2wt.one(2, toy$arr, funcvar = dummy_var, metacond = toy$metacond)
  expect_true(is.matrix(W1))
  expect_equal(dim(W1), c(nlevels(toy$metacond), nlevels(toy$metacond)))

  W_all <- mi4p::rubin2wt.all(toy$arr, funcvar = dummy_var, metacond = toy$metacond, is.parallel = FALSE)
  expect_type(W_all, "list")
  expect_length(W_all, dim(toy$arr)[1])
  expect_true(all(vapply(W_all, is.matrix, logical(1))))
})

test_that("rubin2.all combines within/between with supplied functions", {
  toy <- make_toy_array(P = 3, S = 6, M = 4)
  out <- mi4p::rubin2.all(toy$arr, toy$metacond, funcmean = dummy_mean, funcvar = dummy_var, is.parallel = FALSE)
  expect_type(out, "list")
  expect_length(out, dim(toy$arr)[1])
  expect_true(all(vapply(out, function(m) all(dim(m) == c(2,2)), logical(1))))
})
