test_that("package loads and has a version", {
  expect_true(is.character(as.character(utils::packageVersion("mi4p"))))
})

test_that("key exports exist", {
  exports <- c(
    "MVgen","check.conditions","check.design","make.design","test.design",
    "rubin1.one","rubin1.all","rubin2wt.one","rubin2wt.all",
    "rubin2bt.one","rubin2bt.all","proj_matrix"
  )
  for (fn in exports) {
    expect_true(exists(fn, where=asNamespace("mi4p"), inherits=FALSE), info = fn)
    expect_true(is.function(get(fn, envir=asNamespace("mi4p"))), info = fn)
  }
})
