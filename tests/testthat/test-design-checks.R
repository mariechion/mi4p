test_that("check.conditions validates basic good/bad cases", {
  good <- factor(rep(c("A","B"), each = 3))
  bad_one_level <- factor(rep("A", 6))
  bad_missing <- factor(c("A","A","B","B",NA,"B"))

  res_good <- mi4p::check.conditions(good)
  expect_true(is.list(res_good))
  expect_true(isTRUE(res_good$valid))

  res_bad1 <- mi4p::check.conditions(bad_one_level)
  expect_false(isTRUE(res_bad1$valid))
  expect_match(res_bad1$warn, "at least two conditions", ignore.case = TRUE)

  res_bad2 <- mi4p::check.conditions(bad_missing)
  expect_false(isTRUE(res_bad2$valid))
})

test_that("check.design and make.design work for level-1 design (Condition + Bio.Rep)", {
  # Construct a simple sample table with two conditions and one replicate column
  sTab <- data.frame(
    Sample.name = paste0("X", 1:10),
    Condition   = factor(rep(c("A","B"), each = 5), levels = c("A","B")),
    Bio.Rep     = 1:10,
    stringsAsFactors = FALSE
  )

  chk <- mi4p::check.design(sTab)
  expect_true(is.list(chk))
  expect_true(isTRUE(chk$valid))

  # Build design and check its shape/content
  D <- mi4p::make.design(sTab)
  expect_true(is.matrix(D))
  expect_equal(nrow(D), nrow(sTab))
  expect_equal(ncol(D), nlevels(sTab$Condition))
  expect_true(all(colSums(D) == as.integer(table(sTab$Condition))))
  expect_true(all(D %in% c(0,1)))
})

test_that("check.design catches missing Bio.Rep entries", {
  sTab_bad <- data.frame(
    Sample.name = paste0("X", 1:6),
    Condition   = factor(rep(c("A","B"), each = 3)),
    Bio.Rep     = c(1,2,3,4,NA,6),
    stringsAsFactors = FALSE
  )
  res <- mi4p::check.design(sTab_bad)
  expect_false(isTRUE(res$valid))
  expect_match(res$warn, "not full filled", ignore.case = TRUE)
})
