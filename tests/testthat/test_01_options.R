library(mockery)
context("Package Options zzz.R")

test_that("should not allow to set unkown bin_method, default missing value is NA", {
  # expect_error(old_opts <- condensier_options(bin_method = "blah"))

  funmiss <- condensier:::testmisfun()
  expect_true(funmiss(NA))
  expect_true(is.na(condensier:::get.misval()))

  condensier:::set.misval(condensier:::gvars, NaN)
  expect_true(is.nan(condensier:::gvars$misval))
  condensier:::set.misval(condensier:::gvars, NA)
  expect_true(is.na(condensier:::gvars$misval))

  expect_error(condensier:::checkpkgs("blahblah"))

  warns <- condensier:::GetWarningsToSuppress()

  testdat <- data.frame(a = rnorm(5), b = rep("str", 5), stringsAsFactors=TRUE)
  expect_true(condensier:::CheckExistFactors(data = testdat)%in%"b")
})

