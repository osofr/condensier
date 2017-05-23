context("sample_value")
# ---------------------------------------------------------------------------------
# Fitting continuous exposure by  binning, conditional on covariates
# Overall exposure g0 (sA) is a mixture of 3 normals,
# individual exposure is normal with mu for each observation being a function of (W1,W2,W3), sd = 1;
# ---------------------------------------------------------------------------------

## helper fun
get.density.sAdat <- function(nsamp = 100000) {
  require(simcausal)
  options(simcausal.verbose = FALSE)
  D <- DAG.empty()
  D <-
  D + node("W1", distr = "rbern", prob = 0.5) +
      node("W2", distr = "rbern", prob = 0.3) +
      node("W3", distr = "rbern", prob = 0.3) +
      node("sA.mu", distr = "rconst", const = (0.98 * W1 + 0.58 * W2 + 0.33 * W3)) +
      node("sA", distr = "rnorm", mean = sA.mu, sd = 1)
  D <- set.DAG(D, n.test = 10)
  datO <- sim(D, n = nsamp, rndseed = 12345)
}

nsamp <- 100
datO <- data.table::data.table(get.density.sAdat(nsamp))
newdata <- datO[1:5, c("W1", "W2", "W3", "sA"), with = FALSE]


test_that("type of the outcome variable in DataStore should be correctly identified as continuous", {
  # nodes <- list(Anodes = "sA", Wnodes = c("W1", "W2", "W3"))
  data_store_obj <- DataStore$new(input_data = data.table::data.table(datO), Y = "sA", X = c("W1", "W2", "W3"))
  expect_true(is.list(data_store_obj$type.sVar))
  expect_true("sA" %in% names(data_store_obj$type.sVar))
  expect_true(data_store_obj$type.sVar[["sA"]] %in% "contin")
})


test_that("the input data in DataStore should be accessible via active binding 'dat.sWsA'", {
  # nodes <- list(Anodes = "sA", Wnodes = c("W1", "W2", "W3"))
  data_store_obj <- DataStore$new(input_data = data.table::data.table(datO), Y = "sA", X = c("W1", "W2", "W3"))
  expect_true(nrow(data_store_obj$dat.sWsA) == nrow(datO))
})


test_that("get.dat.sWsA(rowsubset, covars) should return all rows when missing(rowsubset)", {
  # nodes <- list(Anodes = "sA", Wnodes = c("W1", "W2", "W3"))
  data_store_obj <- DataStore$new(input_data = data.table::data.table(datO), Y = "sA", X = c("W1", "W2", "W3"))
  subset_mat <- data_store_obj$get.dat.sWsA(covars = c("W1","W2"))
  expect_true(nrow(subset_mat)==nrow(datO))
  expect_true(all(names(subset_mat) %in% c("W1","W2")))
})

test_that("Subsetting by get.dat.sWsA(rowsubset, covars) should work for integer and logical rowsubset", {
  # nodes <- list(Anodes = "sA", Wnodes = c("W1", "W2", "W3"))
  data_store_obj <- DataStore$new(input_data = data.table::data.table(datO), Y = "sA", X = c("W1", "W2", "W3"))
  subset_mat <- data_store_obj$get.dat.sWsA(rowsubset = c(TRUE, TRUE), covars = c("W1","W2"))
  expect_error(data_store_obj$get.dat.sWsA(rowsubset = c(1.0, 2.2), covars = c("W1","W2")))
  subset_mat_2 <- data_store_obj$get.dat.sWsA(rowsubset = c(1L, 2L), covars = c("W1","W2"))
  expect_true(all.equal(subset_mat, subset_mat_2))
})

test_that("Subsetting by get.dat.sWsA(rowsubset, covars) should throw an error when covars don't exist", {
  data_store_obj <- DataStore$new(input_data = data.table::data.table(datO), Y = "sA", X = c("W1", "W2", "W3"))
  expect_error(subset_mat_error <- data_store_obj$get.dat.sWsA(rowsubset = TRUE, covars = c("W1","W4", "W5")))
})



