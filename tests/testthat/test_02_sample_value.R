context("sample_value")
# ---------------------------------------------------------------------------------
# Fitting continuous exposure by  binning, conditional on covariates
# Overall exposure g0 (sA) is a mixture of 3 normals,
# individual exposure is normal with mu for each observation being a function of (W1,W2,W3), sd = 1;
# ---------------------------------------------------------------------------------
`%+%` <- function(a, b) paste0(a, b)
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
allNA <- function(x) all(is.na(x))

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

test_that("sample_value should work as expected", {
  dens_fit <- fit_density(X = c("W1", "W2", "W3"), Y = "sA", input_data = datO, nbins = 20, bin_estimator = speedglmR6$new())
  set.seed(123456)
  ## no outcome in the input dataset:
  sampledY_1 <- sample_value(dens_fit, newdata[, c("W1", "W2", "W3"), with = FALSE])
  set.seed(123456)
  sampledY_2 <- sample_value(dens_fit, newdata)
  expect_true(all.equal(sampledY_1, sampledY_2))
  set.seed(123456)
  newdata_obj <- DataStore$new(input_data = newdata, Y = "sA", X = c("W1", "W2", "W3"), auto_typing = FALSE)
  sampledY_test <- dens_fit$sampleA(newdata = newdata_obj)
  expect_true(all.equal(sampledY_2, sampledY_test))
})

