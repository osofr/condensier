context("fit_density")
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


nsamp <- 10000
datO <- data.table::data.table(get.density.sAdat(nsamp))
newdata <- datO[1:5, c("W1", "W2", "W3", "sA"), with = FALSE]


test_that("fit_density should work and return object SummariesModel", {
  dens_fit <- fit_density(X = c("W1", "W2", "W3"), Y = "sA", input_data = datO, nbins = 20, bin_estimator = speedglmR6$new())
  expect_true(inherits(dens_fit, "SummariesModel"))
})


test_that("predict_probability should return probabilities", {
  dens_fit <- fit_density(X = c("W1", "W2", "W3"), Y = "sA", input_data = datO, nbins = 20, bin_estimator = speedglmR6$new())
  preds <- predict_probability(dens_fit, newdata)
  expect_true(is.numeric(preds))
})


test_that("new DataStore auto_typing FALSE should not define the variable types", {
  newdata_obj <- DataStore$new(input_data = newdata, Y = "sA", X = c("W1", "W2", "W3"), auto_typing = FALSE)
  expect_true(is.null(newdata_obj$type.sVar))
})


test_that("call to predict should return SummariesModel object", {
  newdata_obj <- DataStore$new(input_data = newdata, Y = "sA", X = c("W1", "W2", "W3"), auto_typing = FALSE)
  dens_fit <- fit_density(X = c("W1", "W2", "W3"), Y = "sA", input_data = datO, nbins = 20, bin_estimator = speedglmR6$new())
  preds <- predict_probability(dens_fit, newdata)
  res <- dens_fit$predict(newdata = newdata_obj)
  expect_true(is(res, "SummariesModel"))
})


test_that("calls to predict_probability and predictAeqa should return identical things", {
  newdata_obj <- DataStore$new(input_data = newdata, Y = "sA", X = c("W1", "W2", "W3"), auto_typing = FALSE)
  dens_fit <- fit_density(X = c("W1", "W2", "W3"), Y = "sA", input_data = datO, nbins = 20, bin_estimator = speedglmR6$new())
  preds <- predict_probability(dens_fit, newdata)
  preds_test <- dens_fit$predictAeqa(newdata = newdata_obj)
  expect_true(is.numeric(preds_test))
  expect_true(all.equal(preds, preds_test))
})






