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
      node("sA", distr = "rnorm", mean = sA.mu, sd = 1) +
      node("sA.cat", distr = "rconst", const = as.integer(sA))
  D <- set.DAG(D, n.test = 10)
  datO <- sim(D, n = nsamp, rndseed = 12345)
}


nsamp <- 10000
datO <- data.table::data.table(get.density.sAdat(nsamp))
ncats <- unique(datO[["sA.cat"]])
print(ncats)
newdata <- datO[1:100, c("W1", "W2", "W3", "sA.cat"), with = FALSE]

test_that("fit_density should work and return object SummariesModel", {
  dens_fit <- fit_density(X = c("W1", "W2", "W3"), Y = "sA.cat", input_data = datO, bin_estimator = speedglmR6$new())
  expect_true(inherits(dens_fit, "SummariesModel"))
})

test_that("predict_probability should return probabilities", {
  dens_fit <- fit_density(X = c("W1", "W2", "W3"), Y = "sA.cat", input_data = datO, bin_estimator = speedglmR6$new())
  preds <- predict_probability(dens_fit, newdata)
  expect_true(is.numeric(preds))
})

test_that("sample_value should return categorical values for categorical outcome", {
  dens_fit <- fit_density(X = c("W1", "W2", "W3"), Y = "sA.cat", input_data = datO, bin_estimator = speedglmR6$new())
  pred_values <- sample_value(dens_fit, newdata)
  expect_true(is.integer(pred_values))
  expect_true(all(pred_values %in% ncats))
  # cbind(newdata, pred_values = pred_values)
})






