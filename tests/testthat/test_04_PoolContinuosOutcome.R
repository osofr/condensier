context("Pooled Continuous Outcome")

library(mockery)
library(condensier)

op <- options(sl3.verbose = TRUE)
# op <- options(condensier.verbose = TRUE)

library("simcausal")
D <- DAG.empty()
D <-
D + node("W1", distr = "rbern", prob = 0.5) +
  node("W2", distr = "rbern", prob = 0.3) +
  node("W3", distr = "rbern", prob = 0.3) +
  node("sA.mu", distr = "rconst", const = (0.98 * W1 + 0.58 * W2 + 0.33 * W3)) +
  node("Z", distr = "rnorm", mean = sA.mu, sd = 1)
D <- set.DAG(D, n.test = 10)
datO <- sim(D, n = 10000, rndseed = 12345)
newdata <- datO[1:5, c("W1", "W2", "W3", "Z")]

test_that("pooled density fitting method works", {
  dens_fit_pooled <- fit_density(
      X = c("W1", "W2", "W3"),
      Y = "Z",
      input_data = datO,
      nbins = 10,
      bin_method = "equal.mass",
      pool = TRUE,
      bin_estimator = speedglmR6$new())
  preds_long <- predict_probability(dens_fit_pooled, newdata)
  #        50%        80%        60%       100%        10%
  # 0.43141041 0.39898270 0.28679683 0.04857514 0.03469636

  ## not yet implemented
  expect_error(sampledY_long <- sample_value(dens_fit_pooled, newdata))
})