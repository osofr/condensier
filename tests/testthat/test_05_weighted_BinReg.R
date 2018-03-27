context("optional 'weights' argument")

# ------------------------------------------------------------------------------
# Fitting weighted densities for censored data
# ------------------------------------------------------------------------------
library(data.table)
library(dplyr)
library(tools)
library(mixtools)
set.seed(48925)

################################################################################
# simulate simple condensier-weights data using recommended method
################################################################################
sim_toy_weights_data <- function() {
  # big sample
  n <- 30000

  W1 <- rbinom(n, 1, 1/2)
  W2 <- rbinom(n, 1, 1/2)
  A <- rnorm(n, mean = 4 * W1 * W2, sd = 1)
  Pi <- plogis(-2 + 5 * W2)
  subsample <- rbinom(n, 1, Pi)

  # induce censoring in observed data
  DeltaA <- A
  DeltaA[subsample == 0] <- NA

  # set up observed data structure
  data_in <- as.data.table(list(W1 = W1, W2 = W2, DeltaA = DeltaA,
                                IPCW = 1 / Pi)) %>%
    as.data.table()

  ## **** UNCOMMENT TO SAVE THIS DATASET ****
  ## WILL SAVE TO CURRENT WORKING DIRECTORY
  notrun.save.toy.data.weights.1 <- function(data_in) {
    save(data_in, compress = TRUE, file = "./data_in_weights.rda",
         compression_level = 9)
    resaveRdaFiles("./data_in_weights.rda", compress = "bzip2")
  }
  return(data_in)
}

# generate and clean up data
data_in <- sim_toy_weights_data() %>%
  dplyr::filter(!is.na(DeltaA)) %>%
  as.data.table()

# perform condensier regression that we care about
fit_dens_out <- fit_density(X = "W1", Y = "DeltaA", weights = "IPCW",
                            input_data = data_in,
                            nbins = 50,
                            bin_method = "equal.mass",
                            bin_estimator = condensier::speedglmR6$new(),
                            parfit = FALSE)
fit_preds <- sample_value(model_fit = fit_dens_out, newdata = data_in)

# sample values should be N(0, 1) for W1 = 0:
test_that("Sample values are N(0, 1) for W1 = 0.", {
  meanA_W1is0 <- abs(mean(fit_preds[data_in$W1 == 0]))
  sdA_W1is0 <- sd(fit_preds[data_in$W1 == 0])

  expect_equal(meanA_W1is0, expected = 0, tol = 0.005)
  expect_equal(sdA_W1is0, expected = 1, tol = 0.02)
})

# sample values should be N(0, 1) for W1 = 0:
test_that("Sample values are a mixture of N(0, 1) and N(4, 1) for W1 = 1.", {
  # find means and SDs of Gaussian mixture using EM
  preds_fit_em <- normalmixEM(fit_preds[data_in$W1 == 1])

  meanA_W1is1 <- unlist(preds_fit_em["mu"], use.names = FALSE)
  sdA_W1is1 <- unlist(preds_fit_em["sigma"], use.names = FALSE)

  expect_equal(meanA_W1is1, expected = c(0, 4), tol = 0.02)
  expect_equal(sdA_W1is1, expected = c(1, 1), tol = 0.08)
})

################################################################################

################################################################################
# SETTING 1 FROM https://github.com/osofr/condensier/pull/12
# @BENKESER COMMENT 14 DECEMBER 2018
################################################################################

# $$W_1 \sim \text{Bin}(p = 0.5)$$
# $$\Delta \mid W_1 = w_1 \sim Bin(Logistic(loc = w_1))$$
# $$A \mid W_1 = w_1 \sim N(\mu = 2 \cdot w_1, \sigma = 1)$$
# $$O = (W_1, \Delta \cdot A)$$
sim_data_set1 <- function(n_obs = 1000, w_prob = 0.5) {
  w <- rbinom(n = n_obs, size = 1, prob = w_prob)
  delta <- rbinom(n = n_obs, size = 1, prob = plogis(w))
  a <- rnorm(n = n_obs, mean = 2 * w, sd = 1)
  data_in <- as.data.table(cbind(a, delta, w, 1 / plogis(w))) %>%
    setnames(., c("A", "Delta", "W", "Weights")) %>%
    dplyr::filter(Delta == 1) %>%
    dplyr::select(-Delta)
  return(data_in)
}

fit1 <- fit_density(X = "W", Y = "A", input_data = data_in,
                    bin_method = "equal.mass", nbins = n_bins,
                    bin_estimator = speedglmR6$new())
fit2 <- fit_density(X = "W", Y = "A", input_data = data_in,
                    bin_method = "equal.mass", nbins = n_bins,
                    bin_estimator = speedglmR6$new(), weights = "Weights")
preds1 <- sample_value(model_fit = fit1, newdata = data_in)
preds2 <- sample_value(model_fit = fit2, newdata = data_in)

################################################################################

################################################################################
# SETTING 2 FROM https://github.com/osofr/condensier/pull/12
# @BENKESER COMMENT 14 DECEMBER 2018
################################################################################

#$$W_1 \sim Binom(p = 0.5)$$
#$$W_2 \sim Binom(p = 0.5)$$
#$$\Delta \mid W_1 = w_1, W_2 = w_2 \sim Bin(Logistic(w_1 + w_2))$$
#$$A | W_1 = w_1, W_2 = w_2 ~ N(\mu = 2 \cdot w_1 \cdot w_2, \sigma = 1)$$
#$$O = (W_1, \Delta \cdot W_2, \Delta \cdot A)$$

sim_data_set2 <- function(n_obs = 1000, w_prob = 0.5) {
  w1 <- rbinom(n = n_obs, size = 1, prob = w_prob)
  w2 <- rbinom(n = n_obs, size = 1, prob = w_prob)
  delta <- rbinom(n = n_obs, size = 1, prob = plogis(w1 + w2))
  a <- rnorm(n = n_obs, mean = 2 * w1 * w2, sd = 1)
  data_in <- as.data.table(cbind(a, delta, w1, w2, 1 / plogis(w1 + w2))) %>%
    setnames(., c("A", "Delta", "W1", "W2", "Weights")) %>%
    dplyr::filter(Delta == 1) %>%
    dplyr::select(-Delta)
  return(data_in)
}

data_in <- sim_data_set2(n_obs = n_samp)
fit1 <- fit_density(X = c(paste0("W", seq_len(2))), Y = "A",
                    input_data = data_in,
                    bin_method = "equal.mass", nbins = n_bins,
                    bin_estimator = speedglmR6$new())
fit2 <- fit_density(X = c(paste0("W", seq_len(2))), Y = "A",
                    input_data = data_in,
                    bin_method = "equal.mass", nbins = n_bins,
                    bin_estimator = speedglmR6$new(), weights = "Weights")
preds1 <- sample_value(model_fit = fit1, newdata = data_in)
preds2 <- sample_value(model_fit = fit2, newdata = data_in)

