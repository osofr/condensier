context("optional 'weights' argument")

# ------------------------------------------------------------------------------
# Fitting weighted densities for censored data
# ------------------------------------------------------------------------------
library(data.table)
library(here)
library(dplyr)
library(tools)
library(mixtools)
set.seed(48925)

# save data for future tests
# uncomment and run to generate and save _new_ test data
save_toy_data_weights <- function(test_data_weights) {
  save(test_data_weights, compress = TRUE,
       file = here("data", "test_data_weights.rda"),
       compression_level = 9)
  resaveRdaFiles(here("data", "test_data_weights.rda"), compress = "bzip2")
}

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
  return(data_in)
}

# generate and clean up data
if (!file.exists(here("data", "test_data_weights.rda"))) {
  toy_data_with_weights <- sim_toy_weights_data()
} else {
  load(here("data", "test_data_weights.rda"))
}

# uncomment and run to save _new_ test data
#save_toy_data_weights(toy_data_with_weights)

# remove missingness
data_in <- test_data_weights %>%
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

preds1 <- predict_probability(model_fit = fit_dens_out, newdata = data_in)


fit_dens_out2 <- fit_density(X = "W1", Y = "DeltaA", weights = data_in[["IPCW"]],
                            input_data = data_in,
                            nbins = 50,
                            bin_method = "equal.mass",
                            bin_estimator = condensier::speedglmR6$new(),
                            parfit = FALSE)
preds2 <- predict_probability(model_fit = fit_dens_out2, newdata = data_in)

test_that(paste("Weights argument as name and as vector are identical",
                "Sample values are N(0, 1) for W1 = 0."), {
  expect_equal(sum(abs(preds1-preds2)), expected = 0, tol = 0.02)
})


# sample values should be N(0, 1) for W1 = 0:
test_that(paste("Weights argument works properly:",
                "Sample values are N(0, 1) for W1 = 0."), {
  meanA_W1is0 <- abs(mean(fit_preds[data_in$W1 == 0]))
  sdA_W1is0 <- sd(fit_preds[data_in$W1 == 0])

  expect_equal(meanA_W1is0, expected = 0, tol = 0.02)
  expect_equal(sdA_W1is0, expected = 1, tol = 0.085)
})

# sample values should be N(0, 1) for W1 = 0:
test_that(paste("Weights argument work properly:",
                "Sample values are mixture of N(0, 1) and N(4, 1) | W1 = 1."), {
  # find means and SDs of Gaussian mixture using EM
  suppressMessages(
    preds_fit_em <- mixtools::normalmixEM(fit_preds[data_in$W1 == 1])
  )

  meanA_W1is1 <- unlist(preds_fit_em["mu"], use.names = FALSE)
  sdA_W1is1 <- unlist(preds_fit_em["sigma"], use.names = FALSE)

  expect_equal(meanA_W1is1, expected = c(0, 4), tol = 0.012)
  expect_equal(sdA_W1is1, expected = c(1, 1), tol = 0.06)
})

