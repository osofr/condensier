context("optional 'weights' argument")

# ------------------------------------------------------------------------------
# Fitting weighted densities for censored data
# ------------------------------------------------------------------------------
library(data.table)
library(dplyr)
set.seed(48925)

# big sample 
n <- 30000

W1 <- rbinom(n, 1, 1/2)
W2 <- rbinom(n, 1, 1/2)
A <- rnorm(n, mean = 4 * W1 * W2, sd = 1)
Pi <- plogis(-2 + 5 * W2)
subsample <- rbinom(n, 1, Pi)

DeltaA <- A
DeltaA[subsample == 0] <- NA

# now say we want density of A | W1
# let's plot full data vs. observed data parameters
#layout(matrix(1:4, 2, 2, byrow = TRUE))
x_seq <- seq(-4, 8, length = 5000)

# full data parameter is N(0,1) | W1 = 0
hist(A[W1 == 0], freq = FALSE)
lines(y = dnorm(x_seq), x = x_seq, col = 2)

# observed data parameter is the same
hist(DeltaA[W1 == 0], freq = FALSE)
lines(y = dnorm(x_seq), x = x_seq, col = 2)

# full data parameter should be an equal mixture of N(0,1) and N(4,1) | W1 = 1
hist(A[W1 == 1], freq = FALSE)
lines(y = 1/2 * (dnorm(x_seq) + dnorm(x_seq, 4, 1)), x = x_seq, col = 2)

# but observed data parameter is not, because of biased sampling
hist(DeltaA[W1 == 1], freq = FALSE)
lines(y = 1/2 * (dnorm(x_seq) + dnorm(x_seq, 4, 1)), x = x_seq, col = 2)

# check whether condensier recovers the proper mixture | W1 = 1
# when using true IPCW weights and regressing DeltaA ~ W1 | DeltaA != NA
data_in <- as.data.table(list(W1 = W1, W2 = W2, A = DeltaA, IPCW = 1 / Pi)) %>%
  dplyr::filter(!is.na(A)) %>%
  as.data.table()

# THIS APPEARS TO BE THE REGRESSION THAT WE CARE ABOUT
fit_dens_out <- fit_density(X = "W1", Y = "A", weights = "IPCW",
                            input_data = data_in,
                            nbins = 50,
                            bin_method = "equal.mass",
                            bin_estimator = condensier::speedglmR6$new(),
                            parfit = FALSE)
fit_preds <- sample_value(model_fit = fit_dens_out, newdata = data_in)
# BUT IT PRODUCES UNEXPECTED RESULTS
hist(fit_preds[W1 == 0], freq = FALSE)
hist(fit_preds[W1 == 1], freq = FALSE)
hist(fit_preds, freq = FALSE)


# THIS IS NOT THE SPECIFIED REGRESSION OF INTEREST
fit_dens_out <- fit_density(X = "W2", Y = "A", weights = "IPCW",
                            input_data = data_in,
                            nbins = 50,
                            bin_method = "equal.mass",
                            bin_estimator = condensier::speedglmR6$new(),
                            parfit = FALSE)
fit_preds <- sample_value(model_fit = fit_dens_out, newdata = data_in)
# BUT THIS PRODUCES THE RESULTS WE CARE ABOUT
hist(fit_preds[W1 == 0], freq = FALSE)
hist(fit_preds[W1 == 1], freq = FALSE)
hist(fit_preds, freq = FALSE)
