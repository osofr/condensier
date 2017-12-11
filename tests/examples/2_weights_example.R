## set up problem
library("condensier")
set.seed(638549)
n_obs <- 1000
n_w <- 1
a1_mean <- 2
a0_mean <- 0

## baseline covariate -- simple, binary
W <- as.numeric(replicate(n_w, rbinom(n_obs, 1, 0.5)))

## set and organize treatment based on baseline W
A1 <- rnorm(length(which(W == 1)), mean = a1_mean, sd = 1)
A0 <- rnorm(length(which(W == 0)), mean = a0_mean, sd = 1)
A <- rep(NA, n_obs)
A[which(W == 0)] <- A0
A[which(W == 1)] <- A1

## random censoring (1 is not censored)
C <- rbinom(n_obs, size = 1, prob = 0.75)

## create outcome
Y <- A + W + rnorm(n_obs)

## fit IPC weights for censoring
ipc_mod <- stats::glm(C ~ W, family = stats::binomial)
ipcw_probs <- stats::predict(ipc_mod, newdata = as.data.frame(cbind(W, C)))
ipc_weights <- C / as.numeric(ipcw_probs)
ipc_weights_out <- ipc_weights[ipc_weights != 0]

## fit conditional densities with IPC weights
data_O <- as.data.frame(cbind(W, A, C, Y)) %>%
    dplyr::filter(C == 1 )
data_O$weights <- ipc_weights_out
fit <- condensier::fit_density(X = "W", Y = "A", input_data = data_O,
                               weights = "weights")
pred <- condensier::predict_probability(model_fit = fit, newdata = data_O)

