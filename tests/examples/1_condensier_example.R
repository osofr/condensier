library("simcausal")
D <- DAG.empty()
D <-
D + node("W1", distr = "rbern", prob = 0.5) +
  node("W2", distr = "rbern", prob = 0.3) +
  node("W3", distr = "rbern", prob = 0.3) +
  node("sA.mu", distr = "rconst", const = (0.98 * W1 + 0.58 * W2 + 0.33 * W3)) +
  node("sA", distr = "rnorm", mean = sA.mu, sd = 1)
D <- set.DAG(D, n.test = 10)
datO <- sim(D, n = 10000, rndseed = 12345)

## Fit conditional density using equal mass bins (same number of obs per bin):
dens_fit <- fit_density(
    X = c("W1", "W2", "W3"),
    Y = "sA",
    input_data = datO,
    nbins = 20,
    bin_method = "equal.mass",
    bin_estimator = speedglmR6$new())

## Wrapper function to predict the conditional probability (likelihood)
## for new observations:
newdata <- datO[1:5, c("W1", "W2", "W3", "sA"), with = FALSE]
preds <- predict_probability(dens_fit, newdata)

## Wrapper function to sample the values from the conditional density fit:
sampledY <- sample_value(dens_fit, newdata)

## Fit conditional density using equal length bins
## (divides the range of support of the outcome to define each bin):
dens_fit <- fit_density(
    X = c("W1", "W2", "W3"),
    Y = "sA",
    input_data = datO,
    nbins = 20,
    bin_method = "equal.len",
    bin_estimator = speedglmR6$new())

## Wrapper function to predict the conditional probability (likelihood)
## for new observations:
newdata <- datO[1:5, c("W1", "W2", "W3", "sA"), with = FALSE]
preds <- predict_probability(dens_fit, newdata)

## Wrapper function to sample the values from the conditional density fit:
sampledY <- sample_value(dens_fit, newdata)


## Fit conditional density using custom bin definitions (argument intrvls):
dens_fit <- fit_density(
    X = c("W1", "W2", "W3"),
    Y = "sA",
    input_data = datO,
    bin_estimator = speedglmR6$new(),
    intrvls = list(sA = seq(-4,4, by = 0.1)))


## Fit conditional density using custom bin definitions and
## pool all bin indicators into a single long-format dataset.
## The pooling results in a single regression that is fit for all bin hazards,
## with a bin indicator added as an additional covariate.
dens_fit <- fit_density(
    X = c("W1", "W2", "W3"),
    Y = "sA",
    input_data = datO,
    bin_estimator = speedglmR6$new(),
    intrvls = list(sA = seq(-4,4, by = 0.1)),
    pool = TRUE)
