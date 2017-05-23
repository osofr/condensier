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

test_that("direct RegressionClass method should still work, means of the resampled values should be close to the observed means", {
  data_store_obj <- DataStore$new(input_data = data.table::data.table(datO), Y = "sA", X = c("W1", "W2", "W3"))
  # Define est_params_list:
  reg.sVars <- list(outvars = c("sA"), predvars = c("W1", "W2", "W3"))
  subset_vars <- lapply(reg.sVars$outvars, function(var) {var})
  sA_class <- data_store_obj$type.sVar[reg.sVars$outvars]
  # Put all est_params in RegressionClass (regression with speedglm package)
  regclass <- RegressionClass$new(outvar.class = sA_class,
                                  outvar = reg.sVars$outvars,
                                  predvars = reg.sVars$predvars,
                                  subset = subset_vars,
                                  nbins=50)
  summeas.g0 <- SummariesModel$new(reg = regclass, data_object = data_store_obj)
  summeas.g0$fit(data = data_store_obj)
  summeas.g0$predict(newdata = data_store_obj)  # summeas.g0$sA_nms

  ## resample the outcome (sA), conditional on observed predictors, using the current conditional density fit for P(sA|X)
  set.seed(123456)
  resampA_20times <- lapply(1:20, function(i) summeas.g0$sampleA(newdata = data_store_obj))
  resampA_combined <- unlist(resampA_20times)
  resampA_average <- rowMeans(do.call("cbind", resampA_20times))
  # summary(resampA_average - datO[["sA"]])
  # print(mean(resampA_average) - mean(datO[["sA"]]))
  expect_true(mean(resampA_average) - mean(datO[["sA"]]) < 0.005)

  # Get P(sA|W) for the observed data (W,sA):
  h_gN <- summeas.g0$predictAeqa(newdata = data_store_obj) # *** DataStore$O.datnetA IS TO BE RENAMED TO $O.O.datnetA for clarity ***

  ## ---------------------------------------------------------------------------------------------------------
  ## Plot density of the entire observed outcome vs. resampled outcomes (resampled 20 times and then combined into one long vector)
  ## ---------------------------------------------------------------------------------------------------------
  set.seed(12345)
  dens_A <- density(datO[["sA"]])
  dens_resampA <- density(resampA_combined)

  ## plot the density fit with observed A vs. resampled A:
  plot(dens_A)
  lines(dens_resampA, col = "blue")

  ## numerically compare the two density fits:
  common_x <- intersect(round(dens_A[["x"]], 2), round(dens_resampA[["x"]], 2))
  idx_A <- which(round(dens_A[["x"]], 2) %in% common_x)
  idx_resampA <- which(round(dens_resampA[["x"]], 2) %in% common_x)
  length(idx_resampA) == length(idx_A)
  ## MSE on the density fit of the resampled A values (compared to density fit of the observed A values)
  MSE_resampA <- mean((dens_A[["y"]][idx_A] - dens_resampA[["y"]][idx_resampA])^2)
  # print("MSE for observed density fit vs. resampled: " %+% MSE_resampA)
  expect_true(MSE_resampA < 1e-5)

  ## ---------------------------------------------------------------------------------------------------------
  ## Same as above but conditional on some FIXED values of predictors (W1,W2,W3)
  ## + adding the predicted discretized probs conditional on same W's
  ## ---------------------------------------------------------------------------------------------------------
  setWvals <- c(W1 = 0, W2 = 0, W3 = 1)
  subs <- (datO$W1==setWvals["W1"] & datO$W2==setWvals["W2"] & datO$W3==setWvals["W3"])

  ## helper fun
  get.setW.sAdat <- function(setWvals, nsamp = 100000) {
    require(simcausal)
    W1val <- setWvals[1]; W2val <- setWvals[2]; W3val <- setWvals[3]
    D <- DAG.empty()
    D <-
    D + node("W1", distr = "rconst", const = .(W1val)) +
        node("W2", distr = "rconst", const = .(W2val)) +
        node("W3", distr = "rconst", const = .(W3val)) +
        node("sA", distr = "rnorm", mean = (0.98 * W1 + 0.58 * W2 + 0.33 * W3), sd = 1)
    D <- set.DAG(D, n.test = 10)
    datWset <- sim(D, n = nsamp, rndseed = 12345)
    setWmat <- as.matrix(data.frame(W1 = W1val, W2 = W2val, W3 = W3val, sA = seq(-4, 4, by = 0.2)))
    return(list(setWsA = datWset$sA, setWmat = setWmat))
  }

  suppressWarnings( setWdat_res <- get.setW.sAdat(setWvals, nsamp))

  ## plot densitity first:
  plot(density(setWdat_res$setWsA))
  lines(datO[subs,][["sA"]], h_gN[subs], type = "p", cex = .3, col = "red")
  set.seed(12345)
  resampA <- summeas.g0$sampleA(newdata = data_store_obj)
  ## add the density of the resamp outcomes:
  lines(density(resampA[subs]), col = "blue")
})
