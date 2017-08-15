context("different binning methods")

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

nsamp <- 10000
datO <- data.table::data.table(get.density.sAdat(nsamp))
newdata <- datO[1:5, c("W1", "W2", "W3", "sA"), with = FALSE]

test_that("different binning methods should give the right mean predicted probability", {

  data_store_obj <- DataStore$new(input_data = data.table::data.table(datO), Y = "sA", X = c("W1", "W2", "W3"))
  reg.sVars <- list(outvars = c("sA"), predvars = c("W1", "W2", "W3"))
  subset_vars <- lapply(reg.sVars$outvars, function(var) {var})
  sA_class <- data_store_obj$type.sVar[reg.sVars$outvars]

  # Put all est_params in RegressionClass (regression with speedglm package)
  # print("fitting h_gN based equal.len intervals (default) and speedglm (default): ")
  regclass <- RegressionClass$new(outvar.class = sA_class,
                                  outvar = reg.sVars$outvars,
                                  predvars = reg.sVars$predvars,
                                  subset = subset_vars,
                                  bin_bymass = FALSE,
                                  bin_bydhist = FALSE,
                                  max_nperbin = 1000
                                  # nbins = 10
                                  )

  summeas.g0 <- SummariesModel$new(reg = regclass, data_object = data_store_obj)
  summeas.g0$fit(data = data_store_obj)

  # Test the coef and summary functions for binoutmodel class:
  out_ContinSummaryModel <- summeas.g0$getPsAsW.models()$`P(sA|sW).1`
  out_BinModels <- out_ContinSummaryModel$getPsAsW.models()
  # print(condensier:::coef.BinOutModel(out_BinModels[[1]]))
  # print(condensier:::summary.BinOutModel(out_BinModels[[2]]))

  # (intrvls <- out_ContinSummaryModel$intrvls)
  # (intrvls.width <- diff(intrvls))
  # length(intrvls.width)
  # ord.sVar <- nodeobjs$data_store_obj$discretize.sVar(name.sVar = "sA", intervals = out_ContinSummaryModel$intrvls)
  # ord.sVar_bw <- intrvls.width[ord.sVar]
  # print(head(cbind(sA = nodeobjs$data_store_obj$dat.sVar[, "sA"], ord.sVar, bw = ord.sVar_bw, nodeobjs$data_store_obj$dat.bin.sVar), 5))
  # print("freq count for transformed ord.sVar: "); print(table(ord.sVar))
  # plot(density(ord.sVar))
  # hist(ord.sVar)
  summeas.g0$predict(newdata = data_store_obj)  # summeas.g0$sA_nms

  # Get P(sA|W) for the observed data (W,sA):
  # SHOULD BE SIMILAR TO THE OBSERVED DENSITY OF s.A (but discretized)
  h_gN <- summeas.g0$predictAeqa(newdata = data_store_obj)

  # print("h_gN fit under speedglm: " %+% mean(h_gN)) # [1] 0.2718823
  expect_true(abs(mean(h_gN)-0.2718823) < 10^-4)
  # ---------------------------------------------------------------------------------------------------------
  # Plot predicted discretized probs conditional on some values of W's
  # ---------------------------------------------------------------------------------------------------------
  setWvals <- c(W1 = 0, W2 = 0, W3 = 1)
  subs <- (datO$W1==setWvals["W1"] & datO$W2==setWvals["W2"] & datO$W3==setWvals["W3"])
  # sum(subs)

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

  suppressWarnings( setWdat_res <- get.setW.sAdat(setWvals, nsamp) )

  ##plot densitity first:
  plot(density(setWdat_res$setWsA))
  lines(datO[subs,][["sA"]], h_gN[subs], type = "p", cex = .3, col = "red")

  ##plot predicted vals first:
  plot(datO[subs][["sA"]], h_gN[subs], type = "p", cex = .3, col = "red")
  lines(density(setWdat_res$setWsA))

  resampA <- summeas.g0$sampleA(newdata = data_store_obj)
  # lines(density(resampA[subs]), col = "blue")

  # ---------------------------------------------------------------------------------------------------------
  # Plot all predicted discretized probs together (without conditioning on particular subset of W's)
  # ---------------------------------------------------------------------------------------------------------
  # plot(datO[,"sA"], h_gN, type = "p", cex = .3, col = "red")
  # lines(density(datO[,"sA"]))

  # ---------------------------------------------------------------------------------------------------------
  # **** Same fit as before but doing regressions with stats::glm.fit ****
  # ---------------------------------------------------------------------------------------------------------
  # print("fitting h_gN based equal.len intervals (default) and glm.fit: ")
  regclass.gml <- RegressionClass$new(bin_estimator = glmR6$new(), outvar.class = sA_class,
                                      outvar = reg.sVars$outvars,
                                      predvars = reg.sVars$predvars,
                                      subset = subset_vars,
                                      bin_bymass = FALSE,
                                      bin_bydhist = FALSE,
                                      max_nperbin = 1000)
  summeas.g0.glm <- SummariesModel$new(reg = regclass.gml, data_object = data_store_obj)
  summeas.g0.glm$fit(data = data_store_obj)

  # Test the coef and summary functions for binoutmodel class:
  out_ContinSummaryModel <- summeas.g0.glm$getPsAsW.models()$`P(sA|sW).1`
  out_BinModels <- out_ContinSummaryModel$getPsAsW.models()
  print(condensier:::coef.BinOutModel(out_BinModels[[1]]))
  print(str(condensier:::summary.BinOutModel(out_BinModels[[2]])))

  summeas.g0.glm$predict(newdata = data_store_obj)
  h_gN.glm <- summeas.g0.glm$predictAeqa(newdata = data_store_obj) # *** DataStore$O.datnetA IS TO BE RENAMED TO $O.O.datnetA for clarity ***

  # print("h_gN fit under glm: " %+% mean(h_gN.glm)) # [1] 0.2718823
  expect_true(abs(mean(h_gN.glm)-0.2718823) < 10^-4)
  # plot densitity first:
  # plot(density(setWdat_res$setWsA))
  # lines(datO[subs][["sA"]], h_gN[subs], type = "p", cex = .3, col = "red")

  # ---------------------------------------------------------------------------------------------------------
  # **** Same fit as before but doing binnning by mass intervals & regressions with speedglm ****
  # ---------------------------------------------------------------------------------------------------------
  # options(tmlenet.verbose = TRUE)
  # gvars$verbose <- TRUE  # set to TRUE after loading all package functions to print all output
  # ls(gvars)

  # print("fitting h_gN based on bin_bymass = TRUE and speedglm (default): ")
  regclass.binmass <- RegressionClass$new(bin_estimator = speedglmR6$new(),
                                          bin_bymass = TRUE,
                                          max_nperbin = 1000,
                                          outvar.class = sA_class,
                                          outvar = reg.sVars$outvars,
                                          predvars = reg.sVars$predvars,
                                          subset = subset_vars)
  summeas.g0 <- SummariesModel$new(reg = regclass.binmass, data_object = data_store_obj)
  summeas.g0$fit(data = data_store_obj)
  summeas.g0$predict(newdata = data_store_obj)
  h_gN <- summeas.g0$predictAeqa(newdata = data_store_obj) # *** DataStore$O.datnetA IS TO BE RENAMED TO $O.O.datnetA for clarity ***
  # mean(h_gN) # [1] 0.2668144
  expect_true(abs(mean(h_gN)-0.2668144) < 10^-4)

  # get get the observed data:
  # data_store_obj$dat.sVar

  # plot densitity first:
  # plot(density(setWdat_res$setWsA))
  # lines(datO[subs][["sA"]], h_gN[subs], type = "p", cex = .3, col = "red")

  # ---------------------------------------------------------------------------------------------------------
  # **** Same fit as before but binning using "dhist" function & regressions with speedglm ****
  # ---------------------------------------------------------------------------------------------------------
  # print("fitting h_gN based on dhist intervals and speedglm (default): ")
  regclass.bidhist <- RegressionClass$new(bin_estimator = speedglmR6$new(),
                                          bin_bymass = FALSE,
                                          bin_bydhist = TRUE,
                                          max_nperbin = 1000,
                                          outvar.class = sA_class,
                                          outvar = reg.sVars$outvars,
                                          predvars = reg.sVars$predvars,
                                          subset = subset_vars)
  summeas.g0 <- SummariesModel$new(reg = regclass.bidhist, data_object = data_store_obj)
  summeas.g0$fit(data = data_store_obj)

  out_ContinSummaryModel <- summeas.g0$getPsAsW.models()$`P(sA|sW).1`
  # print("Intervals for dhist: ")
  # print(out_ContinSummaryModel$intrvls)
  # [1] -1003.5240566    -3.5240566    -1.9562682    -0.8766233    -0.2052710     0.2963970     0.7404754     1.2078705
  # [9]     1.6959939     2.3455552     3.3931981     4.9362651  1004.9362651

  out_BinModels <- out_ContinSummaryModel$getPsAsW.models()
  # print(condensier:::coef.BinOutModel(out_BinModels[[1]]))
  # print(condensier:::summary.BinOutModel(out_BinModels[[2]]))

  summeas.g0$predict(newdata = data_store_obj)

  h_gN <- summeas.g0$predictAeqa(newdata = data_store_obj) # *** DataStore$O.datnetA IS TO BE RENAMED TO $O.O.datnetA for clarity ***
  # print("regclass.bidhist mean(h_gN) under speedglm: " %+% mean(h_gN))
  # [1] 0.276875
  # [1] 0.27687500785635
  expect_true(abs(mean(h_gN)-0.276875) < 10^-4)

})