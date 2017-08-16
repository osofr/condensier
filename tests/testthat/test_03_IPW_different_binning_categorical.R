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
      node("sA", distr = "rnorm", mean = sA.mu, sd = 1) +
      node("sA.cat", distr = "rconst", const = as.integer(sA))
  D <- set.DAG(D, n.test = 10)
  datO <- sim(D, n = nsamp, rndseed = 12345)
}

nsamp <- 10000
datO <- data.table::data.table(get.density.sAdat(nsamp))
ncats <- sort(unique(datO[["sA.cat"]]))
print(ncats)
newdata <- datO[1:100, c("W1", "W2", "W3", "sA.cat"), with = FALSE]

test_that("categorical input variable should be typed as such", {
  data_store_obj <- DataStore$new(input_data = data.table::data.table(datO), Y = "sA.cat", X = c("W1", "W2", "W3"))
  reg.sVars <- list(outvars = c("sA.cat"), predvars = c("W1", "W2", "W3"))
  subset_vars <- lapply(reg.sVars$outvars, function(var) {var})
  sA_class <- data_store_obj$type.sVar[reg.sVars$outvars]
  expect_true(sA_class[["sA.cat"]] %in% "categor")
})

test_that("CategorSummaryModel is engaged appropriately and the cut-offs for the categorical input var are defined correctly", {
  data_store_obj <- DataStore$new(input_data = data.table::data.table(datO), Y = "sA.cat", X = c("W1", "W2", "W3"))
  # Put all est_params in RegressionClass (regression with speedglm package)
  # print("fitting h_gN based equal.len intervals (default) and speedglm (default): ")
  regclass <- RegressionClass$new(outvar.class = list(sA.cat = "categor"),
                                  outvar = "sA.cat",
                                  predvars = c("W1", "W2", "W3"),
                                  subset = list("sA.cat"),
                                  bin_bymass = FALSE,
                                  bin_bydhist = FALSE,
                                  max_nperbin = 1000
                                  # nbins = 10
                                  )
  summeas.g0 <- SummariesModel$new(reg = regclass, data_object = data_store_obj)
  summeas.g0$fit(data = data_store_obj)
  out_CatSummaryModel <- summeas.g0$getPsAsW.models()$`P(sA|sW).1`

  expect_true(is(out_CatSummaryModel, "CategorSummaryModel"))
  expect_true(all.equal(ncats, out_CatSummaryModel$levels))
  expect_true(length(ncats) == out_CatSummaryModel$nbins)

  ## get individual fits for each bin hazard:
  out_BinModels <- out_CatSummaryModel$getPsAsW.models()
  print(condensier:::coef.BinOutModel(out_BinModels[[1]]))
  print(condensier:::summary.BinOutModel(out_BinModels[[2]]))

  summeas.g0$predict(newdata = data_store_obj)  # summeas.g0$sA_nms
  # Get P(sA|W) for the observed data (W,sA):
  # SHOULD BE SIMILAR TO THE OBSERVED MASS PROB FUNCTION OF sA.cat
  h_gN <- summeas.g0$predictAeqa(newdata = data_store_obj)
  expect_true(abs(mean(h_gN)-0.4064357) < 10^-4)

})

# # ---------------------------------------------------------------------------------------------------------
# # Plot predicted discretized probs conditional on some values of W's
# # ---------------------------------------------------------------------------------------------------------
# setWvals <- c(W1 = 0, W2 = 0, W3 = 1)
# subs <- (datO$W1==setWvals["W1"] & datO$W2==setWvals["W2"] & datO$W3==setWvals["W3"])
# # sum(subs)

# ## helper fun
# get.setW.sAdat <- function(setWvals, nsamp = 100000) {
#   require(simcausal)
#   W1val <- setWvals[1]; W2val <- setWvals[2]; W3val <- setWvals[3]
#   D <- DAG.empty()
#   D <-
#   D + node("W1", distr = "rconst", const = .(W1val)) +
#       node("W2", distr = "rconst", const = .(W2val)) +
#       node("W3", distr = "rconst", const = .(W3val)) +
#       node("sA", distr = "rnorm", mean = (0.98 * W1 + 0.58 * W2 + 0.33 * W3), sd = 1) +
#       node("sA.cat", distr = "rconst", const = as.integer(sA))
#   D <- set.DAG(D, n.test = 10)
#   datWset <- sim(D, n = nsamp, rndseed = 12345)
#   setWmat <- as.matrix(data.frame(W1 = W1val, W2 = W2val, W3 = W3val, sA.cat = seq(-3, 4, by = 0.2)))
#   return(list(setWsA = datWset$sA.cat, setWmat = setWmat))
# }

# suppressWarnings( setWdat_res <- get.setW.sAdat(setWvals, nsamp) )

# ##plot densitity first:
# hist(setWdat_res$setWsA, probability = TRUE)
# lines(datO[subs,][["sA.cat"]], h_gN[subs], type = "p", cex = .3, col = "red")
# resampA <- summeas.g0$sampleA(newdata = data_store_obj)
# lines(density(resampA[subs]), col = "blue")

# ##plot predicted vals first:
# plot(datO[subs][["sA.cat"]], h_gN[subs], type = "p", cex = .3, col = "red")
# lines(hist(setWdat_res$setWsA, probability = TRUE))
