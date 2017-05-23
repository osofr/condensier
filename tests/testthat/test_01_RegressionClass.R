# helper function for generating some data
get.testDat <- function(nsamp = 100000) {
  `%+%` <- function(a, b) paste0(a, b)
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

test_that("testing RegressionClass and SummariesModel factorize the joint likelihood as expected", {
  # Tests for RegressionClass:
  reg_test1 <- RegressionClass$new(outvar.class = c(condensier:::gvars$sVartypes$bin, condensier:::gvars$sVartypes$bin),
                                  outvar = c("A1", "A2"),
                                  predvars = c("W1", "W2"))
  class(reg_test1)
  reg_test1$subset
  model1 <- SummariesModel$new(reg = reg_test1)
  # [1] "Init BinOutModel:"
  # [1] "P(A1|W1,W2)"
  # [1] "Init BinOutModel:"
  # [1] "P(A2|A1,W1,W2)"
  expect_true("BinOutModel" %in% class(model1$getPsAsW.models()$`P(sA|sW).1`))

  reg_test2 <- RegressionClass$new(outvar.class = c(condensier:::gvars$sVartypes$bin, condensier:::gvars$sVartypes$bin),
                                  outvar = c("A1", "A2"),
                                  predvars = c("W1", "W2"),
                                  subset = list("A1"))
  class(reg_test2)
  reg_test2$subset
  model2 <- SummariesModel$new(reg = reg_test2)

  expect_error(
    reg_test3 <- RegressionClass$new(outvar.class = c(condensier:::gvars$sVartypes$cont, condensier:::gvars$sVartypes$cont),
                                    outvar = c("sA"), predvars = c("W1", "W2", "W3"),
                                    subset = list(quote(sA==0)))
    )

  reg_test3 <- RegressionClass$new(outvar.class = condensier:::gvars$sVartypes$cont,
                                  outvar = c("sA"), predvars = c("W1", "W2", "W3"),
                                  subset = list("sA"))
                                  # subset = list(quote(sA==0)))

  reg_test_new <- reg_test3$clone()
  expect_true(all.equal(reg_test3, reg_test_new))
  expect_true("RegressionClass" %in% class(reg_test3))
})


test_that("testing that bin intervals are detected with minimal RegressionClass inputs", {
  reg_test3 <- RegressionClass$new(outvar.class = condensier:::gvars$sVartypes$cont,
                                  outvar = c("sA"), predvars = c("W1", "W2", "W3"),
                                  subset = list("sA"))
  nsamp <- 1000
  datO <- get.testDat(nsamp)
  data_obj <- DataStore$new(input_data = datO, Y = "sA", X = c("W1", "W2", "W3"))
  # :
  intervls <- data_obj$detect.sVar.intrvls(reg_test3$outvar,
                                           nbins = reg_test3$nbins,
                                           bin_bymass = reg_test3$bin_bymass,
                                           bin_bydhist = reg_test3$bin_bydhist,
                                           max_nperbin = reg_test3$max_nperbin)
  class(reg_test3$subset[[1]])
  model3 <- SummariesModel$new(reg = reg_test3, data_object = data_obj)
})
