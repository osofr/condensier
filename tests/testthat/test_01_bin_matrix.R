library(mockery)
context("DataStore Bin Matrix")

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

test_that("automatic outcome type recognition", {
  test_mat <- as.matrix(data.frame(a = c(0,1,0,0,1), b = rep(5,5), c = c(1,2,3,4,5), d = rnorm(5)))
  correct.types <- list(a = condensier:::gvars$sVartypes$bin,
                        b = condensier:::gvars$sVartypes$bin,
                        c = condensier:::gvars$sVartypes$cat,
                        d = condensier:::gvars$sVartypes$cont)
  out.types <- condensier:::detect.col.types(test_mat, max_n_cat = 10)
  expect_true(all.equal(correct.types, out.types))
})

test_that("binary outcome (treated as continuous), binning results in same outcome as original binary.", {
  nbins <- 10L
  binvar <- rbinom(n=100, size = 1, prob = 0.3)

  # binning by equal length
  int_bylen <- c(0, 0.1, 0.9, 1)
  ord1 <- findInterval(x = binvar, vec = int_bylen, rightmost.closed = TRUE)
  bins_1 <- condensier:::make.bins_mtx_1(x.ordinal = ord1, nbins = nbins, bin.nms = paste0("B_",1:nbins))
  res1 <- cbind(binvar = binvar, ord1 = ord1, bins_1)

  # binning by equal mass
  int_bymass <- quantile(binvar, prob = c(0, 0.1, 0.9, 1))
  ord2 <- findInterval(x = binvar, vec = int_bymass, rightmost.closed = TRUE)
  bins_2 <- condensier:::make.bins_mtx_1(x.ordinal = ord1, nbins = nbins, bin.nms = paste0("B_",1:nbins))
  res2 <- cbind(binvar = binvar, ord2 = ord2, bins_2)
})


# making bin indicator matrix from ordinal sVar:
test_that("bin indicator matrix for ordinal variable that's in range 1:nbins", {
  ord_vec_1 <- sample(c(1:5), 200, replace=TRUE)
  nbins <- length(unique(ord_vec_1))
  # levels <- sort(unique(ord_vec_1))
  bin.nms <- paste0("B_",(1:nbins))
  bin_ord_vec_1 <- condensier:::make.bins_mtx_1(ord_vec_1, nbins = nbins, bin.nms = bin.nms)
  out_mat_1 <- cbind(ord_vec_1 = ord_vec_1, bin_ord_vec_1)
})


test_that("bin indicator matrix for ordinal variable that's not in standard range 1:nbins", {
  ord_vec_2 <- sample(c(-2,1:5,8,9), 200, replace=TRUE)
  nbins <- length(unique(ord_vec_2))
  levels <- sort(unique(ord_vec_2))
  bin.nms <- "B_"%+%levels
  bin_ord_vec_2 <- condensier:::make.bins_mtx_1(ord_vec_2, nbins = nbins, bin.nms = bin.nms, levels = levels)
  out_mat_2 <- cbind(ord_vec_2 = ord_vec_2, bin_ord_vec_2)
})


test_that("continuous and categorical interval cuttofs are created as expected", {
  nsamp <- 1000
  nbins <- 10
  # oldopts <- condensier_options(max_n_cat = 5, nbins = nbins)
  # ----------------------------------------------------------------------------------------
  # Continuous
  # ----------------------------------------------------------------------------------------
  datO <- get.testDat(nsamp = nsamp)
  data.table::setDT(datO)
  Kmax <- 3L
  datO[, "nF" := sample(1L:Kmax, nrow(datO), replace = TRUE)]
  datNetObs <- DataStore$new(input_data = datO, Y = "sA", X = c("W1", "W2", "W3", "nF"), max_n_cat = 5)

  # datNetObs <- makedat(nsamp=nsamp, Kmax=3)
  obsdat.sW <- datNetObs$dat.sVar
  print("head(obsdat.sW)"); print(head(obsdat.sW,10))

  print("Var types: "); str(datNetObs$type.sVar)
  print("Contin covars names: "); print(datNetObs$names.c.sVar)
  defints1 <- datNetObs$detect.sVar.intrvls("sA",
                                            nbins = nbins,
                                            bin_bymass = TRUE,
                                            bin_bydhist = FALSE,
                                            max_nperbin = 100)
  print("No normalization bin intervals by mass: ");
  print(defints1)

  # ----------------------------------------------------------------------------------------
  # Categorical w ncats < max_n_cat
  # ----------------------------------------------------------------------------------------
  defints2 <- datNetObs$detect.sVar.intrvls("nF",
                                            nbins = nbins,
                                            bin_bymass = TRUE,
                                            bin_bydhist = FALSE,
                                            max_nperbin = 100)

  # Testing ordinals with ncats < nbins get nbins = ncats:
  expect_true((length(defints2)-1) < nbins)
  # Testing all categories in ordinal are represented for categorical with 3 values (nF):
  expect_true(all(unique(datNetObs$dat.sWsA[["nF"]])%in%defints2))
})

