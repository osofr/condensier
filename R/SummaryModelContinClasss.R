# -------------------------------------------------------------------------------------------
# same code in ContinSummaryModel$new and CategorSummaryModel$new replaced with outside function:
# Define subset evaluation for new bins:
# ******************************************************
# NOTE: Subsetting by var name only (which automatically evaluates as !gvars$misval(var)) for speed & memory efficiency
# ******************************************************
# -------------------------------------------------------------------------------------------
def_regs_subset <- function(self) {
  bin_regs <- self$reg$clone() # instead of defining new RegressionClass now cloning parent reg object and then ADDING new SETTINGS
  bin_regs$reg_hazard <- TRUE # don`t add degenerate bins as predictors in each binary regression
  if (!self$reg$pool) {
    add.oldsubset <- TRUE
    new.subsets <- lapply(self$reg$bin_nms,
                              function(var) {
                                res <- var
                                if (add.oldsubset) res <- c(res, self$reg$subset)
                                res
                              })

    new.sAclass <- as.list(rep_len(gvars$sVartypes$bin, self$reg$nbins))
    names(new.sAclass) <- self$reg$bin_nms
    bin_regs$ChangeOneToManyRegresssions(regs_list = list(outvar.class = new.sAclass,
                                                          outvar = self$reg$bin_nms,
                                                          predvars = self$reg$predvars,
                                                          subset = new.subsets))
    bin_regs$subset
  # Same but when pooling across bin indicators:
  } else {
    bin_regs$outvar.class <- gvars$sVartypes$bin
    bin_regs$outvar <- self$outvar
    bin_regs$outvars_to_pool <- self$reg$bin_nms
    if (gvars$verbose)  {
      print("pooled bin_regs$outvar: "); print(bin_regs$outvar)
      print("bin_regs$outvars_to_pool: "); print(bin_regs$outvars_to_pool)
      print("bin_regs$subset: "); print(bin_regs$subset)
    }
  }
  bin_regs$resetS3class()
  return(bin_regs)
}

# -------------------------------------------------------------------------------------------

## ---------------------------------------------------------------------
#' R6 class for fitting and predicting joint probability for a univariate continuous summary measure sA[j]
#'
#' This R6 class defines and fits a conditional probability model \code{P(sA[j]|sW,...)} for a univariate
#'  continuous summary measure \code{sA[j]}. This class inherits from \code{\link{SummariesModel}} class.
#'  Defines the fitting algorithm for a regression model \code{sA[j] ~ sW + ...}.
#'  Reconstructs the likelihood \code{P(sA[j]=sa[j]|sW,...)} afterwards.
#'  Continuous \code{sA[j]} is discretized using either of the 3 interval cutoff methods,
#'  defined via \code{\link{RegressionClass}} object \code{reg} passed to this class constructor.
#'  The fitting algorithm estimates the binary regressions for hazard \code{Bin_sA[j][i] ~ sW},
#'  i.e., the probability that continuous \code{sA[j]} falls into bin \code{i}, \code{Bin_sA[j]_i},
#'  given that \code{sA[j]} does not belong to any prior bins \code{Bin_sA[j]_1, ..., Bin_sA[j]_{i-1}}.
#'  The dataset of discretized summary measures (\code{BinsA[j]_1,...,BinsA[j]_M}) is created
#'  inside the passed \code{data} or \code{newdata} object while discretizing \code{sA[j]} into \code{M} bins.
#'
#' @docType class
#' @format An \code{\link{R6Class}} generator object
#' @keywords R6 class
#' @details
#' \itemize{
#' \item{\code{reg}} - .
#' \item{\code{outvar}} - .
#' \item{\code{intrvls}} - .
#' \item{\code{intrvls.width}} - .
#' \item{\code{bin_weights}} - .
#' }
#' @section Methods:
#' \describe{
#'   \item{\code{new(reg, data_object, DataStore.gstar, ...)}}{...}
#'   \item{\code{fit(data)}}{...}
#'   \item{\code{predict(newdata)}}{...}
#'   \item{\code{predictAeqa(newdata)}}{...}
#' }
#' @section Active Bindings:
#' \describe{
#'   \item{\code{cats}}{...}
#' }
#' @export
ContinSummaryModel <- R6Class(classname = "ContinSummaryModel",
  inherit = SummariesModel,
  portable = TRUE,
  class = TRUE,
  public = list(
    reg = NULL,
    outvar = character(),     # the name of the continous outcome var (sA[j])
    intrvls = NULL,
    intrvls.width = NULL,
    bin_weights = NULL,
    # Define settings for fitting contin sA and then call $new for super class (SummariesModel)
    initialize = function(reg, data_object, DataStore.gstar, ...) {
      self$reg <- reg
      self$outvar <- reg$outvar
      if (is.null(reg$intrvls)) {
        assert_that(is.DataStore(data_object))
        self$intrvls <- data_object$detect.sVar.intrvls(reg$outvar,
                                                      nbins = self$reg$nbins,
                                                      bin_bymass = self$reg$bin_bymass,
                                                      bin_bydhist = self$reg$bin_bydhist,
                                                      max_nperbin = self$reg$max_nperbin)
        if (!missing(DataStore.gstar)) {
          assert_that(is.DataStore(DataStore.gstar))
          gstar.intrvls <- DataStore.gstar$detect.sVar.intrvls(reg$outvar,
                                                      nbins = self$reg$nbins,
                                                      bin_bymass = self$reg$bin_bymass,
                                                      bin_bydhist = self$reg$bin_bydhist,
                                                      max_nperbin = self$reg$max_nperbin)
          self$intrvls <- unique(sort(union(self$intrvls, gstar.intrvls)))
        }
        # Define the number of bins (no. of binary regressions to run),
        # new outvar var names (bin names); all predvars remain unchanged;
        self$reg$intrvls <- self$intrvls
      } else {
        self$intrvls <- self$reg$intrvls
      }
      self$reg$nbins <- length(self$intrvls) - 1
      self$reg$bin_nms <- data_object$bin.nms.sVar(reg$outvar, self$reg$nbins)
      # Save bin widths in reg class (naming the vector entries by bin names):
      self$intrvls.width <- diff(self$intrvls)
      self$intrvls.width[self$intrvls.width <= gvars$tolerr] <- 1
      self$reg$intrvls.width <- self$intrvls.width
      names(self$reg$intrvls.width) <- names(self$intrvls.width) <- self$reg$bin_nms
      if (gvars$verbose)  {
        print("ContinSummaryModel outcome: "%+%self$outvar)
        # print("ContinSummaryModel reg$nbins: " %+% self$reg$nbins)
      }
      bin_regs <- def_regs_subset(self = self)
      super$initialize(reg = bin_regs, no_set_outvar = TRUE, ...)
    },

    # Transforms data for continous outcome to discretized bins sA[j] -> BinsA[1], ..., BinsA[M] and calls $super$fit on that transformed data
    # Gets passed redefined subsets that exclude degenerate Bins (prev subset is defined for names in sA - names have changed though)
    fit = function(data) {
      assert_that(is.DataStore(data))
      # Binirizes & saves binned matrix inside DataStore
      data$binirize.sVar(name.sVar = self$outvar, intervals = self$intrvls, nbins = self$reg$nbins, bin.nms = self$reg$bin_nms)
      if (gvars$verbose) {
        print("performing fitting for continuous outcome: " %+% self$outvar)
        print("freq counts by bin for continuous outcome: "); print(table(data$ord.sVar))
        print("binned dataset: "); print(head(cbind(data$ord.sVar, data$dat.bin.sVar), 5))
      }
      super$fit(data) # call the parent class fit method
      if (gvars$verbose) message("fit for outcome " %+% self$outvar %+% " succeeded...")
      data$emptydat.bin.sVar # wiping out binirized mat in data DataStore object...
      self$wipe.alldat # wiping out all data traces in ContinSummaryModel...
      invisible(self)
    },

    # P(A^s=1|W^s=w^s): uses private$m.fit to generate predictions
    predict = function(newdata) {
      if (missing(newdata)) {
        stop("must provide newdata")
      }
      assert_that(is.DataStore(newdata))
      if (gvars$verbose) print("performing prediction for continuous outcome: " %+% self$outvar)
      # mat_bin doesn't need to be saved (even though its invisibly returned); mat_bin is automatically saved in datnet.sW.sA - a potentially dangerous side-effect!!!
      newdata$binirize.sVar(name.sVar = self$outvar, intervals = self$intrvls, nbins = self$reg$nbins, bin.nms = self$reg$bin_nms)
      super$predict(newdata)
      newdata$emptydat.bin.sVar # wiping out binirized mat in newdata DataStore object...
      invisible(self)
    },

    # Convert contin. sA vector into matrix of binary cols, then call parent class method: super$predictAeqa()
    # Invisibly return cumm. prob P(sA=sa|sW=sw)
    predictAeqa = function(newdata) { # P(A^s=a^s|W^s=w^s) - calculating the likelihood for obsdat.sA[i] (n vector of a`s)
      assert_that(is.DataStore(newdata))
      newdata$binirize.sVar(name.sVar = self$outvar, intervals = self$intrvls, nbins = self$reg$nbins, bin.nms = self$reg$bin_nms)
      if (gvars$verbose) print("performing prediction for categorical outcome: " %+% self$outvar)
      bws <- newdata$get.sVar.bw(name.sVar = self$outvar, intervals = self$intrvls)
      self$bin_weights <- (1 / bws) # weight based on 1 / (sVar bin widths)
      # Option 1: ADJUST FINAL PROB by bw.j TO OBTAIN density at a point f(sa|sw) = P(sA=sa|sW=sw):
      cumprodAeqa <- super$predictAeqa(newdata = newdata) * self$bin_weights
      # Alternative 2: ALso integrate the difference of sA value and its left most bin cutoff: x - b_{j-1} and pass it
      # This is done so that we can integrate the constant hazard all the way to the value of x:
        # * (1 - bw.j.sA_diff*(1/self$bin_weights)*probA1) (discrete)
        # * exp(-bw.j.sA_diff*(1/self$bin_weights)*probA1) (continuous)
      # bw.j.sA_diff <- newdata$get.sVar.bwdiff(name.sVar = self$outvar, intervals = self$intrvls)
      # cumprodAeqa <- super$predictAeqa(newdata = newdata, bw.j.sA_diff = bw.j.sA_diff) * self$bin_weights
      newdata$emptydat.bin.sVar # wiping out binirized mat in newdata object...
      self$bin_weights <- NULL # wiping out self$bin_weights...
      self$wipe.alldat # wiping out all data traces in ContinSummaryModel...
      private$cumprodAeqa <- cumprodAeqa
      return(cumprodAeqa)
    },

    sampleA = function(newdata, ...) {
      assert_that(is.DataStore(newdata))
      # bring the sampled variable back to its original scale / levels:
      dat <- super$sampleA(newdata = newdata)
      rate <- 1
      # TODO: There should be a more elaborate way of sampling here, but it'll do for now
      # If the found value is in one of the tails, make sure the probability of it going to
      # A very distant value is small.
      dat.sampled <- sapply(dat, function(current_dat) {
        if (current_dat == 1) {
          self$intrvls[current_dat+1] -  rexp(length(current_dat), rate)
        } else if (current_dat == (length(self$intrvls) - 1)) {
          self$intrvls[current_dat] +  rexp(length(current_dat), rate)
        } else {
          runif(length(current_dat), self$intrvls[current_dat], self$intrvls[current_dat+1])
        }
      })
      return(dat.sampled)
    }

  ),
  active = list(
    cats = function() {seq_len(self$reg$nbins)}
  )
)
