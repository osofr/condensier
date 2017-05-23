## ---------------------------------------------------------------------
#' R6 class for fitting and predicting joint probability for a univariate categorical summary measure sA[j]
#'
#' This R6 class defines and fits a conditional probability model \code{P(sA[j]|sW,...)} for a univariate
#'  categorical summary measure \code{sA[j]}. This class inherits from \code{\link{SummariesModel}} class.
#'  Defines the fitting algorithm for a regression model \code{sA[j] ~ sW + ...}.
#'  Reconstructs the likelihood \code{P(sA[j]=sa[j]|sW,...)} afterwards.
#'  Categorical \code{sA[j]} is first redefined into \code{length(levels)} bin indicator variables, where
#'  \code{levels} is a numeric vector of all unique categories in \code{sA[j]}.
#'  The fitting algorithm estimates the binary regressions for hazard for each bin indicator, \code{Bin_sA[j][i] ~ sW},
#'  i.e., the probability that categorical \code{sA[j]} falls into bin \code{i}, \code{Bin_sA[j]_i},
#'  given that \code{sA[j]} does not fall in any prior bins \code{Bin_sA[j]_1, ..., Bin_sA[j]_{i-1}}.
#'  The dataset of bin indicators (\code{BinsA[j]_1,...,BinsA[j]_M}) is created
#'  inside the passed \code{data} or \code{newdata} object when defining \code{length(levels)} bins for \code{sA[j]}.
#'
#' @docType class
#' @format An \code{\link{R6Class}} generator object
#' @keywords R6 class
#' @details
#' \itemize{
#' \item{\code{reg}} - .
#' \item{\code{outvar}} - .
#' \item{\code{levels}} - .
#' \item{\code{nbins}} - .
#' }
#' @section Methods:
#' \describe{
#'   \item{\code{new(reg, data_object, ...)}}{...}
#'   \item{\code{fit(data)}}{...}
#'   \item{\code{predict(newdata)}}{...}
#'   \item{\code{predictAeqa(newdata)}}{...}
#' }
#' @section Active Bindings:
#' \describe{
#'   \item{\code{cats}}{...}
#' }
#' @export
CategorSummaryModel <- R6Class(classname = "CategorSummaryModel",
  inherit = SummariesModel,
  portable = TRUE,
  class = TRUE,
  public = list(
    reg = NULL,
    outvar = character(),     # the name of the categorical outcome var (sA[j])
    levels = numeric(),       # all unique values for sA[j] sorted in increasing order
    nbins = integer(),
    # Define settings for fitting cat sA and then call $new for super class (SummariesModel)
    initialize = function(reg, data_object, ...) {
      self$reg <- reg
      self$outvar <- reg$outvar
      # Define the number of bins (no. of binary regressions to run) based on number of unique levels for categorical sVar:
      # all predvars remain unchanged
      if (is.null(reg$levels)) {
        assert_that(is.DataStore(data_object))
        self$levels <- self$reg$levels <- data_object$detect.cat.sVar.levels(reg$outvar)
      } else {
        self$levels <- self$reg$levels
      }
      self$nbins <- self$reg$nbins <- length(self$levels)
      self$reg$bin_nms <- data_object$bin.nms.sVar(reg$outvar, self$reg$nbins)
      if (gvars$verbose)  {
        print("CategorSummaryModel outcome: "%+%self$outvar)
        # print("CategorSummaryModel reg$levels: "); print(self$reg$levels)
        # print("CategorSummaryModel reg$nbins: " %+% self$reg$nbins)
      }
      bin_regs <- def_regs_subset(self = self)
      super$initialize(reg = bin_regs, no_set_outvar = TRUE, ...)
    },

    # Transforms data for categorical outcome to bin indicators sA[j] -> BinsA[1], ..., BinsA[M] and calls $super$fit on that transformed data
    # Gets passed redefined subsets that exclude degenerate Bins (prev subset is defined for names in sA - names have changed though)
    fit = function(data) {
      assert_that(is.DataStore(data))
      # Binirizes & saves binned matrix inside DataStore for categorical sVar
      data$binirize.cat.sVar(name.sVar = self$outvar, levels = self$levels)
      if (gvars$verbose) {
        print("performing fitting for categorical outcome: " %+% self$outvar)
        print("freq counts by bin for categorical outcome: "); print(table(data$get.sVar(self$outvar)))
        print("binned dataset: "); print(head(cbind(sA = data$get.sVar(self$outvar), data$dat.bin.sVar), 5))
      }
      super$fit(data) # call the parent class fit method
      if (gvars$verbose) message("fit for " %+% self$outvar %+% " var succeeded...")
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
      if (gvars$verbose) print("performing prediction for categorical outcome: " %+% self$outvar)
      newdata$binirize.cat.sVar(name.sVar = self$outvar, levels = self$levels)
      super$predict(newdata)
      newdata$emptydat.bin.sVar # wiping out binirized mat in newdata DataStore object...
      invisible(self)
    },

    # Invisibly return cumm. prob P(sA=sa|sW=sw)
    # P(A^s=a^s|W^s=w^s) - calculating the likelihood for obsdat.sA[i] (n vector of a's):
    predictAeqa = function(newdata) {
      assert_that(is.DataStore(newdata))
      if (gvars$verbose) print("performing prediction for categorical outcome: " %+% self$outvar)
      newdata$binirize.cat.sVar(name.sVar = self$outvar, levels = self$levels)
      cumprodAeqa <- super$predictAeqa(newdata = newdata)
      newdata$emptydat.bin.sVar # wiping out binirized mat in newdata object...
      self$wipe.alldat # wiping out all data traces in ContinSummaryModel...
      private$cumprodAeqa <- cumprodAeqa
      return(cumprodAeqa)
    },

    sampleA = function(newdata) {
      assert_that(is.DataStore(newdata))
      # bring the sampled variable back to its original scale / levels:
      sampleA <- self$levels[super$sampleA(newdata = newdata)]
      return(sampleA)
    }
  ),
  active = list(
    cats = function() {seq_len(self$reg$nbins)}
  )
)
