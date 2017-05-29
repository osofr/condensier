#----------------------------------------------------------------------------------
# Classes that control modelling of the multivariate joint probability model P(sA|sW).
#----------------------------------------------------------------------------------
#' @importFrom assertthat assert_that

# S3 generic constructor for the summary model classes:
newsummarymodel <- function(reg, data_object, ...) { UseMethod("newsummarymodel") }
# Summary model constructor for generic regression with multivariate outcome, but one set of predictors
newsummarymodel.generic <- function(reg, data_object, ...) SummariesModel$new(reg = reg, data_object = data_object, ...)
# Summary model constructor for binary outcome sA[j]:
newsummarymodel.binary <- function(reg, ...) BinOutModel$new(reg = reg, ...)
# Summary model constructor for continuous outcome sA[j]:
newsummarymodel.contin <- function(reg, data_object, ...) ContinSummaryModel$new(reg = reg, data_object = data_object, ...)
# Summary model constructor for categorical outcome sA[j]:
newsummarymodel.categor <- function(reg, data_object, ...) CategorSummaryModel$new(reg = reg, data_object = data_object, ...)

## ---------------------------------------------------------------------
#' R6 class for fitting and predicting model P(sA|sW) under g.star or g.0
#'
#' This R6 class Class for defining, fitting and predicting the probability model
#'  \code{P(sA|sW)} under \code{g_star} or under \code{g_0} for summary measures
#'  (\code{sW,sA}). Defines and manages the factorization of the multivariate conditional
#'  probability model \code{P(sA=sa|...)} into univariate regression models
#'  \code{sA[j] ~ sA[j-1] + ... + sA[1] + sW}. The class \code{self$new} method automatically
#'  figures out the correct joint probability factorization into univariate conditional
#'  probabilities based on name ordering provided by (\code{sA_nms}, \code{sW_nms}).
#'  When the outcome variable \code{sA[j]} is binary, this class will automatically call
#'  a new instance of \code{\link{BinOutModel}} class.
#'  Provide \code{self$fit()} function argument \code{data} as a \code{\link{DataStore}} class object.
#'  This data will be used for fitting the model \code{P(sA|sW)}.
#'  Provide \code{self$fit()} function argument \code{newdata} (also as \code{DataStore} class) for predictions of the type
#'  \code{P(sA=1|sW=sw)}, where \code{sw} values are coming from \code{newdata} object.
#'  Finally, provide \code{self$predictAeqa} function \code{newdata} argument
#'  (also \code{DataStore} class object) for getting the likelihood predictions \code{P(sA=sa|sW=sw)}, where
#'  both, \code{sa} and \code{sw} values are coming from \code{newdata} object.
#'
#' @docType class
#' @format An \code{\link{R6Class}} generator object
#' @keywords R6 class
#' @details
#' \itemize{
#' \item{\code{n_regs}} - .
#' \item{\code{parfit_allowed}} - .
#' }
#' @section Methods:
#' \describe{
#'   \item{\code{new(reg, ...)}}{...}
#'   \item{\code{length}}{...}
#'   \item{\code{getPsAsW.models}}{...}
#'   \item{\code{getcumprodAeqa}}{...}
#'   \item{\code{copy.fit(SummariesModel)}}{...}
#'   \item{\code{fit(data)}}{...}
#'   \item{\code{predict(newdata)}}{...}
#'   \item{\code{predictAeqa(newdata, ...)}}{...}
#' }
#' @section Active Bindings:
#' \describe{
#'   \item{\code{wipe.alldat}}{...}
#' }
#' @export
SummariesModel <- R6Class(classname = "SummariesModel",
  portable = TRUE,
  class = TRUE,
  public = list(
    reg = NULL,
    outvar = character(),   # outcome name(s)
    predvars = character(), # names of predictor vars
    n_regs = integer(),        # total no. of reg. models (logistic regressions)
    parfit_allowed = FALSE,    # allow parallel fit of multivar outvar when 1) reg$parfit = TRUE & 2) all.outvar.bin = TRUE
    initialize = function(reg, no_set_outvar = FALSE, ...) {
      self$reg <- reg

      if (!no_set_outvar) self$outvar <- reg$outvar
      self$predvars <- reg$predvars

      self$n_regs <- length(reg$outvar) # Number of sep. logistic regressions to run
      all.outvar.bin <-  all(reg$outvar.class %in% gvars$sVartypes$bin)

      if (reg$parfit & all.outvar.bin & (self$n_regs > 1)) self$parfit_allowed <- TRUE
      #**** NOTE: for ltmle this should be changed to: if (reg$sep_predvars_sets) self$parfit_allowed <- TRUE

      if (gvars$verbose) {
        print("#----------------------------------------------------------------------------------")
        print("New instance of SummariesModel:")
        print("#----------------------------------------------------------------------------------")
        if (reg$sep_predvars_sets) {
          print("Outcomes: "); print(unlist(lapply(reg$outvar, function(outvars) "(" %+% paste(outvars, collapse = ", ") %+% ")")))
          print("Predictors: "); print(unlist(lapply(reg$predvars, function(predvars) "(" %+% paste(predvars, collapse = ", ") %+% ")")))
        } else {
          print("Outcomes: " %+% paste(reg$outvar, collapse = ", "))
          print("Predictors: " %+% paste(reg$predvars, collapse = ", "))
        }
        print("No. of regressions: " %+% self$n_regs)
        print("All outcomes binary? " %+% all.outvar.bin)
        if (self$parfit_allowed) print("Performing parallel fits: " %+% self$parfit_allowed)
        print("#----------------------------------------------------------------------------------")
      }

      # factorize the joint into univariate regressions, by dimensionality of the outcome variable (sA_nms):
      for (k_i in 1:self$n_regs) {
        reg_i <- reg$clone()
        reg_i$ChangeManyToOneRegresssion(k_i, reg)
        # Calling the constructor for the summary model P(sA[j]|\bar{sA}[j-1], sW}), dispatching on reg_i class
        PsAsW.model <- newsummarymodel(reg = reg_i, ...)
        private$PsAsW.models <- append(private$PsAsW.models, list(PsAsW.model))
        names(private$PsAsW.models)[k_i] <- "P(sA|sW)."%+%k_i
      }
      invisible(self)
    },

    length = function(){ base::length(private$PsAsW.models) },
    getPsAsW.models = function() { private$PsAsW.models },  # get all summary model objects (one model object per outcome var sA[j])
    getcumprodAeqa = function() { private$cumprodAeqa },  # get joint prob as a vector of the cumulative prod over j for P(sA[j]=a[j]|sW)


    update = function(newdata) {
      assert_that(is.DatNet.sWsA(newdata))
      #TODO: This should probably be combined with the fit function
      # serial loop over all regressions in PsAsW.models:
      if (!self$parfit_allowed) {
        for (k_i in seq_along(private$PsAsW.models)) {
          private$PsAsW.models[[k_i]]$update(newdata = newdata)
        }
      # parallel loop over all regressions in PsAsW.models:
      } else if (self$parfit_allowed) {
        val <- checkpkgs(pkgs=c("foreach", "doParallel", "matrixStats"))
        mcoptions <- list(preschedule = FALSE)
        # NOTE: Each fitRes[[k_i]] will contain a copy of every single R6 object that was passed by reference ->
        # *** the size of fitRes is 100x the size of private$PsAsW.models ***
        updateRes <- foreach::foreach(k_i = seq_along(private$PsAsW.models), .options.multicore = mcoptions) %dopar% {
          private$PsAsW.models[[k_i]]$update(newdata = newdata)
        }
        # copy the fits one by one from BinOutModels above into private field for BinOutModels
        for (k_i in seq_along(private$PsAsW.models)) {
          # The copy.fit and copy.update functions will probably be the exact same thing. Use the general
          # fit method here.
          private$PsAsW.models[[k_i]]$copy.fit(updateRes[[k_i]])
        }
      }
      invisible(self)
    },

    # The fit function gets called in the following way:
    # 1. from outside, starting the process.
    # 2. then the loop calls the inner summary models.
    # 3. These summmary models then call this function again, but they have a different list of submodels
    # 4. this propagates down, until super is not called anymore, and hence we are at the bottom of the tree
    #   the bin models.
    fit = function(data) {
      assert_that(is.DataStore(data))
      # serial loop over all regressions in PsAsW.models:
      if (!self$parfit_allowed) {
        for (k_i in seq_along(private$PsAsW.models)) {
          private$PsAsW.models[[k_i]]$fit(data = data)
        }
      # parallel loop over all regressions in PsAsW.models:
      } else if (self$parfit_allowed) {
        val <- checkpkgs(pkgs=c("foreach", "doParallel", "matrixStats"))
        mcoptions <- list(preschedule = FALSE)
        # NOTE: Each fitRes[[k_i]] will contain a copy of every single R6 object that was passed by reference ->
        # *** the size of fitRes is 100x the size of private$PsAsW.models ***
        fitRes <- foreach::foreach(k_i = seq_along(private$PsAsW.models), .options.multicore = mcoptions) %dopar% {
          private$PsAsW.models[[k_i]]$fit(data = data)
        }
        # copy the fits one by one from BinOutModels above into private field for BinOutModels
        for (k_i in seq_along(private$PsAsW.models)) {
          private$PsAsW.models[[k_i]]$copy.fit(fitRes[[k_i]])
        }
      }
      invisible(self)
    },

    # P(A^s=1|W^s=w^s): uses private$m.fit to generate predictions
    predict = function(newdata) {
      if (missing(newdata)) {
        stop("must provide newdata")
      }
      assert_that(is.DataStore(newdata))
      # serial loop over all regressions in PsAsW.models:
      if (!self$parfit_allowed) {
        for (k_i in seq_along(private$PsAsW.models)) {
          private$PsAsW.models[[k_i]]$predict(newdata = newdata)
        }
      # parallel loop over all regressions in PsAsW.models:
      } else if (self$parfit_allowed) {
        val <- checkpkgs(pkgs=c("foreach", "doParallel", "matrixStats"))
        mcoptions <- list(preschedule = FALSE)
        # NOTE: Each predRes[[k_i]] will contain a copy of every single R6 object that was passed by reference ->
        # *** the size of fitRes is 100x the size of private$PsAsW.models ***
        predRes <- foreach::foreach(k_i = seq_along(private$PsAsW.models), .options.multicore = mcoptions) %dopar% {
          private$PsAsW.models[[k_i]]$predict(newdata = newdata)
        }
        # copy the predictions one by one from BinOutModels above into private field for BinOutModels
        for (k_i in seq_along(private$PsAsW.models)) {
          private$PsAsW.models[[k_i]]$copy.predict(predRes[[k_i]])
        }
      }
      invisible(self)
    },

    # WARNING: This method cannot be chained together with other methods (s.a, class$predictAeqa()$fun())
    # Uses daughter objects (stored from prev call to fit()) to get predictions for P(sA=obsdat.sA|sW=sw)
    # Invisibly returns the joint probability P(sA=sa|sW=sw), also saves it as a private field "cumprodAeqa"
    # P(A^s=a^s|W^s=w^s) - calculating the likelihood for obsdat.sA[i] (n vector of a's):
    predictAeqa = function(newdata, ...) {
      assert_that(!missing(newdata))
      assert_that(is.DataStore(newdata))
      n <- newdata$nobs
      if (!self$parfit_allowed) {
        cumprodAeqa <- rep.int(1, n)
        # loop over all regressions in PsAsW.models:
        for (k_i in seq_along(private$PsAsW.models)) {
          cumprodAeqa <- cumprodAeqa * private$PsAsW.models[[k_i]]$predictAeqa(newdata = newdata, ...)
        }
      } else if (self$parfit_allowed) {
        val <- checkpkgs(pkgs=c("foreach", "doParallel", "matrixStats"))
        mcoptions <- list(preschedule = TRUE)
        probAeqa_list <- foreach::foreach(k_i = seq_along(private$PsAsW.models), .options.multicore = mcoptions) %dopar% {
          private$PsAsW.models[[k_i]]$predictAeqa(newdata = newdata, ...)
        }
        probAeqa_mat <- do.call('cbind', probAeqa_list)
        cumprodAeqa <- matrixStats::rowProds(probAeqa_mat)
      }
      private$cumprodAeqa <- cumprodAeqa
      return(cumprodAeqa)
    },

    ## Sample new outcome (A), given predictors (X) in newdata
    ## The outcome can be continuous / categorical
    sampleA = function(newdata, ...) {

      assert_that(!missing(newdata))
      assert_that(is.DataStore(newdata))
      n <- newdata$nobs

      # loop over all regressions in PsAsW.models, sample CONDITIONALLY on observations that haven't been put in a specific bin yet
      sampleA_mat <- matrix(0L, nrow = n, ncol = length(private$PsAsW.models))

      for (k_i in seq_along(private$PsAsW.models)) {
        sampleA_newcat <- private$PsAsW.models[[k_i]]$sampleA(newdata = newdata, ...)
        if (k_i == 1L) sampleA_mat[, k_i] <- sampleA_newcat
        # carry forward all previously sampled 1's (degenerate ones a bin a chosen for the first time)
        if (k_i > 1) {
          # if you succeeded at the previous bin, your 1L is carried through till the end:
          sampleA_mat[(sampleA_mat[, k_i - 1] == 1L), k_i] <- 1L
          # if you haven't succeeded at the previous bin, you get a chance to succeed at this category:
          sampleA_mat[(sampleA_mat[, k_i - 1] == 0L), k_i] <- sampleA_newcat[(sampleA_mat[, k_i - 1] == 0L)]
        }
      }

      if (length(private$PsAsW.models) > 1) {
        sampleA_mat[, length(private$PsAsW.models)] <- 1L # make last category a reference category
        sampleA_cat <- rowSums(1L - sampleA_mat) + 1L
      } else {
        sampleA_cat <- as.vector(sampleA_mat)
      }

      return(sampleA_cat)
    }
  ),

  active = list(
    # recursively call all saved daughter model fits and wipe out any traces of saved data
    wipe.alldat = function() {
      for (k_i in seq_along(private$PsAsW.models)) {
        private$PsAsW.models[[k_i]]$wipe.alldat
      }
      return(self)
    }
  ),

  private = list(
    deep_clone = function(name, value) {
      # if value is is an environment, quick way to copy:
      # list2env(as.list.environment(value, all.names = TRUE), parent = emptyenv())
      # if a list of R6 objects, make a deep copy of each:
      if (name == "PsAsW.models") {
        lapply(value, function(PsAsW.model) PsAsW.model$clone(deep=TRUE))
      # to check the value is an R6 object:
      } else if (inherits(value, "R6")) {
        value$clone(deep=TRUE)
      } else {
        # For all other fields, just return the value
        value
      }
    },
    PsAsW.models = list(),
    fitted.pbins = list(),
    cumprodAeqa = NULL
  )
)
