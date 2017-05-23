
## ---------------------------------------------------------------------------------------
#' Fit (multivariate) conditional density
#'
#' @param X A vector containing the names of predictor variables to use for modeling.
#' @param Y A character string name of the column that represent the response variable(s) in the model.
#' @param input_data Input dataset, can be a \code{data.frame} or a \code{data.table}.
#' @param nbins ...
#' @param bin_estimator ....
#' @param verbose Set to \code{TRUE} to print messages on status and information to the console.
#' Turn this on by default using \code{options(condensier.verbose=TRUE)}.
#' @return An R6 object containing the model fit(s).
#' @export
fit_density <- function(
                      X,
                      Y,
                      input_data,
                      nbins = getopt("nbins"),
                      bin_estimator = getopt("bin_estimator"),
                      verbose = getOption("condensier.verbose")
                      ) {

  # gvars$verbose <- verbose
  if (!is.data.table(input_data)) data.table::setDT(input_data)

  ## import the input data into internal storage class
  data_store_obj <- DataStore$new(input_data = data.table::data.table(datO), Y = "sA", X = c("W1", "W2", "W3"))

  # Find the class of the provided variable
  outcome.class <- data_store_obj$type.sVar[Y]

  subset_vars <- lapply(Y, function(var) {var})

  # Put all est_params in RegressionClass
  regclass <- RegressionClass$new(bin_estimator = bin_estimator,
                                  nbins = nbins,
                                  outvar.class = outcome.class,
                                  outvar = Y,
                                  predvars = X,
                                  subset = subset_vars)

  # Create the conditional density, based on the regression just specified and fit it
  conditional_density <- SummariesModel$new(reg = regclass, data_object = data_store_obj)

  conditional_density$fit(data = data_store_obj)

  return(conditional_density)
}

#' @export
predict_probability <- function(model_fit, newdata) {
  assert_that(is(model_fit, "SummariesModel"))
  newdata_obj <- DataStore$new(newdata,
                               Y = model_fit$reg$outvar,
                               X = model_fit$reg$predvars,
                               auto_typing = FALSE)
  return(model_fit$predictAeqa(newdata = newdata_obj))
}

#' @export
sample_value <- function(model_fit, newdata) {
  assert_that(is(model_fit, "SummariesModel"))
  newdata_obj <- DataStore$new(newdata,
                               Y = model_fit$reg$outvar,
                               X = model_fit$reg$predvars,
                               auto_typing = FALSE)
  return(model_fit$sampleA(newdata = newdata_obj))
}


