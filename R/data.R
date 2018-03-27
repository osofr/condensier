#' Mock data set to test density estimation with weights
#'
#' A completely artificial data set from a simulation study designed to test the
#' functionality of the argument "weights" in the function "fit_density", as
#' implemented in pull request #12 to https://github.com/osofr/condensier
#'
#' @format A \code{data.frame} with 30,000 rows and 4 columns, including some
#'  observations marked with \code{NA} for missing data:
#' \describe{
#'   \item{W1}{A binary covariate with a multiplicative effect on the mean of
#'     the treatment. This background covariate does _not_ impact the censoring
#'     mechanism.}
#'   \item{W2}{A binary covariate with a multiplicative effect on the mean of
#'     the treatment. This background covariates directly affects the censoring
#'     mechanism.}
#'   \item{DeltaA}{A continuous value, the realization of an observation from a
#'     Gaussian distribution, with censored observations marked with \code{NA}.}
#'   \item{IPCW}{A _true_ inverse probability of censoring weight computed from
#'     the censoring mechanism used to induce missingness in the data.}
#' }
#
"test_data_weights"
