
#-----------------------------------------------------------------------------
# Global State Vars (can be controlled globally with options(condensier.optname = ))
#-----------------------------------------------------------------------------
gvars <- new.env(parent = emptyenv())
gvars$verbose <- TRUE      # verbose mode (print all messages)
gvars$opts <- list()        # named list of package options that is controllable by the user (condensier_options())
gvars$misval <- NA_integer_ # the default missing value for observations (# gvars$misval <- -.Machine$integer.max)
gvars$misXreplace <- 0L     # the default replacement value for misval that appear in the design matrix
gvars$tolerr <- 10^-12      # tolerance error: assume for abs(a-b) < gvars$tolerr => a = b
gvars$sVartypes <- list(bin = "binary", cat = "categor", cont = "contin")

#' Get Option Value for \code{condensier}
#'
#' @param optname The name of the consdensier option.
#' @return An option value.
#' @seealso \code{\link{condensier_options}}
#' @export
getopt <- function(optname) {
  opt <- gvars$opts
  if (!(optname %in% (names(opt)))) stop(optname %+% ": this options does not exist")
  return(opt[[optname]])
}

#' Print Current Option Settings for \code{condensier}
#' @return Invisibly returns a list of \code{condensier} options.
#' @seealso \code{\link{condensier_options}}
#' @export
print_condensier_opts <- function() {
  print(gvars$opts)
  invisible(gvars$opts)
}

#' (RETIRED) Setting Options for \code{condensier}
#'
#' This function is now retired. Please use \code{\link{fit_density}} directly for tuning parameter set-up.
#' Calling this function now will have no effect. Previously provided additional options that control the estimation algorithm in \code{condensier} package.
#' @param bin_estimator The estimator to use for fitting the binary outcomes (defaults to \code{speedglmR6} which estimates with \code{\link{speedglmR6}})
#'  another default option is \code{\link{glmR6}}.
#' @param bin_method The method for choosing bins when discretizing and fitting the conditional continuous summary
#'  exposure variable \code{sA}. The default method is \code{"equal.len"}, which partitions the range of \code{sA}
#'  into equal length \code{nbins} intervals. Method \code{"equal.mass"} results in a data-adaptive selection of the bins
#'  based on equal mass (equal number of observations), i.e., each bin is defined so that it contains an approximately
#'  the same number of observations across all bins. The maximum number of observations in each bin is controlled
#'  by parameter \code{max_n_bin}. Method \code{"dhist"} uses a mix of the above two approaches,
#'  see Denby and Mallows "Variations on the Histogram" (2009) for more detail.
#' @param parfit Default is \code{FALSE}. Set to \code{TRUE} to use \code{foreach} package and its functions
#'  \code{foreach} and \code{dopar} to perform
#'  parallel logistic regression fits and predictions for discretized continuous outcomes. This functionality
#'  requires registering a parallel backend prior to running \code{condensier} function, e.g.,
#'  using \code{doParallel} R package and running \code{registerDoParallel(cores = ncores)} for integer
#'  \code{ncores} parallel jobs. For an example, see a test in "./tests/RUnit/RUnit_tests_04_netcont_sA_tests.R".
#' @param nbins Set the default number of bins when discretizing a continous outcome variable under setting
#'  \code{bin_method = "equal.len"}.
#'  If left as \code{NA} the total number of equal intervals (bins) is determined by the nearest integer of
#'  \code{nobs}/\code{max_n_bin}, where \code{nobs} is the total number of observations in the input data.
#' @param max_n_cat Max number of unique categories a categorical variable \code{sA[j]} can have.
#' If \code{sA[j]} has more it is automatically considered continuous.
#' @param poolContinVar Set to \code{TRUE} for fitting a pooled regression which pools bin indicators across all bins.
#' When fitting a model for binirized continuous outcome, set to \code{TRUE}
#' for pooling bin indicators across several bins into one outcome regression?
#' @param max_n_bin Max number of observations per 1 bin for a continuous outcome (applies directly when
#'  \code{bin_method="equal.mass"} and indirectly when \code{bin_method="equal.len"}, but \code{nbins = NA}).
#' @return Invisibly returns a list with old option settings.
#' @seealso \code{\link{print_condensier_opts}}
#' @export
condensier_options <- function(bin_estimator = speedglmR6$new(),
                            parfit = FALSE,
                            bin_method = c("equal.len", "equal.mass", "dhist"),
                            nbins = NA,
                            max_n_cat = 20,
                            poolContinVar = FALSE,
                            max_n_bin = 1000
                            ) {
  warning("condensier_options() is now retired. Please use fit_density directly for tuning parameter set-up. Calling this function has no effect.")
  old.opts <- gvars$opts
  bin_method <- bin_method[1L]

  if (bin_method %in% "equal.len") {
  } else if (bin_method %in% "equal.mass") {
  } else if (bin_method %in% "dhist") {
  } else {
    stop("bin_method argument must be either 'equal.len', 'equal.mass' or 'dhist'")
  }

  opts <- list(
    bin_estimator = bin_estimator,
    bin_method = bin_method,
    parfit = parfit,
    nbins = nbins,
    max_n_cat = max_n_cat,
    poolContinVar = poolContinVar,
    max_n_bin = max_n_bin
  )
  # gvars$opts <- opts
  invisible(old.opts)
}

# returns a function (alternatively a call) that tests for missing values in (sA, sW)
testmisfun <- function() {
  if (is.na(gvars$misval)) {
    return(is.na)
  } else if (is.null(gvars$misval)){
    return(is.null)
  } else if (is.integer(gvars$misval)) {
    return(function(x) {x==gvars$misval})
  } else {
    return(function(x) {x%in%gvars$misval})
  }
}

get.misval <- function() {
  gvars$misfun <- testmisfun()
  gvars$misval
}

set.misval <- function(gvars, newmisval) {
  oldmisval <- gvars$misval
  gvars$misval <- newmisval
  gvars$misfun <- testmisfun()    # EVERYTIME gvars$misval HAS CHANGED THIS NEEDS TO BE RESET/RERUN.
  invisible(oldmisval)
}
gvars$misfun <- testmisfun()

# Allows condensier functions to use e.g., getOption("condensier.verbose") to get verbose printing status
.onLoad <- function(libname, pkgname) {
  op <- options()
  op.condensier <- list(
    condensier.verbose = gvars$verbose
  )
  # reset all options to their defaults on load:
  # condensier_options()

  toset <- !(names(op.condensier) %in% names(op))
  if(any(toset)) options(op.condensier[toset])

  invisible()
}

.onAttach <- function(...) {
  packageStartupMessage('condensier')
  packageStartupMessage('The condensier package is still in beta testing. Interpret results with caution.')
  #   packageStartupMessage('Version: ', utils::packageDescription('condensier')$Version)
  #   packageStartupMessage('Package created on ', utils::packageDescription('condensier')$Date, '\n')
  #   packageStartupMessage('Please note this package is still in its early stages of development.
   # Check for updates and report bugs at http://github.com/osofr/condensier.', '\n')
  #   packageStartupMessage('To see the vignette use vignette("condensier_vignette", package="condensier").
  # To see all available package documentation use help(package = "condensier") and ?condensier.', '\n')
  #   packageStartupMessage('To see the latest updates for this version, use news(package = "condensier").', '\n')
}











