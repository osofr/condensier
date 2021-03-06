% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BinOutModelClass.R
\docType{class}
\name{BinDat}
\alias{BinDat}
\title{R6 class for storing the design matrix and binary outcome for a single logistic regression}
\format{An \code{\link{R6Class}} generator object}
\usage{
BinDat
}
\description{
This R6 class can request, store and manage the design matrix Xmat, as well as the binary outcome Bin for the
 logistic regression P(Bin|Xmat).
 Can also be used for converting data in wide format to long when requested,
 e.g., when pooling across binary indicators (fitting one pooled logistic regression model for several indicators)
 The class has methods that perform queries to data storage R6 class DataStore to get appropriate data columns & row subsets
}
\details{
\itemize{
\item{bin_names} - Names of the bins.
\item{ID} - Vector of observation IDs, \code{1:n}, used for pooling.
\item{pooled_bin_name} - Original name of the continuous covariate that was discretized into bins and then pooled.
\item{nbins} - Number of bins.
\item{outvar} - Outcome name.
\item{predvars} - Predictor names.
\item{pool} - Perform pooling of bins?
\item{outvars_to_pool} - Outcome bin indicators to pool?
\item{subset_expr} - Defines the subset which would be used for fitting this model (logical, expression or indices).
\item{subset_idx} - Subset \code{subset_expr} converted to logical vector.
}
}
\section{Methods}{

\describe{
  \item{\code{new(reg)}}{Uses \code{reg} R6 \code{\link{RegressionClass}} object to instantiate a new storage container for a
  design matrix and binary outcome.}
  \item{\code{show()}}{ Print information on outcome and predictor names used in this regression model}
  \item{\code{newdata()}}{...}
  \item{\code{define.subset_idx(...)}}{...}
  \item{\code{setdata()}}{...}
  \item{\code{setdata.long()}}{...}
}
}

\section{Active Bindings}{

\describe{
  \item{\code{emptydata}}{...}
  \item{\code{emptyY}}{...}
  \item{\code{emptySubset_idx}}{...}
  \item{\code{emptyN}}{...}
  \item{\code{getXmat}}{...}
  \item{\code{getY}}{...}
}
}

\keyword{R6}
\keyword{class}
