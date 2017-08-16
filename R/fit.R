
## ---------------------------------------------------------------------------------------
#' Fit (multivariate) conditional density
#'
#' @param X A vector containing the names of predictor variables to use for modeling.
#' @param Y A character string name of the column that represent the response variable(s) in the model.
#' @param input_data Input dataset, can be a \code{data.frame} or a \code{data.table}.
#' @param bin_estimator The estimator to use for fitting the binary outcomes (defaults to \code{speedglmR6} which estimates with \code{\link{speedglmR6}})
#'  another default option is \code{\link{glmR6}}.
#' @param bin.method The method for choosing bins when discretizing and fitting the conditional continuous summary
#'  exposure variable \code{sA}. The default method is \code{"equal.len"}, which partitions the range of \code{sA}
#'  into equal length \code{nbins} intervals. Method \code{"equal.mass"} results in a data-adaptive selection of the bins
#'  based on equal mass (equal number of observations), i.e., each bin is defined so that it contains an approximately
#'  the same number of observations across all bins. The maximum number of observations in each bin is controlled
#'  by parameter \code{maxNperBin}. Method \code{"dhist"} uses a mix of the above two approaches,
#'  see Denby and Mallows "Variations on the Histogram" (2009) for more detail.
#' @param parfit Default is \code{FALSE}. Set to \code{TRUE} to use \code{foreach} package and its functions
#'  \code{foreach} and \code{dopar} to perform
#'  parallel logistic regression fits and predictions for discretized continuous outcomes. This functionality
#'  requires registering a parallel backend prior to running \code{condensier} function, e.g.,
#'  using \code{doParallel} R package and running \code{registerDoParallel(cores = ncores)} for integer
#'  \code{ncores} parallel jobs. For an example, see a test in "./tests/RUnit/RUnit_tests_04_netcont_sA_tests.R".
#' @param nbins Set the default number of bins when discretizing a continous outcome variable under setting
#'  \code{bin.method = "equal.len"}.
#'  If left as \code{NA} the total number of equal intervals (bins) is determined by the nearest integer of
#'  \code{nobs}/\code{maxNperBin}, where \code{nobs} is the total number of observations in the input data.
#' @param maxncats Max number of unique levels for cat outcome.
#' If the outcome has more levels it is automatically considered continuous.
#' @param poolContinVar Set to \code{TRUE} for fitting a pooled regression which pools bin indicators across all bins.
#' When fitting a model for binirized continuous outcome, set to \code{TRUE}
#' for pooling bin indicators across several bins into one outcome regression?
#' @param maxNperBin Max number of observations per 1 bin for a continuous outcome (applies directly when
#'  \code{bin.method="equal.mass"} and indirectly when \code{bin.method="equal.len"}, but \code{nbins = NA}).
#' @param verbose Set to \code{TRUE} to print messages on status and information to the console.
#' Turn this on by default using \code{options(condensier.verbose=TRUE)}.
#'
#' @section Details:
#' **********************************************************************
#'
#' \itemize{
#' \item As described in the following section, the first step is to construct an estimator \eqn{P_{g_N}(sA | sW)}
#'    for the common (in \code{i}) conditional density \eqn{P_{g_0}(sA | sW)} for common (in \code{i}) unit-level summary
#'    measures (\code{sA},\code{sW}).
#'
#' \item The same fitting algorithm is applied to construct an estimator \eqn{P_{g^*_N}(sA^* | sW^*)} of the common (in \code{i})
#'    conditional density \eqn{P_{g^*}(sA^* | sW^*)} for common (in \code{i}) unit-level summary measures (\code{sA^*},\code{sW^*})
#'    implied by the user-supplied stochastic intervention \code{f_gstar1} or \code{f_gstar2} and the observed distribution of \code{W}.
#'
#' \item These two density estimators form the basis of the IPTW estimator,
#'    which is evaluated at the observed N data points \eqn{O_i=(sW_i, sA_i, Y_i), i=1,...,N} and is given by
#'    \deqn{\psi^{IPTW}_n = \sum_{i=1,...,N}{Y_i \frac{P_{g^*_N}(sA^*=sA_i | sW=sW_i)}{P_{g_N}(sA=sA_i | sW=sW_i)}}.}
#' }
#'
#' @section Modeling \code{P(sA|sW)} for summary measures \code{(sA,sW)}:
#' **********************************************************************
#'
#' Non-parametric
#'  estimation of the common \strong{unit-level} multivariate joint conditional probability model \code{P_g0(sA|sW)},
#'  for unit-level summary measures \code{(sA,sW)} generated from the observed exposures and baseline covariates
#'  \eqn{(A,W)=(A_i,W_i : i=1,...,N)} (their joint density given by \eqn{g_0(A|W)Q(W)}), is performed by first
#'  constructing the dataset of N summary measures, \eqn{(sA_i,sW_i : i=1,...,N)}, and then fitting the usual i.i.d. MLE
#'  for the common density \code{P_g0(sA|sW)} based on the pooled N sample of these summary measures.
#'
#'  Note that \code{sA} can be multivariate and any of its components \code{sA[j]} can be either binary, categorical
#'  or continuous.
#'  The joint probability model for \code{P(sA|sA)} = \code{P(sA[1],...,sA[k]|sA)} can be factorized as
#'  \code{P(sA[1]|sA)} * \code{P(sA[2]|sA, sA[1])} * ... * \code{P(sA[k]|sA, sA[1],...,sA[k-1])},
#'  where each of these conditional probability models is fit separately, depending on the type of the outcome variable
#'  \code{sA[j]}.
#'
#'  If \code{sA[j]} is binary, the conditional probability \code{P(sA[j]|sW,sA[1],...,sA[j-1])} is evaluated via logistic
#'  regression model.
#'  When \code{sA[j]} is continuous (or categorical), its range will be fist partitioned into \code{K} bins and the
#'  corresponding \code{K}
#'  bin indicators (\code{B_1,...,B_K}), where each bin indicator \code{B_j} is then used as an outcome in a
#'  separate logistic regression model with predictors given by \code{sW, sA[1],...,sA[k-1]}.
#'  Thus, the joint probability \code{P(sA|sW)} is defined by such a tree of binary logistic regressions.
#'
#' For simplicity, we now suppose \code{sA} is continuous and univariate and we describe here an algorithm for fitting
#'  \eqn{P_{g_0}(sA | sW)} (the algorithm
#'  for fitting \eqn{P_{g^*}(sA^* | sW^*)} is equivalent, except that exposure \code{A} is replaced with exposure \code{A^*}
#'  generated under \code{f_gstar1} or \code{f_gstar2} and
#'  the predictors \code{sW} from the regression formula \code{hform.g0} are replaced with predictors \code{sW^*}
#'  specified by the regression formula \code{hform.gstar}).
#'
#' \enumerate{
#' \item Generate a dataset of N observed continuous summary measures (\code{sA_i}:i=1,...,N) from observed
#'  ((\code{A_i},\code{W_i}):i=1,...,N).
#'
#' \item Divide the range of \code{sA} values into intervals S=(i_1,...,i_M,i_{M+1}) so that any observed data point
#'    \code{sa_i} belongs to one interval in S, namely,
#'    for each possible value sa of \code{sA} there is k\\in{1,...,M}, such that, i_k < \code{sa} <= i_{k+1}.
#'    Let the mapping B(sa)\\in{1,...,M} denote a unique interval in S for sa, such that, i_{B(sa)} < sa <= i_{B(sa)+1}.
#'    Let bw_{B(sa)}:=i_{B(sa)+1}-i_{B(sa)} be the length of the interval (bandwidth) (i_{B(sa)},i_{B(sa)+1}).
#'    Also define the binary indicators b_1,...,b_M, where b_j:=I(B(sa)=j), for all j <= B(sa) and b_j:=NA for all j>B(sa).
#'    That is we set b_j to missing ones the indicator I(B(sa)=j) jumps from 0 to 1.
#'    Now let \code{sA} denote the random variable for the observed summary measure for one unit
#'    and denote by (B_1,...,B_M) the corresponding random indicators for \code{sA} defined as B_j := I(B(\code{sA}) = j)
#'    for all j <= B(\code{sA}) and B_j:=NA for all j>B(\code{sA}).
#'
#' \item For each j=1,...,M, fit the logistic regression model for the conditional probability P(B_j = 1 | B_{j-1}=0, sW), i.e.,
#'    at each j this is defined as the conditional probability of B_j jumping from 0 to 1 at bin j, given that B_{j-1}=0 and
#'    each of these logistic regression models is fit only among the observations that are still at risk of having B_j=1 with B_{j-1}=0.
#'
#' \item Normalize the above conditional probability of B_j jumping from 0 to 1 by its corresponding interval length (bandwidth) bw_j to
#'    obtain the discrete conditional hazards h_j(sW):=P(B_j = 1 | (B_{j-1}=0, sW) / bw_j, for each j.
#'    For the summary measure \code{sA}, the above conditional hazard h_j(sW) is equal to P(\code{sA} \\in (i_j,i_{j+1}) | \code{sA}>=i_j, sW),
#'    i.e., this is the probability that \code{sA} falls in the interval (i_j,i_{j+1}), conditional on sW and conditional on the fact that
#'    \code{sA} does not belong to any intervals before j.
#'
#' \item  Finally, for any given data-point \code{(sa,sw)}, evaluate the discretized conditional density for P(\code{sA}=sa|sW=sw) by first
#'    evaluating the interval number k=B(sa)\\in{1,...,M} for \code{sa} and then computing \\prod{j=1,...,k-1}{1-h_j(sW))*h_k(sW)}
#'    which is equivalent to the joint conditional probability that \code{sa} belongs to the interval (i_k,i_{k+1}) and does not belong
#'    to any of the intervals 1 to k-1, conditional on sW.
#'  }
#'
#' The evaluation above utilizes a discretization of the fact that any continuous density f of random variable X can be written as f_X(x)=S_X(x)*h_X(x),
#'  for a continuous density f of X where S_X(x):=P(X>x) is the survival function for X, h_X=P(X>x|X>=x) is the hazard function for X; as well as the fact that
#'  the discretized survival function S_X(x) can be written as a of the hazards for s<x: S_X(x)=\\prod{s<x}h_X(x).
#'
#' @section Three methods for defining bin (interval) cuttoffs for a continuous one-dimenstional summary measure \code{sA[j]}:
#' **********************************************************************
#'
#' There are 3 alternative methods to defining the bin cutoffs S=(i_1,...,i_M,i_{M+1}) for a continuous summary measure
#'  \code{sA}. The choice of which method is used along with other discretization parameters (e.g., total number of
#'  bins) is controlled via the tmlenet_options() function. See \code{?tmlenet_options} argument \code{bin.method} for
#'  additional details.
#'
#' Approach 1 (\code{equal.len}): equal length, default.
#'
#' *********************
#'
#' The bins are defined by splitting the range of observed \code{sA} (sa_1,...,sa_n) into equal length intervals.
#'  This is the dafault discretization method, set by passing an argument \code{bin.method="equal.len"} to
#'  \code{tmlenet_options} function prior to calling \code{tmlenet()}. The intervals will be defined by splitting the
#'  range of (sa_1,...,sa_N) into \code{nbins} number of equal length intervals, where \code{nbins} is another argument
#'  of \code{tmlenet_options()} function. When \code{nbins=NA} (the default setting) the actual value of \code{nbins}
#'  is computed at run time by taking the integer value (floor) of \code{n/maxNperBin},
#'  for \code{n} - the total observed sample size and \code{maxNperBin=1000} - another argument of
#'  \code{tmlenet_options()} with the default value 1,000.
#'
#' Approach 2 (\code{equal.mass}): data-adaptive equal mass intervals.
#'
#' *********************
#'
#' The intervals are defined by splitting the range of \code{sA} into non-equal length data-adaptive intervals that
#'  ensures that each interval contains around
#'  \code{maxNperBin} observations from (sa_j:j=1,...,N).
#'  This interval definition approach can be selected by passing an argument \code{bin.method="equal.mass"} to
#'  \code{tmlenet_options()} prior to calling \code{tmlenet()}.
#'  The method ensures that an approximately equal number of observations will belong to each interval, where that number
#'  of observations for each interval
#'  is controlled by setting \code{maxNperBin}. The default setting is \code{maxNperBin=1000} observations per interval.
#'
#' Approach 3 (\code{dhist}): combination of 1 & 2.
#'
#' *********************
#'
#' The data-adaptive approach dhist is a mix of Approaches 1 & 2. See Denby and Mallows "Variations on the Histogram"
#'  (2009)). This interval definition method is selected by passing an argument \code{bin.method="dhist"} to
#'  \code{tmlenet_options()}  prior to calling \code{tmlenet()}.
#'
#' @return An R6 object of class \code{SummariesModel} containing the conditional density fit(s).
#' @example tests/examples/1_condensier_example.R
#' @export
fit_density <- function(
                      X,
                      Y,
                      input_data,
                      bin.method = c("equal.mass", "equal.len", "dhist"),
                      nbins = NA_integer_,
                      maxncats = 20,
                      poolContinVar = FALSE,
                      maxNperBin = 1000,
                      parfit = FALSE,
                      bin_estimator = speedglmR6$new(),
                      verbose = getOption("condensier.verbose")
                      ) {

  curr.gvars <- gvars$verbose
  gvars$verbose <- verbose

  bin.method <- bin.method[1L]
  if (bin.method %in% "equal.len") {
  } else if (bin.method %in% "equal.mass") {
  } else if (bin.method %in% "dhist") {
  } else { stop("bin.method argument must be either 'equal.len', 'equal.mass' or 'dhist'") }

  nbins <- nbins[1L]

  # opts <- list(
  #   bin_estimator = bin_estimator,
  #   bin.method = bin.method,
  #   parfit = parfit,
  #   nbins = nbins,
  #   maxncats = maxncats,
  #   poolContinVar = poolContinVar,
  #   maxNperBin = maxNperBin
  # )
  # gvars$opts <- opts

  if (!is.data.table(input_data)) data.table::setDT(input_data)

  ## import the input data into internal storage class
  data_store_obj <- DataStore$new(input_data = input_data, Y = Y, X = X, maxncats = maxncats)

  # Find the class of the provided variable
  outcome.class <- data_store_obj$type.sVar[Y]

  # subset_vars <- lapply(Y, function(var) {var})

  # Put all est_params in RegressionClass
  regclass <- RegressionClass$new(bin_estimator = bin_estimator,
                                  nbins = nbins,
                                  outvar.class = outcome.class,
                                  outvar = Y,
                                  predvars = X,
                                  parfit = parfit,
                                  bin_bymass = bin.method%in%"equal.mass",
                                  bin_bydhist = bin.method%in%"dhist",
                                  max_nperbin = maxNperBin,
                                  pool_cont = poolContinVar)
                                  # subset = subset_vars)

  # Create the conditional density, based on the regression just specified and fit it
  conditional_density <- SummariesModel$new(reg = regclass, data_object = data_store_obj)

  # print("conditional_density$reg$subset"); print(conditional_density$reg$subset)
  conditional_density$fit(data = data_store_obj)

  gvars$verbose <- curr.gvars
  return(conditional_density)
}

## ---------------------------------------------------------------------------------------
#' Predict probability (likelihood) for cond. density fit for new data
#'
#' @param model_fit An R6 object of class \code{SummariesModel} returned by \code{\link{fit_density}}.
#' @param newdata A \code{data.table} or \code{data.frame} containing new predictors and IMPORTANTLY the outcomes.
#' @return A numeric vector containing the likelihood predictions for new observation.
#' When \code{bin.method = "equal.mass"} the output is a named vector, with the names corresponding to the
#' specific bin percentile of each individual outcome.
#' @example tests/examples/1_condensier_example.R
#' @export
predict_probability <- function(model_fit, newdata) {
  assert_that(is(model_fit, "SummariesModel"))
  newdata_obj <- DataStore$new(newdata,
                               Y = model_fit$reg$outvar,
                               X = model_fit$reg$predvars,
                               auto_typing = FALSE)

  return(model_fit$predictAeqa(newdata = newdata_obj))
}

## ---------------------------------------------------------------------------------------
#' Sample values from the existing conditional density fit for given new data
#'
#' @param model_fit An R6 object of class \code{SummariesModel} returned by \code{\link{fit_density}}.
#' @param newdata A \code{data.table} or \code{data.frame} containing new predictors / covariates, the outcomes are not needed and wont be used.
#' @return A numeric vector containing the sampled predictions for new observation.
#' @example tests/examples/1_condensier_example.R
#' @export
sample_value <- function(model_fit, newdata) {
  assert_that(is(model_fit, "SummariesModel"))
  newdata_obj <- DataStore$new(newdata,
                               Y = model_fit$reg$outvar,
                               X = model_fit$reg$predvars,
                               auto_typing = FALSE)
  return(model_fit$sampleA(newdata = newdata_obj))
}


