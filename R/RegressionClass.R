
## ---------------------------------------------------------------------
#' R6 class that defines regression models evaluating P(sA|sW), for summary measures (sW,sA)
#'
#' This R6 class defines fields and methods that controls all the parameters for non-parametric
#'  modeling and estimation of multivariate joint conditional probability model \code{P(sA|sW)} for summary measures \code{(sA,sW)}.
#'  Note that \code{sA} can be multivariate and any component of \code{sA[j]} can be either binary, categorical or continuous.
#'  The joint probability for \code{P(sA|sA)} = \code{P(sA[1],...,sA[k]|sA)} is first factorized as
#'  \code{P(sA[1]|sA)} * \code{P(sA[2]|sA, sA[1])} * ... * \code{P(sA[k]|sA, sA[1],...,sA[k-1])},
#'  where each of these conditional probability models is defined by a new instance of a \code{\link{SummariesModel}} class
#'  (and a corresponding instance of the \code{RegressionClass} class).
#'  If \code{sA[j]} is binary, the conditional probability \code{P(sA[j]|sW,sA[1],...,sA[j-1])} is evaluated via logistic regression model.
#'  When \code{sA[j]} is continuous (or categorical), its estimation will be controlled by a new instance of
#'  the \code{\link{ContinSummaryModel}} class (or the \code{\link{CategorSummaryModel}} class), as well as the accompanying new instance of the
#'  \code{RegressionClass} class. The range of continuous \code{sA[j]} will be fist partitioned into \code{K} bins and the corresponding \code{K}
#'  bin indicators (\code{B_1,...,B_K}), with \code{K} new instances of \code{\link{SummariesModel}} class, each instance defining a
#'  single logistic regression model for one binary bin indicator outcome \code{B_j} and predictors (\code{sW, sA[1],...,sA[k-1]}).
#'  Thus, the first instance of \code{RegressionClass} and \code{SummariesModel} classes will automatically
#'  spawn recursive calls to new instances of these classes until the entire tree of binary logistic regressions that defines
#'  the joint probability \code{P(sA|sW)} is build.
#'
#' @docType class
#' @format An \code{\link{R6Class}} generator object
#' @keywords R6 class
#' @details
#' \itemize{
#' \item{\code{sep_predvars_sets}} - Logical indicating the type of regression to run,
#'    if \code{TRUE} fit the joint P(\code{outvar}|\code{predvars}) (default),
# '   if \code{FALSE}, fit P(\code{outvar[1]}|\code{predvars[[1]]})*...*P(\code{outvar[K]}|\code{predvars[[K]}].
#'    More specifically, if \code{FALSE} (default), use the same predictors in \code{predvars} (vector of names) for all nodes in \code{outvar};
#'    when \code{TRUE} uses separate sets in \code{predvars} (must be a named list of character vectors) for fitting each node in \code{outvar}.
#' \item{\code{outvar.class}} - Character vector indicating a class of each outcome var: \code{bin} / \code{cont} / \code{cat}.
#' \item{\code{outvar}} - Character vector of regression outcome variable names.
#' \item{\code{predvars}} - Either a pooled character vector of all predictors (\code{sW}) or a vector of regression-specific predictor names.
#'      When \code{sep_predvars_sets=TRUE}, this must be a named list of predictor names, the list names corresponding to each node name in \code{outvar},
#'      and each list item being a vector specifying the regression predictors for a specific outcome in \code{outvar}.
#' \item{{reg_hazard}} - Logical, if TRUE, the joint probability model P(outvar | predvars) is factorized as
#'    \\prod_{j}{P(outvar[j] | predvars)} for each j outvar (for fitting hazard).
#' \item{\code{subset}} - Subset expression (later evaluated to logical vector in the envir of the data).
#' \item{\code{ReplMisVal0}} - Logical, if TRUE all gvars$misval among predicators are replaced with with gvars$misXreplace (0).
#' \item{\code{nbins}} - Integer number of bins used for a continuous outvar, the intervals are defined inside
#'  \code{ContinSummaryModel$new()} and then saved in this field.
#' \item{\code{bin_nms}} - Character vector of column names for bin indicators.
#' \item{\code{bin_estimator}} - A subclass of \code{\link{logisfitR6}}. This is the algorithm used to fit the bins. The
#'    default options are \code{\link{glmR6}} (fit the logistic regression model using \code{glm}), and \code{\link{speedglmR6}},
#'    with which we use use \code{speedglm}.
#' \item{\code{parfit}} - Logical, if TRUE then use parallel \code{foreach::foreach} loop to fit and predict binary logistic
#'    regressions (requires registering back-end cluster prior to calling the fit/predict functions)..
#' \item{\code{bin_bymass}} - Logical, for continuous outvar, create bin cutoffs based on equal mass distribution.
#' \item{\code{bin_bydhist}} - Logical, if TRUE, use dhist approach for bin definitions.  See Denby and Mallows "Variations on the
#'    Histogram" (2009)) for more..
#' \item{\code{max_nperbin}} - Integer, maximum number of observations allowed per one bin.
#' \item{\code{pool}} - Logical, pool binned continuous outvar observations across bins and only fit only regression model
#'    across all bins (adding bin_ID as an extra covaraite)..
#' \item{\code{outvars_to_pool}} - Character vector of names of the binned continuous outvars, should match \code{bin_nms}.
#' \item{\code{intrvls.width}} - Named numeric vector of bin-widths (\code{bw_j : j=1,...,M}) for each each bin in \code{self$intrvls}.
#'    When \code{sA} is not continuous, \code{intrvls.width} IS SET TO 1. When sA is continuous and this variable \code{intrvls.width}
#'    is not here, the intervals are determined inside \code{ContinSummaryModel$new()} and are assigned to this variable as a list,
#'    with \code{names(intrvls.width) <- reg$bin_nms}. Can be queried by \code{BinOutModel$predictAeqa()} as: \code{intrvls.width[outvar]}.
#' \item{\code{intrvls}} - Numeric vector of cutoffs defining the bins or a named list of numeric intervals for \code{length(self$outvar) > 1}.
#' \item{\code{cat.levels}} - Numeric vector of all unique values in categorical outcome variable.
#'    Set by \code{\link{CategorSummaryModel}} constructor.
#' }
#' @section Methods:
#' \describe{
#'   \item{\code{new(sep_predvars_sets = FALSE,
#'                   outvar.class = gvars$sVartypes$bin,
#'                   outvar, predvars, subset, intrvls,
#'                   ReplMisVal0 = TRUE,
#'                   bin_estimator = getopt("bin_estimator"),
#'                   parfit = getopt("parfit"),
#'                   nbins = getopt("nbins"),
#'                   bin_bymass = getopt("bin_method")%in%"equal.mass",
#'                   bin_bydhist = getopt("bin_method")%in%"dhist",
#'                   max_nperbin = getopt("max_n_bin"),
#'                   pool = getopt("poolContinVar")}}{Uses the arguments to instantiate an object of R6 class and define the future regression model.}
#'   \item{\code{ChangeManyToOneRegresssion(k_i, reg)}}{ Take a clone of a parent \code{RegressionClass} (\code{reg}) for \code{length(self$outvar)} regressions
#'    and set self to a single univariate \code{k_i} regression for outcome \code{self$outvar[[k_i]]}.}
#'   \item{\code{ChangeOneToManyRegresssions(regs_list)}}{ Take the clone of a parent \code{RegressionClass} for univariate (continuous outvar) regression
#'     and set self to \code{length(regs_list)} bin indicator outcome regressions.}
#'   \item{\code{resetS3class()}}{...}
#' }
#' @section Active Bindings:
#' \describe{
#'   \item{\code{S3class}}{...}
#'   \item{\code{get.reg}}{...}
#' }
#' @export
RegressionClass <- R6Class("RegressionClass",
  class = TRUE,
  portable = TRUE,
  public = list(
    sep_predvars_sets = logical(), # The type of regression to run, either fitting the joint P(outvar|predvars) (default) or P(outvar[1]|predvars[[1]])*...*P(outvar[K]|predvars[[K]]
                                    # if FALSE (default), use the same predictors in predvars (vector of names) for all nodes in outvar;
                                    # if TRUE use separate sets in predvars (named list of character vectors) for fitting each node in outvar.
    outvar.class = character(),    # vector for classes of the outcome vars: bin / cont / cat
    outvar = character(),          # vector of regression outcome variable names
    predvars = character(),        # either a:
                                    # (1) pool of all predictor names (sW);
                                    # (2) regression-specific predictor names; or
                                    # (3) list (of length(outvar)) of several predicator vectors, each such set is used as a separate regression for a node name in outvar
    reg_hazard = FALSE,            # If TRUE, the joint P(outvar|predvars) is factorized as \prod_{j}{P(outvar[j] | predvars)} for each j outvar (for fitting hazard)
    subset = NULL,                 # subset expression (later evaluated to logical vector in the envir of the data)
    ReplMisVal0 = TRUE,            # if TRUE all gvars$misval among predicators are replaced with with gvars$misXreplace (0)
    nbins = NULL,                  # actual nbins used, for cont. outvar, defined in ContinSummaryModel$new()
    bin_nms = NULL,                # column names for bin indicators
    bin_estimator = NULL,          # The algorithm to use for fitting the regressions
    parfit = logical(),            # TRUE for fitting binary regressions in parallel
    bin_bymass = logical(),        # for cont outvar, create bin cutoffs based on equal mass distribution?
    bin_bydhist = logical(),       # if TRUE, use dhist approach for bin definitions
    max_nperbin = integer(),       # maximum n observations allowed per binary bin
    pool = logical(),         # Pool binned cont outvar obs into long format (adding bin_ID as a covaraite)
    outvars_to_pool = character(), # Names of the binned continuous sVars, should match bin_nms
    intrvls.width = 1L,            # Named vector of bin-widths (bw_j : j=1,...,M) for each each bin in self$intrvls
                                   # When sA is not continuous, intrvls.width IS SET TO 1.
                                   # When sA is continuous, intrvls.width is SET TO self$intrvls.width INSIDE ContinSummaryModel$new() with names(intrvls.width) <- reg$bin_nms
                                   # CAN BE QUERIED BY BinOutModel$predictAeqa() as: intrvls.width[outvar]
    intrvls = NULL,                # Vector of numeric cutoffs defining the bins or a named list of numeric intervals (for length(self$outvar) > 1)
    levels = NULL,
    weights = NULL,
    # family = NULL,               # (NOT IMPLEMENTED) to run w/ other than "binomial" family
    # form = NULL,                 # (NOT IMPLEMENTED) reg formula, if provided run using the usual glm / speedglm functions
    # Adding ReplMisVal0 = TRUE below for case sA = (netA, sA[j]) with sA[j] continuous, was causing an error otherwise
    initialize = function(sep_predvars_sets = FALSE,
                          outvar.class = gvars$sVartypes$bin,
                          outvar, predvars, subset, intrvls = NULL,
                          ReplMisVal0 = TRUE,
                          bin_estimator = speedglmR6$new(),
                          parfit = FALSE,
                          nbins = NA_integer_,
                          bin_bymass = TRUE,
                          bin_bydhist = FALSE,
                          max_nperbin = 1000,
                          pool = FALSE,
                          weights = weights
                          ) {

      assert_that(length(outvar.class) == length(outvar))
      self$sep_predvars_sets <- sep_predvars_sets
      self$outvar.class <- outvar.class
      self$outvar <- outvar
      self$predvars <- predvars

      if (self$sep_predvars_sets) {
        assert_that(is.list(predvars))
        assert_that(is.list(outvar))
        assert_that(length(predvars) == length(outvar))
        if (length(names(predvars)) == 0) stop("when predvars is a list it must be named according to node names in outvar")
        assert_that(all(names(predvars) %in% names(outvar)))
      }

      self$ReplMisVal0 <- ReplMisVal0
      self$bin_estimator <- bin_estimator
      self$parfit <- parfit
      self$nbins <- nbins
      self$bin_bymass <- bin_bymass
      self$bin_bydhist <- bin_bydhist
      self$max_nperbin <- max_nperbin
      self$pool <- pool
      self$weights <- weights

      if (!missing(intrvls)) {
        if (!is.null(intrvls)) {
          assert_that(is.list(intrvls))
          assert_that(length(outvar) == length(intrvls))
          assert_that(all(names(intrvls) %in% outvar))
          self$intrvls <- intrvls
        }
      } else {
        self$intrvls <- NULL
      }

      self$levels <- NULL

      n_regs <- length(outvar)

      if (!missing(subset)) {
        self$subset <- subset
        if (length(subset) < n_regs) {
          self$subset <- rep_len(subset, n_regs)
        } else if (length(subset) > n_regs) {
          # ... TO FINISH ...
          # increase n_regs to all combinations of (n_regs x subset)
          # stop("subset cannot be longer than length(outvar)")
          if (!is.logical(subset)) stop("not implemented")
          # ... TO FINISH ...
        }
      } else {
        # self$subset <- rep_len(list(TRUE), n_regs)
        self$subset <- rep_len(list(NULL), n_regs)
      }
    },

    # take the clone of a parent RegressionClass (reg) for length(self$outvar) regressions
    # and set self to a single univariate k_i regression for outcome self$outvar[[k_i]]
    ChangeManyToOneRegresssion = function(k_i, reg) {
      self$resetS3class()

      assert_that(!missing(k_i))
      if (missing(reg)) stop("reg must be also specified when k_i is specified")
      assert_that(is.count(k_i))
      assert_that(k_i <= length(reg$outvar))

      n_regs <- length(reg$outvar)
      self$outvar.class <- reg$outvar.class[[k_i]] # Class of the outcome var: binary, categorical, continuous:
      self$outvar <- reg$outvar[[k_i]]             # An outcome variable that is being modeled:

      # modeling separate regressions (with different predictors for each outcome in outvar)
      if (reg$sep_predvars_sets) {
        self$sep_predvars_sets <- FALSE # (1) set the children objects to sep_predvars_sets <- FALSE
        if (!(names(reg$predvars)[k_i] %in% names(reg$outvar)[k_i])) stop("the names of list items in predvars must be in the same order as the node names in outvar")
        predvars_1 <- reg$predvars[[names(reg$outvar)[k_i]]]
        predvars_2 <- reg$predvars[[k_i]]
        if (!identical(predvars_1, predvars_2)) stop("fatal error: look up of reg$predvars[[...]] by outvar and by index k_i returning different results")
        self$predvars <- reg$predvars[[k_i]] # (2) extract specific predictors from the named list reg$predvars for each outvar

      # modeling bin hazard indicators, no need to condition on previous outcomes as they will all be degenerate
      } else if (self$reg_hazard) {
        self$predvars <- reg$predvars # Predictors

      # factorization of the joint prob P(A,B,C|D):=P(A|D)*P(B|A,D)*P(C|A,B,D)
      } else {
        self$predvars <- c(reg$outvar[-c(k_i:n_regs)], reg$predvars) # Predictors
      }

      # The subset becomes a list when a single RegressionClass specifies several regression models.
      # Each list item specifies the appropriate subset for regression idnex k_i.
      # On the other hand, if subset is a vector of variable names, all of those variables will be used for
      # choosing the subsets for all n_regs regressions.
      if (is.list(reg$subset)) {
        self$subset <- reg$subset[[k_i]]
      }
      if (is.list(reg$intrvls)) {
        outvar_idx <- which(names(reg$intrvls) %in% self$outvar)
        self$intrvls <- reg$intrvls[[outvar_idx]]
      }

      # Setting the self class for S3 dispatch on SummaryModel type
      if (reg$sep_predvars_sets) {
        self$S3class <- "generic" # Multivariate/Univariate regression at the top level, need to do another round of S3 dispatch on SummaryModel
      } else if (length(self$outvar)==1) {
        self$S3class <- self$outvar.class # Set class on outvar.class for S3 dispatch...
      } else if (length(self$outvar) > 1){
        stop("can't define a univariate regression for an outcome of length > 1")
      } else {
        stop("can't have an outcome with no class type")
      }

      return(invisible(self))
    },

    # take the clone of a parent RegressionClass for univariate (cont outvar) regression
    # and set self to length(regs_list) bin indicator outcome regressions
    ChangeOneToManyRegresssions = function(regs_list) {
      self$outvar.class <- regs_list$outvar.class # Vector of class(es) of outcome var(s): binary, categorical, continuous
      self$outvar <- regs_list$outvar # An outcome variable that is being modeled:
      self$predvars <- regs_list$predvars
      self$subset <- regs_list$subset
      return(invisible(self))
    },

    resetS3class = function() class(self) <- c("RegressionClass", "R6")

  ),

  active = list(
    # For S3 dispatch on newsummarymodel():
    S3class = function(newclass) {
      if (!missing(newclass)) {
        if (length(class(self)) > 2) stop("S3 dispatch class on RegressionClass has already been set")

        if (length(newclass) > 1) stop("cannot set S3 class on RegressionClass with more than one outvar variable")

        class(self) <- c(class(self), newclass)
      } else {
        return(class(self))
      }
    },

    get.reg = function() {
      list(outvar.class = self$outvar.class,
          outvar = self$outvar,
          predvars = self$predvars,
          subset = self$subset, 
          pool = self$pool,
          weights = self$weights)
    }
  )
)