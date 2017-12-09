#-----------------------------------------------------------------------------
# DataStore CLASS STRUCTURE:
#-----------------------------------------------------------------------------
# Contains Methods for:
  # *) detecting sVar types (detect.col.types);
  # *) normalizing continous sVar (normalize_sVar)
  # *) defining interval cuttoffs for continuous sVar (define.intervals)
  # *) turning continuous sVar into categorical (discretize.sVar)
  # *) creating binary indicator matrix for continous/categorical sVar (binirize.sVar, binirize.cat.sVar)
  # *) creating design matrix (Xmat) based on predvars and row subsets (evalsubst)

is.integerish <- function (x) is.integer(x) || (is.numeric(x) && all(x == as.integer(x)))

#' @importFrom stats as.formula glm na.exclude rbinom
NULL

## ---------------------------------------------------------------------
# DETECTING VECTOR TYPES
# sVartypes <- list(bin = "binary", cat = "categor", cont = "contin")
## ---------------------------------------------------------------------
detect.col.types <- function(sVar_mat, Yvars, max_n_cat){
  detect_vec_type <- function(vec) {
    vec_nomiss <- vec[!gvars$misfun(vec)]
    nvals <- length(unique(vec_nomiss))
    if (nvals <= 2L) {
      sVartypes$bin
    } else if ((nvals <= max_n_cat) && (is.integerish(vec_nomiss))) {
      sVartypes$cat
    } else {
      sVartypes$cont
    }
  }

  assert_that(is.integerish(max_n_cat) && (max_n_cat > 1))
  sVartypes <- gvars$sVartypes

  if (missing(Yvars)) Yvars <- colnames(sVar_mat)

  # for matrix:
  if (is.matrix(sVar_mat)) {
    return(as.list(apply(sVar_mat[, Yvars], 2, detect_vec_type)))
  # for data.table:
  } else if (is.data.table(sVar_mat)) {
    return(as.list(sVar_mat[, Yvars, with = FALSE][, lapply(.SD, detect_vec_type)]))
  } else {
    stop("unrecognized input data class: " %+% class(sVar_mat))
  }
}

## ---------------------------------------------------------------------
# Normalizing / Defining bin intervals / Converting contin. to ordinal / Converting ordinal to bin indicators
## ---------------------------------------------------------------------
normalize <- function(x) {
  if (abs(max(x) - min(x)) > gvars$tolerr) { # Normalize to 0-1 only when x is not constant
    return((x - min(x)) / (max(x) - min(x)))
  } else {  # What is the thing to do when x constant? Set to abs(x), abs(x)/x or 0???
    return(x)
  }
}
normalize_sVar <- function(sVar_vec) {
  nonmiss_idx <- !gvars$misfun(sVar_vec)
  if (sum(nonmiss_idx) > 0) {
    sVar_vec[nonmiss_idx] <- normalize(sVar_vec[nonmiss_idx])
  }
  sVar_vec
}
normalize_matsVar <- function(sVar_mat) apply(sVar_mat, 2, normalize_sVar)

# Define bin cutt-offs for continuous x:
define.intervals <- function(x, nbins, bin_bymass, bin_bydhist, max_nperbin) {
  x <- x[!gvars$misfun(x)]  # remove missing vals
  nvals <- length(unique(x))
  if (is.na(nbins)) {
    if (is.na(max_nperbin)) stop("must specify max_n_bin when setting nbins = NA")
    nbins <- max(5L, as.integer(length(x) / max_nperbin))
  }
  # if nbins is too high, for ordinal, set nbins to n unique obs and cancel quantile based interval defns
  if (nvals < nbins) {
    nbins <- nvals
    bin_bymass <- FALSE
  }
  if (abs(max(x) - min(x)) > gvars$tolerr) {  # when x is not constant
    if ((bin_bymass) & !is.na(max_nperbin)) {
      if ((length(x) / max_nperbin) > nbins) nbins <- as.integer(length(x) / max_nperbin)
    }
    intvec <- seq.int(from = min(x), to = max(x) + 1, length.out = (nbins + 1)) # interval type 1: bin x by equal length intervals of 0-1
    # intvec <- c(min(x), sort(runif(n = nbins, min = min(x), max = max(x))), max(x) + 1)
  } else {  # when x is constant, force the smallest possible interval to be at least [0,1]
    intvec <- seq.int(from = min(0L, min(x)), to = max(1L, max(x)), length.out = (nbins + 1))
  }
  if (bin_bymass) {
    intvec <- quantile(x = x, probs = normalize(intvec)) # interval type 2: bin x by mass (quantiles of 0-1 intvec as probs)
    intvec[1] <- intvec[1] - 0.01
    intvec[length(intvec)] <- intvec[length(intvec)] + 0.01
  } else if (bin_bydhist) {
    intvec <- dhist(x, plot = FALSE, nbins = nbins)$xbr
    intvec[1] <- intvec[1] - 0.01
    intvec[length(intvec)] <- intvec[length(intvec)] + 0.01
  }
  # adding -Inf & +Inf as leftmost & rightmost cutoff points to make sure all future data points end up in one of the intervals:
  # intvec <- c(-Inf, min(intvec)-0.01, intvec)
  # intvec <- c(min(intvec) - 0.1, intvec)
  intvec <- c(min(intvec)-1000L, intvec, max(intvec)+1000L)
  # intvec <- c(min(intvec) - 9999999L, min(intvec) - 0.1, intvec, max(intvec) + 0.1, max(intvec) + 9999999L)
  return(intvec) # return(list(intbylen = intvec, intbymass = intvecq))
}

# Turn any x into ordinal (1, 2, 3, ..., nbins) for a given interval cutoffs (length(intervals)=nbins+1)
make.ordinal <- function(x, intervals) findInterval(x = x, vec = intervals, rightmost.closed = TRUE)

# Make dummy indicators for ordinal x (sA[j])
# Approach: creates B_j that jumps to 1 only once and stays 1 (degenerate) excludes reference category (last)
make.bins_mtx_1 <- function(x.ordinal, nbins, bin.nms, levels = 1:nbins) {
  n <- length(x.ordinal)
  new.cats <- 1:nbins
  dummies_mat <- matrix(1L, nrow = n, ncol = length(new.cats))
  for(cat in new.cats[-length(new.cats)]) {
    subset_Bj0 <- x.ordinal > levels[cat]
    dummies_mat[subset_Bj0, cat] <- 0L
    subset_Bjmiss <- x.ordinal < levels[cat]
    dummies_mat[subset_Bjmiss, cat] <- gvars$misval
  }
  dummies_mat[, new.cats[length(new.cats)]] <- gvars$misval
  colnames(dummies_mat) <- bin.nms
  dummies_mat
}


## ---------------------------------------------------------------------
#' R6 class for storing and managing the combined summary measures \code{sW} & \code{sA} from DatNet classes.
#'
#' This class inherits from \code{DatNet} and extends its methods to handle a single matrix dataset of
#'  all summary measures \code{(sA,sW)}
#'  The class \code{DataStore} is the only way to access data in the entire package.
#'  Contains methods for combining, subsetting, discretizing & binirizing summary measures \code{(sW,sA)}.
#'  For continous sVar this class provides methods for detecting / setting bin intervals,
#'  normalization, disretization and construction of bin indicators.
#'  The pointers to this class get passed on to \code{SummariesModel} functions: \code{$fit()},
#'  \code{$predict()} and \code{$predictAeqa()}.
#'
#' @docType class
#' @format An \code{\link{R6Class}} generator object
#' @keywords R6 class
#' @details
#' \itemize{
#' \item{\code{datnetW}} - .
#' \item{\code{datnetA}} - .
#' \item{\code{active.bin.sVar}} - Currently discretized continous \code{sVar} column in data matrix \code{mat.sVar}.
#' \item{\code{mat.bin.sVar}} - Matrix of the binary indicators for discretized continuous covariate \code{active.bin.sVar}.
#' \item{\code{ord.sVar}} - Ordinal (categorical) transformation of a continous covariate \code{sVar}.
#' \item{\code{YnodeVals}} - .
#' \item{\code{det.Y}} - .
#' \item{\code{p}} - .
#' }
#' @section Methods:
#' \describe{
#'   \item{\code{new(datnetW, datnetA, YnodeVals, det.Y, ...)}}{...}
#'   \item{\code{addYnode(YnodeVals, det.Y)}}{...}
#'   \item{\code{evalsubst(subsetexpr, subsetvars)}}{...}
#'   \item{\code{get.dat.sWsA(rowsubset = TRUE, covars)}}{...}
#'   \item{\code{get.outvar(rowsubset = TRUE, var)}}{...}
#'   \item{\code{copy.sVar.types()}}{...}
#'   \item{\code{bin.nms.sVar(name.sVar, nbins)}}{...}
#'   \item{\code{pooled.bin.nm.sVar(name.sVar)}}{...}
#'   \item{\code{detect.sVar.intrvls(name.sVar, nbins, bin_bymass, bin_bydhist, max_nperbin)}}{...}
#'   \item{\code{detect.cat.sVar.levels(name.sVar)}}{...}
#'   \item{\code{discretize.sVar(name.sVar, intervals)}}{...}
#'   \item{\code{binirize.sVar(name.sVar, intervals, nbins, bin.nms)}}{...}
#'   \item{\code{binirize.cat.sVar(name.sVar, levels)}}{...}
#'   \item{\code{get.sVar.bw(name.sVar, intervals)}}{...}
#'   \item{\code{get.sVar.bwdiff(name.sVar, intervals)}}{...}
#'   \item{\code{make.dat.sWsA(p = 1, f.g_fun = NULL, sA.object = NULL)}}{...}
#' }
#' @section Active Bindings:
#' \describe{
#'    \item{\code{dat.sWsA}}{...}
#'    \item{\code{dat.bin.sVar}}{...}
#'    \item{\code{emptydat.bin.sVar}}{...}
#'    \item{\code{names.sWsA}}{...}
#'    \item{\code{nobs}}{...}
#'    \item{\code{noNA.Ynodevals}}{...}
#'    \item{\code{nodes}}{...}
#' }
#' @importFrom assertthat assert_that is.count is.flag
#' @export
DataStore <- R6Class(classname = "DataStore",
  portable = TRUE,
  class = TRUE,
  public = list(
    max_n_cat = NULL,           # Max number of unique levels for cat outcome. If the outcome has more levels it is automatically considered continuous.
    mat.sVar = NULL,           # Matrix storing all evaluated sVars, with named columns
    type.sVar = NULL,          # named list with sVar types: list(names.sVar[i] = "binary"/"categor"/"contin"), can be overridden
    norm.c.sVars = FALSE,      # flag = TRUE if want to normalize continous covariates
    active.bin.sVar = NULL,    # name of active binarized cont sVar, changes as fit/predict is called (bin indicators are temp. stored in mat.bin.sVar)
    mat.bin.sVar = NULL,       # temp storage mat for bin indicators on currently binarized continous sVar (from self$active.bin.sVar)
    ord.sVar = NULL,           # Ordinal (cat) transform for continous sVar
    weights = NULL,

    initialize = function(input_data, X, Y, auto_typing = TRUE, max_n_cat = 20, weights = NULL, ...) {
      if (!is.data.table(input_data)) data.table::setDT(input_data)
      self$dat.sVar <- input_data
      self$nodes <- list(X = X, Y = Y)

      self$max_n_cat <- max_n_cat
      self$weights <- weights

      ## We need to evaluate the types of OUTCOME VARIABLES ONLY (the types of predictors are irrelevant)
      if (auto_typing) self$def.types.sVar(Y = Y) # Define the type of each Y: bin, cat or cont

      # normalize continuous and non-missing sVars, overwrite their columns in mat.sVar with normalized [0,1] vals
      if (self$norm.c.sVars) self$norm_c_sVars()

      return(invisible(self))
    },

    # Define the type (class) of each summary measure: bin, cat or cont
    # type.sVar acts as a flag: only detect types when !is.null(type.sVar)
    # otherwise can pass type.sVar = list(sVar = NA, ...) or a value type.sVar = NA/gvars$sVartypes$bin/etc
    def.types.sVar = function(type.sVar = NULL, Y) {
      # Detect the type of each sVar[i]: gvars$sVartypes$bin,  gvars$sVartypes$cat, gvars$sVartypes$cont
      if (is.null(type.sVar)) {
        self$type.sVar <- detect.col.types(self$dat.sVar, Y, self$max_n_cat)
      } else {
        n.sVar <- length(self$names.sVar)
        len <- length(type.sVar)
        assert_that((len == n.sVar) || (len == 1L))
        if (len == n.sVar) {
          assert_that(is.list(type.sVar))
          assert_that(all(names(type.sVar) %in% self$names.sVar))
        } else {
          assert_that(is.string(type.sVar))
          type.sVar <- as.list(rep(type.sVar, n.sVar))
          names(type.sVar) <- self$names.sVar
        }
        self$type.sVar <- type.sVar
      }
      invisible(self)
    },

    # Normalize continuous sVars # This could be memory-costly
    norm_c_sVars = function() {
      names.c.sVar <- self$names.c.sVar
      if (length(names.c.sVar) == 0L) return(invisible(self))

      if (self$norm.c.sVars && (length(names.c.sVar) > 0)) {
        for (name.c.sVar in names.c.sVar) {
          self$mat.sVar[, name.c.sVar] <- normalize_sVar(self$mat.sVar[, name.c.sVar])
        }
      }
      invisible(self)
    },

    # Eval the expression (in the environment of the data.frame "data" + global constants "gvars"):
    evalsubst = function(subsetexpr, subsetvars) {
      if (missing(subsetexpr)) {
        assert_that(!missing(subsetvars))
        res <- rep.int(TRUE, self$nobs)
        for (subsetvar in subsetvars) {
          # *) find the var of interest (in self$dat.sWsA or self$dat.bin.sVar), give error if not found
          sVar.vec <- self$get.outvar(var = subsetvar)
          assert_that(!is.null(sVar.vec))
          # *) reconstruct correct expression that tests for missing values
          res <- res & (!gvars$misfun(sVar.vec))
        }
        return(which(res))
      # ******************************************************
      # NOTE: Below is currently not being used, all subsetting now is done with subsetvars above, for speed & memory efficiency
      # ******************************************************
      } else {
        if (is.logical(subsetexpr)) {
          return(which(subsetexpr))
        } else {
          # ******************************************************
          # THIS WAS A BOTTLENECK: for 500K w/ 1000 bins: 4-5sec
          # REPLACING WITH env that is made of data.frames instead of matrices
          # ******************************************************
          # eval.env <- c(data.frame(self$dat.sWsA), data.frame(self$dat.bin.sVar), as.list(gvars))
          # res <- try(eval(subsetexpr, envir = eval.env, enclos = baseenv())) # to evaluate vars not found in data in baseenv()
          stop("disabled for memory/speed efficiency")
          return(res)
        }
      }
    },

    # return a covar matrix which will be used as a design matrix for BinOutModelClass
    get.dat.sWsA = function(rowsubset, covars) {
      if (!missing(covars)) {

        if (length(unique(colnames(self$dat.sWsA))) < length(colnames(self$dat.sWsA))) {
          warning("repeating column names in the final data set; please check for duplicate summary measure / node names")
        }

        # columns to select from main design matrix (in the same order as listed in covars):
        sel.sWsA <- intersect(covars, colnames(self$dat.sWsA))

        if (missing(rowsubset)) {
          rowsubset <- 1:nrow(self$dat.sWsA)
        } else if (is.logical(rowsubset)) {
          rowsubset <- which(rowsubset)
        } else if (!is.integer(rowsubset)) {
          stop("rowsubset must be logical or integer subset of rows")
        }

        if (is.matrix(self$dat.sWsA)) {
          dfsel <- self$dat.sWsA[rowsubset, sel.sWsA, drop = FALSE] # data stored as matrix
        } else if (is.data.table(self$dat.sWsA)) {
          dfsel <- self$dat.sWsA[rowsubset, sel.sWsA, with = FALSE] # data stored as data.table
        } else {
          stop("self$dat.sWsA is of unrecognized class: " %+% class(self$dat.sWsA))
        }

        # columns to select from binned continuous/cat var matrix (if it was previously constructed):
        if (!is.null(self$dat.bin.sVar)) {
          sel.binsA <- intersect(covars, colnames(self$dat.bin.sVar))
        } else {
          sel.binsA <- NULL
        }
        if (length(sel.binsA)>0) {
          if (ncol(dfsel)==0) {
            dfsel <- as.data.table(self$dat.bin.sVar[rowsubset, sel.binsA, drop = FALSE])
          } else {
            dfsel <- cbind(dfsel, self$dat.bin.sVar[rowsubset, sel.binsA, drop = FALSE])
          }
        }

        found_vars <- covars %in% colnames(dfsel)
        if (!all(found_vars)) stop("some covariates can't be found (perhaps not declared as summary measures (def_sW(...) or def_sW(...))): "%+%
                                    paste(covars[!found_vars], collapse=","))
        return(dfsel)
      } else {
        return(self$dat.sWsA[rowsubset, , drop = FALSE])
      }
    },

    get.outvar = function(rowsubset = TRUE, var) {
      # if (length(self$nodes) < 1) stop("DataStore$nodes list is empty!")

      if (var %in% self$names.sWsA) {
        out <- self$dat.sWsA[rowsubset, var, with = FALSE]
      } else if (var %in% colnames(self$dat.bin.sVar)) {
        out <- self$dat.bin.sVar[rowsubset, var]
      } else {
        stop("requested variable " %+% var %+% " does not exist in DataStore!")
      }

      if ((is.list(out) || is.data.table(out)) && (length(out)>1)) {
        stop("selecting regression outcome covariate resulted in more than one column: " %+% var)
      } else if (is.list(out) || is.data.table(out)) {
        return(out[[1]])
      } else {
        return(out)
      }

    },

    get.wts = function(rowsubset = TRUE) {
      if (!is.null(self$weights)) {
        return(self$get.outvar(rowsubset, self$weights))
      } else {
        return(NULL)
      }
    },

    # Need to find a way to over-ride nbins for categorical vars (allowing it to be set to more than gvars$max_n_cat)!
    # Return names of bin indicators for sVar:
    bin.nms.sVar = function(name.sVar, nbins) { name.sVar%+%"_"%+%"B."%+%(1L:nbins) },
    pooled.bin.nm.sVar = function(name.sVar) { name.sVar %+% "_allB.j" },
    detect.sVar.intrvls = function(name.sVar, nbins, bin_bymass, bin_bydhist, max_nperbin) {
      tol.int <- 0.001
      int <- define.intervals(x = self$get.sVar(name.sVar), nbins = nbins, bin_bymass = bin_bymass, bin_bydhist = bin_bydhist, max_nperbin = max_nperbin)
      diffvec <- diff(int)
      if (sum(abs(diffvec) < tol.int) > 0) {
      # if (length(unique(int)) < length(int)) {
        if (gvars$verbose) {
          # message("No. of categories for " %+% name.sVar %+% " was collapsed from " %+%
          #         (length(int)-1) %+% " to " %+% (length(unique(int))-1) %+% " due to too few obs.")
          message("No. of categories for " %+% name.sVar %+% " was collapsed from " %+%
                  (length(int)-1) %+% " to " %+% (length(int[diffvec >= tol.int])-1) %+% " due to too few obs.")
          print("old intervals: "); print(as.numeric(int))
        }
        # Just taking unique interval values is insufficient
        # Instead need to drop all intervals that are "too close" to each other based on some tol value
        # remove all intervals (a,b) where |b-a| < tol.int, but always keep the very first interval (int[1])
        int <- c(int[1], int[2:length(int)][abs(diffvec) >= tol.int])
        # int <- unique(int)
        if (gvars$verbose) {
          print("new intervals: "); print(as.numeric(int))
        }
      }
      return(int)
    },

    detect.cat.sVar.levels = function(name.sVar) {
      levels <- sort(unique(self$get.sVar(name.sVar)))
      return(levels)
    },

    # create a vector of ordinal (categorical) vars out of cont. sVar vector:
    discretize.sVar = function(name.sVar, intervals) {
      self$ord.sVar <- make.ordinal(x = self$get.sVar(name.sVar), intervals = intervals)
      invisible(self$ord.sVar)
    },

    # return matrix of bin indicators for continuous sVar:
    # change name to:
    # binirize.cont.sVar = function(name.sVar, intervals, nbins, bin.nms) {
    binirize.sVar = function(name.sVar, intervals, nbins, bin.nms) {
      self$active.bin.sVar <- name.sVar
      self$dat.bin.sVar <- make.bins_mtx_1(x.ordinal = self$discretize.sVar(name.sVar, intervals), nbins = nbins, bin.nms = bin.nms)
      invisible(self$dat.bin.sVar)
    },

    # return matrix of bin indicators for ordinal sVar:
    binirize.cat.sVar = function(name.sVar, levels) {
      nbins <- length(levels)
      bin.nms <- self$bin.nms.sVar(name.sVar, nbins)
      self$active.bin.sVar <- name.sVar
      self$dat.bin.sVar <- make.bins_mtx_1(x.ordinal = self$get.sVar(name.sVar), nbins = nbins, bin.nms = bin.nms, levels = levels)
      invisible(self$dat.bin.sVar)
    },

    # return the bin widths vector for the discretized continuous sVar (self$ord.sVar):
    get.sVar.bw = function(name.sVar, intervals) {
      if (!(self$active.bin.sVar %in% name.sVar)) stop("current discretized sVar name doesn't match name.sVar in get.sVar.bin.widths()")
      if (is.null(self$ord.sVar)) stop("sVar hasn't been discretized yet")
      intrvls.width <- diff(intervals)
      intrvls.width[intrvls.width <= gvars$tolerr] <- 1
      ord.sVar_bw <- intrvls.width[self$ord.sVar]
      return(ord.sVar_bw)
    },

   # return the bin widths vector for the discretized continuous sVar (self$ord.sVar):
    get.sVar.bwdiff = function(name.sVar, intervals) {
      if (!(self$active.bin.sVar %in% name.sVar)) stop("current discretized sVar name doesn't match name.sVar in get.sVar.bin.widths()")
      if (is.null(self$ord.sVar)) stop("sVar hasn't been discretized yet")
      ord.sVar_leftint <- intervals[self$ord.sVar]
      diff_bw <- self$get.sVar(name.sVar) - ord.sVar_leftint
      return(diff_bw)
    },

    fixmiss_sVar_mat = function() {
      self$dat.sVar[gvars$misfun(self$dat.sVar)] <- gvars$misXreplace
      invisible(self)
    },

    fixmiss_sVar_DT = function() {
      # see http://stackoverflow.com/questions/7235657/fastest-way-to-replace-nas-in-a-large-data-table
      dat.sVar <- self$dat.sVar
      for (j in names(dat.sVar))
        set(dat.sVar, which(gvars$misfun(dat.sVar[[j]])), j , gvars$misXreplace)
      invisible(self)
    },

    fixmiss_sVar = function() {
      if (is.matrix(self$dat.sVar)) {
        self$fixmiss_sVar_mat()
      } else if (is.data.table(self$dat.sVar)) {
        self$fixmiss_sVar_DT()
      } else {
        stop("self$dat.sVar is of unrecognized class")
      }
    },

    # --------------------------------------------------
    # Methods for directly handling one continous/categorical sVar in self$mat.sVar;
    # No checking of incorrect input is performed, use at your own risk!
    # --------------------------------------------------
    norm.sVar = function(name.sVar) { normalize_sVar(self$dat.sVar[, name.sVar]) },  # return normalized 0-1 sVar
    set.sVar = function(name.sVar, new.sVar) { self$mat.sVar[, .(name.sVar) := eval(new.sVar)]},
    get.sVar = function(name.sVar) {
      x <- self$dat.sVar[, name.sVar, with=FALSE]
      if (is.list(x) || is.data.table(x) || is.data.frame(x)) x <- x[[1]]
      return(x)
    },
    set.sVar.type = function(name.sVar, new.type) { self$type.sVar[[name.sVar]] <- new.type },
    get.sVar.type = function(name.sVar) { if (missing(name.sVar)) { self$type.sVar } else { self$type.sVar[[name.sVar]] } }
  ),

  active = list(
    dat.bin.sVar = function(dat.bin.sVar) {
      if (missing(dat.bin.sVar)) {
        return(self$mat.bin.sVar)
      } else {
        assert_that(is.matrix(dat.bin.sVar))
        self$mat.bin.sVar <- dat.bin.sVar
      }
    },

    nobs = function() { nrow(self$dat.sVar) },
    ncols.sVar = function() { length(self$names.sVar) },

    names.sWsA = function() { self$names.sVar },
    names.sVar = function() { colnames(self$dat.sVar) },
    names.c.sVar = function() { names(self$type.sVar[self$type.sVar %in% gvars$sVartypes$cont]) }, # names of cont sVars

    dat.sWsA = function() { self$mat.sVar }, # NO LONGER NEEDED, using only self$dat.sVar, KEPT FOR COMPATIBILITY
    dat.sVar = function(dat.sVar) {
      if (missing(dat.sVar)) {
        return(self$mat.sVar)
      } else {
        assert_that(is.matrix(dat.sVar) | is.data.table(dat.sVar))
        self$mat.sVar <- dat.sVar
      }
    },

    emptydat.sVar = function() { self$mat.sVar <- NULL },         # wipe out mat.sVar

    # wipe out binirized mat.sVar:
    emptydat.bin.sVar = function() {
      self$mat.bin.sVar <- NULL
      self$active.bin.sVar <- NULL
    },

    nodes = function(nodes) {
      if (missing(nodes)) {
        return(private$.nodes)
      } else {
        assert_that(is.list(nodes))
        if (length(self$nodes)>0) message("warning: overwriting non-empty self$nodes in DataStore")
        private$.nodes <- nodes
      }
    }
  ),

  private = list(
    .nodes = list()           # names of the nodes in the data (Anode, Ynode, nFnode, etc..)
  )
)