#################################################################################
########################### NETWORK TMLE R PACKAGE ##############################
# authors: Oleg Sofrygin <sofrygin@berkeley.edu> and Mark van der Laan <laan@berkeley.edu>
#################################################################################

# @title tmlenet-package
# @docType package
# @author Oleg Sofrygin, Mark J. van der Laan

#' @useDynLib condensier
#' @import R6
#' @importFrom Rcpp sourceCpp
#' @importFrom grDevices nclass.FD nclass.Sturges nclass.scott
#' @importFrom graphics axis barplot hist par text
#' @importFrom methods is
#' @importFrom stats approx quasibinomial binomial coef glm.control glm.fit plogis predict qlogis qnorm quantile rnorm terms var predict glm.control
#' @importFrom utils data head str

#-----------------------------------------------------------------------------
# Class Membership Tests
#-----------------------------------------------------------------------------
is.DataStore <- function(DataStore) "DataStore"%in%class(DataStore)

# #-----------------------------------------------------------------------------
# # ALL NETWORK VARIABLE NAMES MUST BE CONSTRUCTED BY CALLING THIS FUNCTION.
# # In the future might return the network variable (column vector) itself.
# # Helper function that for given variable name (varnm) and friend index (fidx)
# # returns the characeter name of that network variable varnm[fidx],
# # for fidx = 0 (var itself), ..., kmax. fidx can be a vector, in which case a
# # character vector of network names is returned. If varnm is also a vector, a
# # character vector for all possible combinations of (varnm x fidx) is returned.
# # OUTPUT format: Varnm_net.j:
# #-----------------------------------------------------------------------------
# netvar <- function(varnm, fidx) {
#   cstr <- function(varnm, fidx) {
#     slen <- length(fidx)
#     rstr <- vector(mode = "character", length = slen)
#     netidxstr <- ! (fidx %in% 0L)
#     rstr[netidxstr] <- stringr::str_c('_netF', fidx[netidxstr])  # vs. 1
#     # rstr[netidxstr] <- str_c('.net.', fidx[netidxstr])  # vs. 2
#     return(stringr::str_c(varnm, rstr))
#   }
#   if (length(varnm) > 1) {
#     return(unlist(lapply(varnm, cstr, fidx)))
#   } else {
#     return(cstr(varnm, fidx))
#   }
# }
# # Examples:
# # netvar("A", (0:5))
# # netvar("A", c(0:5, 0, 3))
# # netvar(c("A", "W"), c(0:5, 0, 3))
# # netvar(c("A", "W"), c(0:5, 0, 3))

#-----------------------------------------------------------------------------
# General utilities / Global Vars
#-----------------------------------------------------------------------------
`%+%` <- function(a, b) paste0(a, b)

checkpkgs <- function(pkgs) {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(pkg %+% " package needed for this function to work. Please install it.", call. = FALSE)
    }
  }
}

# Bound g(A|W) probability within supplied bounds
bound <- function(x, bounds){
  x[x<min(bounds)] <- min(bounds)
  x[x>max(bounds)] <- max(bounds)
  return(x)
}

#if warning is in ignoreWarningList, ignore it; otherwise post it as usual
SuppressGivenWarnings <- function(expr, warningsToIgnore) {
  h <- function (w) {
    if (w$message %in% warningsToIgnore) invokeRestart( "muffleWarning" )
  }
  withCallingHandlers(expr, warning = h )
}

GetWarningsToSuppress <- function(update.step=FALSE) {
  warnings.to.suppress <- c("glm.fit: fitted probabilities numerically 0 or 1 occurred",
                            "prediction from a rank-deficient fit may be misleading")
  if (update.step) {
    warnings.to.suppress <- c(warnings.to.suppress, "glm.fit: algorithm did not converge")
  }
  return(warnings.to.suppress)
}

# returns NULL if no factors exist, otherwise return the name of the factor variable(s)
CheckExistFactors <- function(data) {
  testvec <- unlist(lapply(data, is.factor))
  if (any(testvec)) {
    return(names(data)[which(testvec)])
  } else {
    return(NULL)
  }
}

# throw exception if 1) varname doesn't exist; 2) more than one varname is matched
CheckVarNameExists <- function(data, varname) {
  idvar <- names(data) %in% varname
  if (sum(idvar) < 1) stop("variable name " %+% varname %+% " not found in data input")
  if (sum(idvar) > 1) stop("more than one column in the input data has been matched to name "
                            %+% varname %+% ". Consider renaming some of the columns: " %+%
                            paste0(names(data)[idvar], collapse=","))
  return(invisible(NULL))
}


# Return the left hand side variable of formula f as a character
LhsVars <- function(f) {
  f <- as.formula(f)
  return(as.character(f[[2]]))
}
# Return the right hand side variables of formula f as a character vector
RhsVars <- function(f) {
  f <- as.formula(f)
  return(all.vars(f[[3]]))
}

#---------------------------------------------------------------------------------
# Estimate h_bar under g_0 and g* given observed data and vector of c^Y's data is an DataStore object
#---------------------------------------------------------------------------------
# get_all_ests <- function(estnames, DatNet.ObsP0, est_params_list) {
#   #---------------------------------------------------------------------------------
#   # unified estimator naming used throughout the package:
#   # c("TMLE", "h_IPTW", "MLE")
#   #---------------------------------------------------------------------------------
#   # DatNet.ObsP0$det.Y             # TRUE/FALSE for deterministic Y's
#   # DatNet.ObsP0$noNA.Ynodevals    # actual observed Y's
#   # m.Q.init$getoutvarnm           # reg outvar name (Ynode)
#   # DatNet.ObsP0$YnodeVals         # visible Y's with NA for det.Y
#   # m.Q.init$getoutvarval          # Yvals used in prediction (with det.Y obs set to NA)
#   # m.Q.init$getsubset             # valid subset (!det.Y)
#   # m.Q.init$reg                   # regression class (Qreg)

#   # making sure the current dataset consists of the OBSERVED Anodes and sA (summaries)
#   if (!DatNet.ObsP0$Odata$curr_data_A_g0) {
#     if (is.null(DatNet.ObsP0$Odata$A_g0_DT))
#       stop("Can't recover the initial observed A (exposures), as those were over-written and not backed-up")
#     DatNet.ObsP0$Odata$swapAnodes()
#     if (!DatNet.ObsP0$Odata$restored_sA_Vars) DatNet.ObsP0$datnetA$make.sVar(sVar.object = est_params_list$sA)
#   }

#   nodes <- DatNet.ObsP0$nodes
#   Y <- DatNet.ObsP0$noNA.Ynodevals # actual observed Y`s
#   determ.Q <- DatNet.ObsP0$det.Y
#   m.Q.init <- est_params_list$m.Q.init

#   QY.init <- DatNet.ObsP0$noNA.Ynodevals # getting all node vals, inc. deterministic
#   QY.init[!DatNet.ObsP0$det.Y] <- m.Q.init$predict(newdata = DatNet.ObsP0)$getprobA1[!DatNet.ObsP0$det.Y] # getting predictions P(Y=1) for non-DET Y
#   off <- qlogis(QY.init)  # offset log(x/[1-x])

#   #************************************************
#   # h^*/h_N clever covariate:
#   #************************************************
#   fit.hbars_t <- system.time(fit.hbars.res <- fit.hbars(DatNet.ObsP0 = DatNet.ObsP0, est_params_list = est_params_list)) # fit the clever covariate
#   DatNet.gstar <- fit.hbars.res$DatNet.gstar
#   m.h.fit <- fit.hbars.res$m.h.fit
#   h_wts <- fit.hbars.res$h_gstar_h_gN

#   #************************************************
#   # IPTW_h estimator:
#   #************************************************
#   h_IPTW <- Y
#   h_IPTW[!determ.Q] <- Y[!determ.Q] * h_wts[!determ.Q]
#   h_IPTW <- mean(h_IPTW)

#   # message("h_IPTW: "); message(h_IPTW)
#   # if (h_IPTW > 0.5) browser()

#   #************************************************
#   # Get a TMLE update:
#   #************************************************
#   # epsilon (intercept or coefficient for h): tmle.obj$m.Q.star
#   # boot_idx <- seq.int(DatNet.ObsP0$nobs)
#   tmle.obj <- tmle.update(estnames = estnames, Y = Y, off = off, h_wts = h_wts, determ.Q = determ.Q)
#   # print("tmle.obj$m.Q.star.coef"); print(tmle.obj$m.Q.star.coef)

#   #************************************************
#   # Run Monte-Carlo (MC) evaluation for all plug-in estimators (TMLE & Gcomp), under stochastic intervention g^*:
# 	#************************************************
#   MC_fit_params <- append(est_params_list, list(m.Q.star = tmle.obj$m.Q.star.coef))

#   # make sure the current dataset constits of INTERVENED Anodes and sA (under g.star):
#   if (DatNet.ObsP0$datnetW$Odata$curr_data_A_g0) {
#     DatNet.gstar$datnetA$Odata$swapAnodes()
#     if (!DatNet.gstar$datnetA$Odata$restored_sA_Vars)
#       DatNet.gstar$datnetA$make.sVar(sVar.object = est_params_list$new.sA) # just recreate the summaries sA based on current Anode values in the data
#   }

#   # generate new A's under f.gstar, then re-create new summaries sA:
#   # DatNet.gstar$make.dat.sWsA(p = 1, f.g_fun = est_params_list$f.gstar, sA.object = sA, DatNet.ObsP0 = DatNet.ObsP0)
#   # DatNet.gstar$datnetA$Odata$OdataDT
#   # DatNet.gstar$datnetA$Odata$A_g0_DT
#   # DatNet.gstar$datnetA$Odata$sA_g0_DT
#   # both should be the same, as both are pointing to the same data.table:
#   # head(DatNet.ObsP0$dat.sVar)
#   # head(DatNet.gstar$dat.sVar)

#   MC.Eval.psi <- get.MCS_ests(DatNet.ObsP0 = DatNet.ObsP0,
#                               DatNet.gstar = DatNet.gstar,
#                               MC_fit_params = MC_fit_params,
#                               m.h.fit = m.h.fit)

#   # -----------------------------------------------------------------
#   # Mean estimates
#   # -----------------------------------------------------------------
#   MCS_res <- MC.Eval.psi$psi_est_mean
#   ests <- c(TMLE = MCS_res[["TMLE"]],
#             h_IPTW = h_IPTW, # IPTW estimator based on h - clever covariate:
#             MLE = MCS_res[["MLE"]])

#   ests_mat <- matrix(0L, nrow = length(ests), ncol = 1)
#   ests_mat[, 1] <- ests
#   rownames(ests_mat) <- names(ests); colnames(ests_mat) <- "estimate"

#   wts_mat <- matrix(0L, nrow = DatNet.ObsP0$nobs, ncol = 1)
#   colnames(wts_mat) <- c("h_wts")
#   wts_mat[, "h_wts"] <- h_wts

#   # Components of asymptotic variance for tmle_net (E_g^*[\bar{Q}_0(A^s,W^s|W^s)]-\psi_0):
#   # SUBSTRACT overall estimate of psi_0 from fW_i i-specific components
#   fWi_mat <- matrix(0L, nrow = DatNet.ObsP0$nobs, ncol = 1)
#   colnames(fWi_mat) <- c("fWi_Qinit")
#   fWi_mat[,"fWi_Qinit"] <- MCS_res[agrep("fWi_init_", names(MCS_res), max.distance=list(insertions=0, substitutions=0))]

#   QY_mat <- matrix(0L, nrow = DatNet.ObsP0$nobs, ncol = 2)
#   colnames(QY_mat) <- c("QY.init", "QY.star")
#   QY_mat[,] <- cbind(QY.init, tmle.obj$QY.star)


#   if (gvars$verbose)  {
#     print("time spent fitting new fit.hbars.res:"); print(fit.hbars_t)
#     if ("TMLE_A" %in% estnames) {
#       parsubmodel_fits <- rbind(tmle.obj$m.Q.star.coef)
#       rownames(parsubmodel_fits) <- c("epsilon (clever covariate coefficient)")
#     } else if ("TMLE_B" %in% estnames) {
#       parsubmodel_fits <- rbind(tmle.obj$m.Q.star.coef)
#       rownames(parsubmodel_fits) <- c("alpha (intercept)")
#     }
#     print("new parsubmodel_fits: "); print(parsubmodel_fits)
#     print(
#           c(
#           fWi_init = mean(fWi_mat[,"fWi_Qinit"] - ests["TMLE"])
#           ));
#     print("new MC.ests mat: "); print(ests_mat)
#   }


#   # -----------------------------------------------------------------
#   # instance of an R6 object mcEvalPsi:
#   # -----------------------------------------------------------------
#   psi.evaluator <- MC.Eval.psi$psi.evaluator

#   return(list( ests_mat = ests_mat,
#                wts_mat = wts_mat,
#                fWi_mat = fWi_mat,
#                QY_mat = QY_mat,
#                # var_tmleB_boot = var_tmleB_boot,
#                psi.evaluator = psi.evaluator, # for par. bootstrap
#                m.Q.init = m.Q.init,           # for par. bootstrap
#                m.h.fit = m.h.fit,             # for par. bootstrap
#                DatNet.gstar = DatNet.gstar,   # for par. bootstrap
#                sW = est_params_list$sW,       # for par. bootstrap
#                sA = est_params_list$sA,       # for par. bootstrap
#                f.gstar = est_params_list$f.gstar, # for par. bootstrap
#                new.sA = est_params_list$new.sA, # for par. bootstrap
#                MC_fit_params = MC_fit_params,   # for Monte-Carlo eval of the iid W EIC
#                h_g0_SummariesModel = m.h.fit$summeas.g0,
#                h_gstar_SummariesModel = m.h.fit$summeas.gstar
#               ))
# }


#---------------------------------------------------------------------------------
# (NEW INTERFACE FOR SPECIFYING regressions for hform.g0, hform.gstar & Qform)
#---------------------------------------------------------------------------------
get_vars_fromlist <- function(varname, sVar.map) {
  if (varname %in% names(sVar.map)) {
    as.vector(sVar.map[[varname]])
  } else {
    varname
  }
}

# Parse the formulas for summary measure names and create a map to actual covariate names in sA & sW
process_regform <- function(regform, sW.map = NULL, sA.map = NULL) {
  # Getting predictors (sW names):
  regformterms <- terms(regform)
  sW.names <- attributes(regformterms)$term.labels
  sW.names.alt <- colnames(attributes(regformterms)$factors)
  assert_that(all(sW.names == sW.names.alt))

  # Getting outcomes (sA names):
  out.var <- rownames(attributes(regformterms)$factors)[1] # character string
  out.vars.form <- as.formula(". ~ " %+% out.var)
  out.vars.terms <- terms(out.vars.form)
  sA.names <- attributes(out.vars.terms)$term.labels

  outvars <- unlist(lapply(sA.names, get_vars_fromlist, sA.map))
  predvars <- unlist(lapply(sW.names, get_vars_fromlist, sW.map))
  return(list(outvars = outvars, predvars = predvars))
}

# When several reg forms are specified (multivariate Anodes), process outvars into one vector and process predvars in a named list of vectors
process_regforms <- function(regforms, sW.map = NULL, sA.map = NULL) {
  if (length(regforms)==0L) {
    default.reg <- list(outvars =  list(as.vector(unlist(sA.map))), predvars = list(as.vector(unlist(sW.map))))
    names(default.reg[["outvars"]]) <- names(default.reg[["predvars"]]) <- paste0(default.reg$outvars[[1]], collapse="+")
    return(default.reg)
  } else {
    outvars <- vector(mode="list", length=length(regforms))
    predvars <- vector(mode="list", length=length(regforms))
    for (idx in seq_along(regforms)) {
      res <- process_regform(as.formula(regforms[[idx]]), sW.map = sW.map, sA.map = sA.map)
      outvars[[idx]] <- res$outvars
      predvars[[idx]] <- res$predvars
      names(outvars)[idx] <- names(predvars)[idx] <- paste0(outvars[[idx]], collapse="+")
    }
    return(list(outvars = outvars, predvars = predvars))
  }
}

# Evaluate Summary Measures sA and sW
#
# Take input data, create a network matrix (when input network matrix not provided) and evaluate the summary measures
#  previously defined with functions \code{\link{def_sW}} and \code{\link{def_sA}}.
#  This function is called internally by \code{tmlenet} for the evaluation of the summary measures.
#  The R6 class object named \code{DatNet.ObsP0} that is returned by this function can be supplied as an input to the
#  \code{tmlenet} function.
#  When \code{DatNet.ObsP0} is used as an input to \code{tmlenet}, the rest of the input arguments already provided to
#  this function can be omitted from the \code{tmlenet} function call.
# @param data Same as \code{\link{tmlenet}} input argument.
# @param Kmax Same as \code{\link{tmlenet}} input argument.
# @param sW Same as \code{\link{tmlenet}} input argument.
# @param sA Same as \code{\link{tmlenet}} input argument.
# @param IDnode (Optional) Same as \code{\link{tmlenet}} input argument.
# @param NETIDnode (Optional) Same as \code{\link{tmlenet}} input argument.
# @param sep Optional friend ID character separator for friends listed in \code{NETIDnode} column of \code{data}, default is \code{' '};
#  same as \code{\link{tmlenet}} input argument \code{optPars$sep}.
# @param NETIDmat (Optional) Same as \code{\link{tmlenet}} input argument.
# @param verbose Set to \code{TRUE} to print messages on status and information to the console.
#  Turn this on by default using \code{options(tmlenet.verbose=TRUE)}.
# @return A named list that contains:
#  \itemize{
#  \item \code{sW.matrix} - Matrix of evaluated summary measures for \code{sW}.
#  \item \code{sA.matrix} - Matrix of evaluated summary measures for \code{sA}.
#  \item \code{NETIDmat} - Network ID matrix that can be used for \code{NETIDmat} input argument to \code{tmlenet}.
#  \item \code{DatNet.ObsP0} - R6 object of class \code{\link{DataStore}} that stores all the summary measures and the network information.
#    This object be passed to \code{\link{tmlenet}} as an argument, in which case the arguments already provided to \code{eval.summaries} no
#    longer need to be specified to \code{tmlenet}.
#  }
# @seealso \code{\link{tmlenet}} for estimation of network effects and \code{\link{def_sW}} for defining the summary measures.
# @example tests/examples/3_eval.summaries_examples.R
# @export
# eval.summaries <- function(data, Kmax, sW, sA, IDnode = NULL, NETIDnode = NULL, sep = ' ', NETIDmat = NULL,
#                             verbose = getOption("tmlenet.verbose")) {
#   iid_data_flag <- FALSE  # set to true if no network is provided (will run iid TMLE)
#   nFnode = "nF"
#   #----------------------------------------------------------------------------------
#   # SOME INPUT CHECKS
#   #----------------------------------------------------------------------------------
#   assert_that(is.data.frame(data))
#   assert_that(is.integerish(Kmax))
#   Kmax <- as.integer(Kmax)
#   assert_that(is.DefineSummariesClass(sW))
#   assert_that(is.DefineSummariesClass(sA))
#   nobs <- nrow(data)

#   # Check no factors exist in the input data:
#   check1 <- CheckExistFactors(data)
#   if (!is.null(check1)) stop("found factor column(s) in the input data, consider removing or recoding such column(s) as strings: "
#                             %+% paste0(check1, collapse=','))

#   if (is.null(NETIDnode) && is.null(NETIDmat)) {
#     message("No network (friends) specified by NETIDnode or NETIDmat args, assuming the input data is i.i.d.")
#     nFnode <- NULL
#     iid_data_flag <- TRUE
#     if (missing(Kmax)) Kmax <- 1 # need to allow Kmax = 0
#   }

#   #---------------------------------------------------------------------------------
#   # Saving the observed data as a data.table in its own object class OdataDT
#   #---------------------------------------------------------------------------------
#   OdataDT_R6 <- OdataDT$new(Odata = data, nFnode = nFnode, iid_data_flag = iid_data_flag)

#   #----------------------------------------------------------------------------------
#   # Create an object with model estimates, data & network information that is passed on to estimation algorithm(s)
#   #----------------------------------------------------------------------------------
#   netind_cl <- NetIndClass$new(nobs = nobs, Kmax = Kmax)
#   if (!is.null(NETIDnode)) {
#     assert_that(is.character(NETIDnode))
#     # Net_str <- as.character(data[, NETIDnode])
#     Net_str <- as.character(OdataDT_R6$OdataDT[[NETIDnode]])
#     OdataDT_R6$OdataDT[, NETIDnode:=NULL]
#     if (!is.null(IDnode)) {
#       assert_that(is.character(IDnode))
#       # IDs_str <- as.character(data[, IDnode])
#       IDs_str <- as.character(OdataDT_R6$OdataDT[[IDnode]])
#       OdataDT_R6$OdataDT[, IDnode:=NULL]
#     } else {
#       IDs_str <- NULL
#     }
#     netind_cl$makeNetInd.fromIDs(Net_str = Net_str, IDs_str = IDs_str, sep = sep)
#   } else if (!is.null(NETIDmat)) {
#     assert_that(is.matrix(NETIDmat))
#     netind_cl$NetInd <- NETIDmat
#     netind_cl$make.nF()
#   }

#   if (verbose) {
#     message("evaluated the network ID matrix: "); print(head(netind_cl$NetInd))
#     message("evaluated and added to summary measures the number of friends for each observation (nF): "); print(head(netind_cl$nF))
#   }

#   #----------------------------------------------------------------------------------
#   # Test parsing and evaluating the summary measures (in class DefineSummariesClass):
#   #----------------------------------------------------------------------------------
#   # Testing the evaluation of summary measures:
#   # # sW.matrix <- sW$eval.nodeforms(data.df = data, netind_cl = netind_cl)
#   # sW.matrix <- sW$eval.nodeforms(data.df = OdataDT_R6$OdataDT, netind_cl = netind_cl)
#   # # sA.matrix <- sA$eval.nodeforms(data.df = data, netind_cl = netind_cl)
#   # sA.matrix <- sA$eval.nodeforms(data.df = OdataDT_R6$OdataDT, netind_cl = netind_cl)
#   # if (verbose) {
#   #   print("sample matrix of sW summary measurs: : "); print(head(sW.matrix))
#   #   print("sample matrix of sA summary measurs: "); print(head(sA.matrix))
#   #   print("map of sW names to column names: "); print(sW$sVar.names.map)
#   #   print("map of sA names to column names: "); print(sA$sVar.names.map)
#   # }

#   #---------------------------------------------------------------------------------
#   # BUILDING OBSERVED sW & sA: (obsdat.sW - a dataset (matrix) of n observed summary measures sW)
#   # ALL EVALUATIONS HAPPEN BY REFERENCE AND ONLY OdataDT_R6$OdataDT is modified, no copies are made
#   #---------------------------------------------------------------------------------
#   datnetW <- DatNet$new(netind_cl = netind_cl, nFnode = nFnode)
#   datnetW$make.sVar(Odata = OdataDT_R6, sVar.object = sW)
#   datnetW$fixmiss_sVar() # permanently replace NA values in sW with 0
#   datnetA <- DatNet$new(netind_cl = netind_cl)
#   datnetA$make.sVar(Odata = OdataDT_R6, sVar.object = sA)
#   DatNet.ObsP0 <- DataStore$new(Odata = OdataDT_R6, datnetW = datnetW, datnetA = datnetA)$make.dat.sWsA()
#   OdataDT_R6$sW <- sW
#   OdataDT_R6$sA <- sA
#   return(list(NETIDmat = netind_cl$NetInd, DatNet.ObsP0 = DatNet.ObsP0))
# }

# tmlenet <- function(DatNet.ObsP0, data, Kmax, sW, sA,
#                     intervene1.sA, f_gstar1 = NULL, intervene2.sA = NULL, f_gstar2 = NULL,
#                     # Anodes,
#                     Ynode,
#                     Qform = NULL, hform.g0 = NULL, hform.gstar = NULL,
#                     # estimators = c("tmle", "iptw", "gcomp"),
#                     # AnodeDET = NULL,
#                     # YnodeDET = NULL,
#                     NETIDmat = NULL,
#                     IDnode = NULL, NETIDnode = NULL,
#                     verbose = getOption("tmlenet.verbose"),
#                     optPars = list(
#                       bootstrap.var = FALSE,
#                       n.bootstrap = 100,
#                       boot.nodes = NULL,
#                       boot.form = NULL,
#                       alpha = 0.05,
#                       lbound = 0.005,
#                       family = "binomial", # NOT YET IMPLEMENTED
#                       # n_MCsims = 1,
#                       # n_MCsims = ifelse(!missing(data),ceiling(sqrt(nrow(data))),10),
#                       runTMLE = c("tmle.intercept", "tmle.covariate"),
#                       YnodeDET = NULL,
#                       # f_gstar2 = NULL,
#                       sep = ' ',
#                       f_g0 = NULL,
#                       h_g0_SummariesModel = NULL,
#                       h_gstar_SummariesModel = NULL)
#                     ) {

#   oldverboseopt <- getOption("tmlenet.verbose")
#   options(tmlenet.verbose = verbose)
#   gvars$verbose <- verbose
#   n_MCsims <- 1
#   #----------------------------------------------------------------------------------
#   # ADDITIONAL ARGUMENTS (removed from input args of tmlenet())
#   #----------------------------------------------------------------------------------
#   AnodeDET <- NULL # removed from the input, since this is not implemented
#   iid_data_flag <- FALSE  # set to true if no network is provided (usual iid TMLE)
#   Q.SL.library <- c("SL.glm", "SL.step", "SL.glm.interaction") # NOT USED
#   g.SL.library <- c("SL.glm", "SL.step", "SL.glm.interaction") # NOT USED
#   max_npwt <- 50 # NOT USED YET
#   h_logit_sep_k <- FALSE # NOT USED YET
#   bootstrap.var <- ifelse(is.null(optPars$bootstrap.var), FALSE, optPars$bootstrap.var)
#   n.bootstrap <- ifelse(is.null(optPars$n.bootstrap), 100, optPars$n.bootstrap)
#   boot.nodes <- if(is.null(optPars$boot.nodes)) {NULL} else {optPars$boot.nodes}
#   boot.form <- if(is.null(optPars$boot.form)) {NULL} else {optPars$boot.form}
#   alpha <- ifelse(is.null(optPars$alpha), 0.05, optPars$alpha)
#   lbound <- ifelse(is.null(optPars$lbound), 0.005, optPars$lbound)
#   family <- ifelse(is.null(optPars$family), "binomial", optPars$family)
#   sep <- ifelse(is.null(optPars$sep), ' ', optPars$sep)
#   assert_that(is.character(sep) && length(sep)==1L)
#   n_MCsims <- ifelse(is.null(optPars$n_MCsims), 1, optPars$n_MCsims)
#   f_g0 <- if(is.null(optPars$f_g0)) {NULL} else {optPars$f_g0}
#   if (!is.null(f_g0)) assert_that(is.function(f_g0))
#   YnodeDET <- if(is.null(optPars$YnodeDET)) {NULL} else {optPars$YnodeDET}
#   # f_gstar2 <- if(is.null(optPars$f_gstar2)) {NULL} else {optPars$f_gstar2}
#   # if (!is.null(f_gstar2)) assert_that(is.function(f_gstar2) || is.vector(f_gstar2))
#   # if TRUE, only evaluate the intercept-based TMLE (TMLE_B), if FALSE, evaluate only the covariate-based TMLE (TMLE_A)
#   runTMLE <- optPars$runTMLE[1]
#   if (is.null(runTMLE) || (runTMLE[1] %in% "tmle.intercept")) {
#     onlyTMLE_B <- TRUE
#   } else if (runTMLE %in% "tmle.covariate") {
#     onlyTMLE_B <- FALSE
#   } else {
#     stop("optPars[['runTMLE']] argument must be either 'tmle.intercept' or 'tmle.covariate'")
#   }

#   if (verbose) {
#     message("Running tmlenet with the following settings from tmlenet_options(): ");
#     str(gvars$opts)
#     message("Running tmlenet with the following settings from optPars arg of tmlenet(): ");
#     str(optPars)
#   }

#   #----------------------------------------------------------------------------------
#   # DETERMINING INTERNAL / EXTERNAL ESTIMATOR NAMES THAT WILL BE EVALUATED
#   #----------------------------------------------------------------------------------
#   # onlyTMLE_B <- TRUE
#   assert_that(assertthat::is.flag(onlyTMLE_B))
#   estnames.internal <- c("TMLE_A", "TMLE_B", "h_IPTW", "MLE")
#   # estnames.internal <- c("TMLE_A", "TMLE_B", "TMLE_g_IPTW", "h_IPTW", "g_IPTW", "MLE")
#   names(estnames.internal) <- estnames.internal
#   estnames.out <- c("tmle", "h_iptw", "gcomp")
#   if (onlyTMLE_B) {
#     estnames.internal <- estnames.internal[-which(estnames.internal%in%"TMLE_A")]
#   } else {
#     estnames.internal <- estnames.internal[-which(estnames.internal%in%"TMLE_B")]
#   }
#   names(estnames.out) <- estnames.internal
#   estnames.internal <- as.list(estnames.internal)

#   #----------------------------------------------------------------------------------
#   # MONTE-CARLO SIMULATION PARAMETERS
#   #----------------------------------------------------------------------------------
#   nQ.MCsims <- as.integer(n_MCsims)  # number of times to sample MC sim for Q (each of size n)
#   ng.MCsims <- as.integer(n_MCsims)  # number of times to sample MC sim for h (each of size n)
#   assert_that(is.count(nQ.MCsims))
#   assert_that(is.count(ng.MCsims))
#   max.err_est <- 0.1    # maximum percent error for MCS estimators

#   #----------------------------------------------------------------------------------
#   # PARAMETERS FOR ESTIMATING h under g0 & gstar
#   #----------------------------------------------------------------------------------
#   f.g0 <- f_g0
#   f.g0_args <- NULL
#   h_user_fcn <- NULL
#   h_user <- !(is.null(h_user_fcn))

#   #----------------------------------------------------------------------------------
#   # Perform some input checks
#   # Create an object with data & network information that is passed on to estimation algorithm(s)
#   # Build observed sW & sA: (obsdat.sW - a dataset (matrix) of n observed summary measures sW)
#   #----------------------------------------------------------------------------------
#   if (missing(Ynode)) Ynode <- NULL
#   if (is.null(Ynode) && is.null(Qform)) stop("Either Qform or Ynode must be specified")
#   if (!is.null(Qform) && !is.null(Ynode)) {
#     message("Since both Ynode and Qform are specified, the left-hand side of Qform will be ignored, with outcome being set to Ynode: " %+%
#       Ynode)
#   }
#   if (!is.null(Qform) && is.null(Ynode)) {
#     Ynode <- LhsVars(Qform)[1]
#     message("Setting the Ynode to: " %+% Ynode)
#   }

#   if (missing(DatNet.ObsP0)) {
#     # time_evalsumm <- system.time(
#       eval.summ.res <- eval.summaries(data = data, Kmax = Kmax, sW = sW, sA = sA,
#                                       IDnode = IDnode, NETIDnode = NETIDnode,
#                                       sep = sep, NETIDmat = NETIDmat, verbose = FALSE)
#       # )
#     # print("time eval.summaries(...): "); print(time_evalsumm)
#     DatNet.ObsP0 <- eval.summ.res$DatNet.ObsP0
#     data <- DatNet.ObsP0$datnetW$Odata
#     sW <- data$sW
#     sA <- data$sA
#     Kmax <- DatNet.ObsP0$Kmax
#     # sW <- DatNet.ObsP0$datnetW$sVar.object
#     # sA <- DatNet.ObsP0$datnetA$sVar.object
#   }

#   intervene1.sA <- update.intervention.sA(new.sA = intervene1.sA, sA = sA)
#   data$intervene1.sA <- intervene1.sA
#   if (!is.null(intervene2.sA)) {
#     intervene2.sA <- update.intervention.sA(new.sA = intervene2.sA, sA = sA)
#     data$intervene2.sA <- intervene2.sA
#     # making sure exactly the same Anodes are used for both interventions:
#     if (any(intervene1.sA$Anodes != intervene2.sA$Anodes))
#       stop("intervention node names have to be identical for intervene1.sA and intervene2.sA")
#   }

#   # new version of nodes:
#   node_l <- list(nFnode = DatNet.ObsP0$datnetW$nFnode,
#                   Anodes = intervene1.sA$Anodes,  # Note that Anodes might be different in intervene1.sA & intervene2.sA
#                   AnodeDET = AnodeDET,
#                   Ynode = Ynode, YnodeDET = YnodeDET)

#   data$nodes <- node_l
#   data$backupAnodes(Anodes = node_l$Anodes, sA = sA)

#   nobs <- DatNet.ObsP0$nobs
#   DatNet.ObsP0$nodes <- node_l

#   #----------------------------------------------------------------------------------
#   # Defining (and checking) Deterministic Y and A node flags:
#   #----------------------------------------------------------------------------------
#   if (is.null(AnodeDET)) {
#     determ.g <- rep_len(FALSE, nobs)
#   } else {
#     CheckVarNameExists(data$OdataDT, AnodeDET)
#     # determ.g <- (data[, AnodeDET] == 1)
#     setkeyv(data$OdataDT, cols=AnodeDET)
#     data$OdataDT[list(1),AnodeDET, with=FALSE]
#   }
#   if (is.null(YnodeDET)) {
#     determ.Q <- rep_len(FALSE, nobs)
#   } else {
#     CheckVarNameExists(data$OdataDT, YnodeDET)
#     # determ.Q <- (data[, YnodeDET] == 1)
#     setkeyv(data$OdataDT, cols=YnodeDET)
#     data$OdataDT[list(1),YnodeDET, with=FALSE]
#   }

#   for (Anode in node_l$Anodes) CheckVarNameExists(data$OdataDT, Anode)
#   # CheckVarNameExists(data, node_l$Anodes)
#   for (Ynode in node_l$Ynode) CheckVarNameExists(data$OdataDT, Ynode)
#   # CheckVarNameExists(data, node_l$Ynode)

#   #----------------------------------------------------------------------------------
#   # NOTE: YnodeVals = obsYvals, det.Y = determ.Q need to be added to DatNet.ObsP0 after returned from eval.summaries()
#   #----------------------------------------------------------------------------------
#   obsYvals <- data$OdataDT[[node_l$Ynode]]
#   DatNet.ObsP0$addYnode(YnodeVals = obsYvals, det.Y = determ.Q)

#   #----------------------------------------------------------------------------------
#   # (OPTIONAL) ADDING DETERMINISTIC/DEGENERATE Anode FLAG COLUMNS TO sA:
#   # cancelled adding DET nodes to sVar since all sVar are automatically get added to A ~ predictors + DETnodes...
#   #----------------------------------------------------------------------------------
#   # obsdat.sW <- O.datnetW$add_deterministic(Odata = data, userDETcol = "determ.g")$dat.sVar
#   # Testing NA for visible det.Y and true observed Y as protected:
#   # determ.Q <- c(FALSE, FALSE, FALSE, rep(TRUE, length(determ.Q)-3))
#   # length(determ.Q) == length(obsYvals)
#   # DatNet.ObsP0 <- DataStore$new(datnetW = datnetW, datnetA = datnetA, YnodeVals = obsYvals, det.Y = determ.Q)$make.dat.sWsA()
#   # print("DatNet.ObsP0: "); print(DatNet.ObsP0)
#   # print(head(cbind(DatNet.ObsP0$noNA.Ynodevals, DatNet.ObsP0$YnodeVals, data[,node_l$Ynode]), 100))

#   #----------------------------------------------------------------------------------
#   # Optional regressions specs:
#   #----------------------------------------------------------------------------------
#   # Q.sVars <- process_regform(as.formula(Qform), sW.map = c(sW$sVar.names.map, sA$sVar.names.map), sA.map = node_l$Ynode)
#   Q.sVars <- process_regforms(Qform, sW.map = c(sW$sVar.names.map, sA$sVar.names.map), sA.map = node_l$Ynode)
#   # h.g0.sVars <- process_regform(as.formula(hform.g0[[1]]), sW.map = sW$sVar.names.map, sA.map = sA$sVar.names.map)
#   h.g0.sVars <- process_regforms(hform.g0, sW.map = sW$sVar.names.map, sA.map = sA$sVar.names.map)

#   if (!is.null(hform.gstar)) {
#     # h.gstar.sVars <- process_regform(as.formula(hform.gstar), sW.map = sW$sVar.names.map, sA.map = sA$sVar.names.map)
#     h.gstar.sVars <- process_regforms(hform.gstar, sW.map = sW$sVar.names.map, sA.map = sA$sVar.names.map)
#   } else {
#     h.gstar.sVars <- h.g0.sVars
#   }

#   if (verbose) {
#     print("Input regression(s) Qform (E(Y|sA,sW)): "); res <- lapply(Qform, function(Qform) print(as.formula(Qform),showEnv=FALSE))
#     print("Derived regression(s) from Qform:"); str(Q.sVars)
#     print("Input regression(s) hform.g0 (P(sA|sW) under g0): "); res <- lapply(hform.g0, function(hform.g0) print(as.formula(hform.g0),showEnv=FALSE))
#     print("Derived regression(s) from hform.g0: "); str(h.g0.sVars)
#     print("Input regression(s) hform.gstar (P(sA|sW) under g.star): "); res <- lapply(hform.gstar, function(hform.gstar) print(as.formula(hform.gstar),showEnv=FALSE))
#     print("Derived regression(s) from hform.gstar: "); str(h.gstar.sVars)
#   }

#   #-----------------------------------------------------------
#   # Defining and fitting regression for Y ~ sW + sA:
#   #-----------------------------------------------------------
#   check.Qpreds.exist <- unlist(lapply(Q.sVars$predvars, function(PredName) PredName %in% DatNet.ObsP0$names.sVar))
#   if (!all(check.Qpreds.exist)) stop("the following predictors in Qform regression could not be located among the summary measures: " %+%
#                                     paste0(Q.sVars$predvars[!check.Qpreds.exist], collapse = ","))

#   if (verbose) {
#     for (idx in seq_along(node_l$Ynode)) {
#       message("================================================================")
#       message("fitting E(Y|sA,sW):= ", "P(" %+% node_l$Ynode[idx] %+% "=1 | " %+% paste(Q.sVars$predvars[[idx]], collapse = ",") %+% ")")
#       message("================================================================")
#     }
#   }

#   Qreg <- RegressionClass$new(outvar = node_l$Ynode,
#                               predvars = Q.sVars$predvars[[1]],
#                               subset = !determ.Q, ReplMisVal0 = TRUE)

#   m.Q.init <- BinOutModel$new(glm = FALSE, reg = Qreg)$fit(data = DatNet.ObsP0)

#   #----------------------------------------------------------------------------------
#   # Create an object with model estimates, data & network information that is passed on to estimation procedure
#   #----------------------------------------------------------------------------------
#   # 1) define parameters for MC estimation of the substitution estimators
#   # 2) define parameters for estimation of the efficient weights h(A^s|W^s)
#   est_obj <- list(
#                   estnames = estnames.internal,
#                   lbound = lbound[1],
#                   max.err_eps = max.err_est,  # error tolerance for the mean/var M.C. estimate
#                   m.Q.init = m.Q.init,
#                   f.g0 = f.g0,
#                   sW = sW,
#                   sA = sA,
#                   Q.sVars = Q.sVars,
#                   h.g0.sVars = h.g0.sVars,
#                   h.gstar.sVars = h.gstar.sVars,
#                   nQ.MCsims = nQ.MCsims,
#                   ng.MCsims = ng.MCsims,
#                   h_g0_SummariesModel = optPars$h_g0_SummariesModel,
#                   h_gstar_SummariesModel = optPars$h_gstar_SummariesModel,
#                   # Cap the prop weights scaled at max_npwt (for =50 -> translates to max 10% of total weight for n=500 and 5% for n=1000):
#                   max_npwt = max_npwt # NOT IMPLEMENTED
#                   # h_logit_sep_k = h_logit_sep_k, # NOT IMPLEMENTED
#                   # h_user = h_user, # NOT IMPLEMENTED
#                   # h_user_fcn = h_user_fcn, # NOT IMPLEMENTED
#                   )

#   est_obj_g1 <- append(est_obj,
#                       list(
#                         f.gstar = f_gstar1,
#                         new.sA = intervene1.sA
#                         )
#                       )

#   if (!is.null(intervene2.sA) || !is.null(f_gstar2)) {
#     est_obj_g2 <- append(est_obj,
#                       list(
#                         f.gstar = f_gstar2,
#                         new.sA = intervene2.sA
#                         )
#                       )
#   }

#   #----------------------------------------------------------------------------------
#   # Running MC evaluation for substitution TMLE ests
#   #----------------------------------------------------------------------------------
#   tmle_g1_out <- get_all_ests(estnames = estnames.internal, DatNet.ObsP0 = DatNet.ObsP0, est_params_list = est_obj_g1)

#   if (!is.null(intervene2.sA) || !is.null(f_gstar2)) {
#     tmle_g2_out <- get_all_ests(estnames = estnames.internal, DatNet.ObsP0 = DatNet.ObsP0, est_params_list = est_obj_g2)
#   } else {
#     tmle_g2_out <- NULL
#   }

#   iidEIC.eval <- TRUE
#   if (iidEIC.eval) {
#     # MC.tmle.eval.t <- system.time(
#       MC.tmle.eval <- MCeval_fWi(n.MC = n.bootstrap, DatNet.ObsP0 = DatNet.ObsP0, tmle_g1_out = tmle_g1_out, tmle_g2_out = tmle_g2_out)
#     # )
#     # tmle_g1_out$ests_mat["TMLE",] <- mean(MC.tmle.eval$EY_gstar1)
#     # if (!is.null(tmle_g2_out)) tmle_g2_out$ests_mat["TMLE",] <- mean(MC.tmle.eval$EY_gstar2)
#   } else {
#     MC.tmle.eval <- list(EY_gstar1 = NA, EY_gstar2 = NA, ATE = NA)
#   }
#   # print("MC.tmle.eval.t"); print(MC.tmle.eval.t)

#   if (bootstrap.var) {
#     # ------------------------------------------------------------------------------------------
#     # IID BOOSTRAP FOR THE TMLE (DOES'T WORK FOR DEP DATA):
#     # ------------------------------------------------------------------------------------------
#     # var_tmleB_boot <- iid_bootstrap_tmle(n.boot, estnames, DatNet.ObsP0, tmle_g_out, QY_mat, wts_mat)

#     # ------------------------------------------------------------------------------------------
#     # PARAMETRIC BOOSTRAP TMLE variance estimate:
#     # ------------------------------------------------------------------------------------------
#     var_tmleB_boot <- par_bootstrap_tmle(n.boot = n.bootstrap, boot.nodes = boot.nodes, boot.form = boot.form,
#                                          estnames = estnames.internal, DatNet.ObsP0 =  DatNet.ObsP0,
#                                          tmle_g1_out = tmle_g1_out, tmle_g2_out = tmle_g2_out)
#   } else {
#     var_tmleB_boot <- list(EY_gstar1 = NA, EY_gstar2 = NA, ATE = NA)
#   }

#   #----------------------------------------------------------------------------------
#   # Create output list (estimates, as. variances, CIs)
#   #----------------------------------------------------------------------------------
#   EY_gstar1 <- make_EYg_obj(estnames = estnames.internal, estoutnames = estnames.out, alpha = alpha,
#                             DatNet.ObsP0 = DatNet.ObsP0, tmle_g_out = tmle_g1_out,
#                             MC.tmle.eval = MC.tmle.eval$EY_gstar1,
#                             var_tmleB_boot = var_tmleB_boot$EY_gstar1)

#   EY_gstar2 <- NULL; ATE <- NULL

#   if (!is.null(intervene2.sA) || !is.null(f_gstar2)) {
#     EY_gstar2 <- make_EYg_obj(estnames = estnames.internal, estoutnames = estnames.out, alpha = alpha,
#                               DatNet.ObsP0 = DatNet.ObsP0, tmle_g_out=tmle_g2_out,
#                               MC.tmle.eval = MC.tmle.eval$EY_gstar2,
#                               var_tmleB_boot = var_tmleB_boot$EY_gstar2)

#     ATE <- make_EYg_obj(estnames = estnames.internal, estoutnames = estnames.out, alpha = alpha,
#                         DatNet.ObsP0 = DatNet.ObsP0, tmle_g_out = tmle_g1_out, tmle_g2_out = tmle_g2_out,
#                         MC.tmle.eval = MC.tmle.eval$ATE,
#                         var_tmleB_boot = var_tmleB_boot$ATE)
# 	}

# 	tmlenet.res <- list(EY_gstar1 = EY_gstar1, EY_gstar2 = EY_gstar2, ATE = ATE)
# 	class(tmlenet.res) <- c(class(tmlenet.res), "tmlenet")

#   options(tmlenet.verbose = oldverboseopt)
#   gvars$verbose <- oldverboseopt
# 	return(tmlenet.res)
# }