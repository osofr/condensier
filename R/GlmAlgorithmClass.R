#' logisfitR6
#'
#' @docType class
#' @export
logisfitR6 <- R6Class("logisfitR6",
  public =
    list(
         lmclass = NULL,
         fitfunname = NULL,
         initialize = function() {
           self$get_validity
         },

         predict.long = function(datsum_obj, m.fit) {
           if (gvars$verbose) print(paste("calling predict.long for", self$fitfunname))
           X_mat <- datsum_obj$getXmat
           Y_vals <- datsum_obj$getY
           assert_that(!is.null(datsum_obj));assert_that(!is.null(X_mat));
           assert_that(!is.null(datsum_obj$subset_idx))
           assert_that(nrow(X_mat)==length(Y_vals))

           pAout <- rep.int(gvars$misval, datsum_obj$n)
           if (sum(datsum_obj$subset_idx > 0)) {
             # -----------------------------------------------------------------
             # OBTAINING PREDICTIONS FOR LONG FORMAT P(Ind_j = 1 | Bin_j, W) BASED ON EXISTING POOLED FIT:
             # -----------------------------------------------------------------
             probA1 <- private$do.predict(X_mat = X_mat, m.fit = m.fit)
             # -----------------------------------------------------------------
             # GETTING ID-BASED PREDICTIONS (n) as cumprod of P(Ind_j = 1 | Bin_j, W) for j = 1, ..., K
             # -----------------------------------------------------------------
             ProbAeqa_long <- as.vector(probA1^(Y_vals) * (1L - probA1)^(1L - Y_vals))
             res_DT <- data.table(ID = datsum_obj$ID, ProbAeqa_long = ProbAeqa_long)
             res_DT <- res_DT[, list(cumprob = cumprod(ProbAeqa_long)), by = ID]
             data.table::setkeyv(res_DT, c("ID")) # sort by ID
             res_DT_short <- res_DT[unique(res_DT[, key(res_DT), with = FALSE]), mult = 'last']
             ProbAeqa <- res_DT_short[["cumprob"]]
             # print("res_DT: "); print(res_DT)
             # print("res_DT w/ cumprob: "); print(res_DT)
             # print("res_DT_short: "); print(res_DT_short)
             # print("length(ProbAeqa): " %+% length(ProbAeqa))
             # print("head(ProbAeqa, 50)"); print(head(ProbAeqa, 50))
             pAout[datsum_obj$subset_idx] <- ProbAeqa
           }
           return(pAout)
         },

         # Generic prediction fun for logistic regression coefs, predicts P(A = 1 | newXmat)
         # No need for S3 for now, until need different pred. funs for different classes
         # Does not handle cases with deterministic Anodes in the original data..
         predict = function(datsum_obj, m.fit) {
            if (gvars$verbose) print(paste("calling predict for", self$fitfunname))
            X_mat <- datsum_obj$getXmat
            assert_that(!is.null(X_mat))
            assert_that(!is.null(datsum_obj$subset_idx))
            assert_that(!is.null(m.fit))
            # Set to default missing value for A[i] degenerate/degerministic/misval:
            # Alternative, set to default replacement val: pAout <- rep.int(gvars$misXreplace, newdatsum_obj$n)
            pAout <- rep.int(gvars$misval, datsum_obj$n)
            if (sum(datsum_obj$subset_idx > 0)) {
              result <- ifelse(m.fit$has_failed,
                               private$default_predict(X_mat, m.fit),
                               private$do.predict(X_mat, m.fit))
              pAout[datsum_obj$subset_idx] <- result
            }
            return(pAout)
         },

         fit = function(datsum_obj, fn = private$do.fit, ...) {
           if (gvars$verbose) print(paste("calling fit / update for", self$fitfunname))
           X_mat <- datsum_obj$getXmat
           Y_vals <- datsum_obj$getY
           has_failed <- FALSE

           # X_mat has 0 rows: return NA's and avoid throwing exception:
           if (nrow(X_mat) == 0L) {
             m.fit <- rep.int(NA_real_, ncol(X_mat))
             has_failed <- TRUE
           } else {
            m.fit <- try(m.fit <- fn(X_mat, Y_vals, ...), silent = TRUE)

            # If an algorithm fails, we fall back to the fallback function (which is generally the glm function.
            # Cases where this happens include empty bins, constant bins, or bins with a single outcome
            if (inherits(m.fit, "try-error")) {
              if (gvars$verbose) message(self$fitfunname, "failed, falling back on stats:glm.fit; ", m.fit)
              m.fit <- private$fallback_function(X_mat, Y_vals)
              has_failed <- TRUE
            }

           }
           fit <- list(coef = m.fit,
                       linkfun = "logit_linkinv",
                       fitfunname = self$fitfunname,
                       has_failed = has_failed)

           if (gvars$verbose) print(fit$coef)
           class(fit) <- c(class(fit), c(self$lmclass))
           return(fit)
         },

         update = function(datsum_obj, m.fit) {
           if (m.fit$has_failed) {
            # TODO: We have to deal with this in some way. If the calling algorithm is online, and the method had
            # failes earlier, this will cause trouble. GLM is not online at the moment, and will therefore just refit.
            # which is pretty subobtimal. 
            self$fit(datsum_obj, fn = private$do.fit)
           }
           self$fit(datsum_obj, fn = private$do.update, m.fit = m.fit)
          }
         ),
  active =
    list(
         get_validity = function() {
           errors = character()
           if(is.null(self$fitfunname)){
             errors <- c(errors, 'Please define a fitfunname')
           }
           if(is.null(self$lmclass)){
             errors <- c(errors, 'Please define a lmclass')
           }

           # Defines fit
           # Defines predict
           if(length(errors)!= 0) stop(errors)
           TRUE
         }
        ),
  private =
    list(
        fallback_function = NULL,

        do.fit = function(X_mat, Y_vals) {
          stop('Override this function in a subclass')
        },

        default_predict = function(X_mat, m.fit) {
          eta <- X_mat[,!is.na(m.fit$coef), drop = FALSE] %*% m.fit$coef[!is.na(m.fit$coef)]
          result = tryCatch({
            match.fun(FUN = m.fit$linkfun)(eta)
          }, error = function(e) {
            # If for some reason we are not able to call the fast function, fall back to the slow function
            if (m.fit$linkfun == 'logit_linkinv') {
              qlogis(eta)
            } else {
              stop(e)
            }
          })
          return(result)
        },

        do.predict = function(X_mat, m.fit) {
          private$default_predict(X_mat, m.fit)
        },

        do.update = function() {
          stop('The used method does not support updating!')
        }

    )
)

#' glmR6
#'
#' @docType class
#' @export
glmR6 <- R6Class("glmR6",
  inherit = logisfitR6,
  public =
    list(
        fitfunname='glm',
        lmclass='glmR6',
        initialize = function() { }

        ),
  active =
    list(
         get_fit_function = function() {
           return(private$do.fit)
         }
        ),
  private =
    list(
        do.fit = function(X_mat, Y_vals) {
          ctrl <- glm.control(trace = FALSE)
          # ctrl <- glm.control(trace = FALSE, maxit = 1000)
          # SuppressGivenWarnings({
          #   return(stats::glm.fit(x = X_mat, y = Y_vals, family = binomial() , control = ctrl))
          # }, GetWarningsToSuppress())
          suppressWarnings({
            return(stats::glm.fit(x = X_mat, y = Y_vals, family = binomial() , control = ctrl)$coef)
          })
        }

        #do.predict = function(X_mat, m.fit) {
          #eta <- X_mat[,!is.na(m.fit$coef), drop = FALSE] %*% m.fit$coef[!is.na(m.fit$coef)]
          #match.fun(FUN = m.fit$linkfun)(eta)
        #}
    )
)

#' speedglmR6
#'
#' @docType class
#' @export
speedglmR6 <- R6Class("speedglmR6",
  inherit = logisfitR6,
  public =
    list(
        fitfunname='speedglm',
        lmclass='speedglmS3',
        initialize = function() { }
        ),
  active =
    list(
        ),
  private =
    list(
        do.fit = function(X_mat, Y_vals) {
          # , maxit=1000
            suppressWarnings(
              m.fit <- speedglm::speedglm.wfit(X = X_mat, y = Y_vals, family = binomial(), trace = FALSE, method='Cholesky')$coef
            )

        }
    )
)
