#' sl3_wrapper_logisfitR6
#'
#' @docType class
#' @export
sl3_wrapper_logisfitR6 <- R6Class("sl3_wrapper_logisfitR6",
  inherit = logisfitR6,
  public = list(
    lmclass = NULL,
    fitfunname = NULL,
    # sl3_lrnr = NULL,

    initialize = function(sl3_lrnr) {
      assert_that(is(sl3_lrnr,"Lrnr_base"))
      # self$fitfunname <- sl3_lrnr$name
      self$lmclass <- paste0("sl3::", class(sl3_lrnr)[1L])
      self$fitfunname <- paste0("sl3::", class(sl3_lrnr)[1L], "$train")
      private$sl3_lrnr <- sl3_lrnr
      self$get_validity
    },

    fit = function(datsum_obj) {
      if (gvars$verbose) print(paste("calling glm.generic for", self$fitfunname))
      X_mat <- datsum_obj$getXDT
      Y_vals <- datsum_obj$getY
      dataDT <- cbind(X_mat, Y = Y_vals)
      sl3_lrnr <- private$sl3_lrnr
      if (nrow(dataDT) > 0) {
        task <- sl3::sl3_Task$new(dataDT, covariates=colnames(X_mat), outcome=colnames(dataDT)[ncol(dataDT)])
        sl3_lrnr <- sl3_lrnr$train(task)
      }
      if (gvars$verbose) print(sl3_lrnr)
      return(sl3_lrnr)
    },

    predict.long = function(datsum_obj, m.fit) {
     stop("not implemented")
    },

    predict = function(datsum_obj, m.fit) {
      if (gvars$verbose) print(paste("calling predict for", self$fitfunname))
      X_mat <- datsum_obj$getXDT
      assert_that(!is.null(X_mat)); assert_that(!is.null(datsum_obj$subset_idx))
      pAout <- rep.int(gvars$misval, datsum_obj$n)
      if (sum(datsum_obj$subset_idx > 0)) {
        new_task <- sl3::sl3_Task$new(X_mat, covariates=colnames(X_mat), outcome=NULL)
        ## Learner hasn't been trained,
        ## means we are making predictionsin for a last degenerate bin.
        ## These predictions play no role and aren't used.
        if (is(m.fit, "Lrnr_base") && !m.fit$is_trained) {
          pAout[datsum_obj$subset_idx] <- 0.5
        } else {
          pAout[datsum_obj$subset_idx] <- m.fit$predict(new_task)[[1]]
        }
      }
      return(pAout)
    }
  ),
private = list(
  sl3_lrnr = NULL
  )
)