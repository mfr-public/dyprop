#' @include dyprop_init.R
NULL

#' The dyprop Class
#'
#' An S4 class to represent Dynamic Proportionality results, extending the
#' base \code{propr} class to the temporal domain.
#'
#' @slot pseudotime A numeric vector. The ordered pseudotime or time-series coordinates
#'  corresponding to the columns of the count matrix.
#' @slot Phi_Array A 3D array (Genes x Genes x Timepoints). Stores the instantaneous
#'  Dynamic Instability metric (\eqn{\Phi(t)}) calculated by the vectorized scanner.
#' @slot Rho_Array A 3D array (Genes x Genes x Timepoints). Stores the instantaneous
#'  Dynamic Coupling metric (\eqn{\rho(t)}).
#' @slot events A data.frame. The catalog of detected transitions, containing
#'  columns for \code{GeneA}, \code{GeneB}, \code{Event_Class} (Switch/Decoupling),
#'  \code{Tau}, \code{Epsilon}, and \code{FDR}.
#' @slot glmm_fits A list. Stores the post-hoc Generalized Linear Mixed Model objects
#'  for the subset of features validated in Phase III.
#' @slot dictionary_meta A list. Stores metadata about the basis function dictionary
#'  used for scanning (e.g., \code{epsilon_min}, \code{grid_resolution}).
#'
#' @inheritParams propr::propr
#' @seealso \code{\link[propr]{propr}}
#'
#' @exportClass dyprop
setClass(
  "dyprop",
  contains = "propr",
  slots = c(
    pseudotime = "numeric",
    Phi_Array = "array",
    Rho_Array = "array",
    events = "data.frame",
    glmm_fits = "list",
    dictionary_meta = "list"
  ),
  prototype = list(
    pseudotime = numeric(0),
    Phi_Array = array(0, dim = c(0, 0, 0)),
    Rho_Array = array(0, dim = c(0, 0, 0)),
    events = data.frame(),
    glmm_fits = list(),
    dictionary_meta = list()
  )
)

#' Validity Check for dyprop
#'
#' @param object A dyprop object.
#' @return TRUE if valid, otherwise a character string describing the error.
setValidity("dyprop",
  function(object) {
    msg <- NULL
    
    # 1. Inherited Dimension Check
    # Ensure the new arrays match the dimensions of the parent logratio matrix
    n_features <- ncol(object@logratio)
    n_samples <- nrow(object@logratio)
    
    # Check Phi_Array Dimensions (Genes x Genes x Time)
    if (any(dim(object@Phi_Array)[1:2] != n_features) && length(object@Phi_Array) > 0) {
      msg <- c(msg, "Phi_Array dimensions [1:2] do not match number of features in @logratio.")
    }

    # Check Rho_Array Dimensions
    if (any(dim(object@Rho_Array)[1:2] != n_features) && length(object@Rho_Array) > 0) {
      msg <- c(msg, "Rho_Array dimensions [1:2] do not match number of features in @logratio.")
    }
    
    # 2. Pseudotime Consistency
    # The length of pseudotime vector must match the number of samples (rows in counts)
    if (length(object@pseudotime) > 0 && length(object@pseudotime) != n_samples) {
      msg <- c(msg, "Length of @pseudotime vector does not match number of samples (rows) in @counts.")
    }
    
    # 3. Events Structure
    # Ensure events dataframe has the mandatory columns for the Decision Tree
    required_cols <- c("GeneA", "GeneB", "Event_Class", "Tau", "Epsilon")
    if (nrow(object@events) > 0 && !all(required_cols %in% colnames(object@events))) {
      msg <- c(msg, paste("slot @events missing required columns:", paste(required_cols, collapse=", ")))
    }

    if (is.null(msg)) TRUE else msg
  }
)