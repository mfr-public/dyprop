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
#' @seealso \code{\link[propr]{propr}}
#'
#' @exportClass dyprop
setClass(
  "dyprop",
  contains = "propr",
  slots = c(
    pseudotime = "numeric",
    design = "data.frame",
    metric_cache = "environment",
    events = "data.frame",
    glmm_fits = "list",
    dictionary_meta = "list",
    fdr_cutoff = "list"
  ),
  prototype = list(
    pseudotime = numeric(0),
    design = data.frame(),
    events = data.frame(),
    glmm_fits = list(),
    dictionary_meta = list(),
    fdr_cutoff = list()
  )
)

#' Validity Check for dyprop
#' @name dyprop-validity
#'
#' @param object A dyprop object.
#' @return TRUE if valid, otherwise a character string describing the error.
setValidity(
  "dyprop",
  function(object) {
    msg <- NULL

    # 1. Inherited Dimension Check
    n_features <- ncol(object@logratio)
    n_samples <- nrow(object@logratio)

    # 2. Pseudotime Consistency
    # The length of pseudotime vector must match the number of samples (rows in counts)
    if (length(object@pseudotime) > 0 && length(object@pseudotime) != n_samples) {
      msg <- c(msg, "Length of @pseudotime vector does not match number of samples (rows) in @counts.")
    }

    # 3. Events Structure
    # Ensure events dataframe has the mandatory columns for the Decision Tree
    required_cols <- c("GeneA", "GeneB", "Event_Class", "Tau", "Epsilon")
    if (nrow(object@events) > 0 && !all(required_cols %in% colnames(object@events))) {
      msg <- c(msg, paste("slot @events missing required columns:", paste(required_cols, collapse = ", ")))
    }

    if (is.null(msg)) TRUE else msg
  }
)
