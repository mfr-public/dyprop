#' @include AllClasses.R
NULL

#' Estimate False Discovery Rate (FDR)
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function was deprecated in dyprop v0.3.0. The "Topological Vector Paradigm"
#' introduced analytical likelihood tests (LRT) via Generalized Linear Mixed Models (GLMMs)
#' to formally compute statistical significance, making empirical null permutation distributions
#' (Circular Phase Shifts) mathematically obsolete and computationally redundant.
#'
#' @param object A \code{dyprop} object.
#' @param ... Additional arguments (ignored).
#'
#' @return The original \code{dyprop} object without modification.
#'
#' @export
estimateFDR <- function(object, ...) {
    warning("estimateFDR() is deprecated as of dyprop v0.3.0. The pipeline now utilizes analytical GLMM Likelihood Ratio Tests (LRT) in validateGLMM() for superior precision without permutation overhead. Returning object unmodified.", call. = FALSE)
    return(object)
}
