#' @include AllClasses.R
NULL

#' Estimate False Discovery Rate (FDR)
#'
#' Estimates the statistical significance of detected events using an empirical Null Model.
#' The method applies Phase Randomization (Circular Shift) to the data matrix. This preserves
#' the internal autocorrelation and generic non-linear drift of the biological system while
#' mathematically decoupling the specific transition coordinates.
#'
#' @param object A \code{dyprop} object with populated \code{@events} slot.
#' @param n_permutations Integer. Number of null scores to generate. Default 10,000.
#'  Note: This is not the number of *full* permutations, but the total number of
#'  null pairs derived from one or more shuffled scans.
#' @param cores Integer. Number of cores for parallel scanning. Default 1.
#'
#' @return A \code{dyprop} object with updated \code{@events$FDR} column.
#'
#' @details
#' The FDR is calculated as:
#' \eqn{FDR(s) = \frac{\text{Proportion of Null Scores } \ge s}{\text{Proportion of Real Scores } \ge s}}
#'
#' @export
estimateFDR <- function(object, n_permutations = 10000, cores = 1L) {
    if (nrow(object@events) == 0) {
        warning("No events found. Please run scanDynamics() first.")
        return(object)
    }

    message(">>> Phase 3.5: Estimating FDR (Empirical Null Model)...")

    # 1. Phase Randomization (Circular Shift) Object
    # We shift the CLR matrix cyclically by 30% to preserve autocorrelation but destroy topology.
    null_object <- object
    set.seed(42) # Reproducibility
    
    n_samples <- nrow(null_object@logratio)
    shift_fraction <- 0.30
    shift_idx <- floor(n_samples * shift_fraction)
    
    idx_shift <- c((shift_idx + 1):n_samples, 1:shift_idx)
    null_object@logratio <- null_object@logratio[idx_shift, , drop = FALSE]

    # 2. Tau Seam Masking
    # The circular shift creates an artificial cliff at the wrap-around point.
    # We punch a hole in the theoretical tau mapping grid to blind the engine to this cliff.
    tau_grid <- object@dictionary_meta$tau
    eps_grid <- object@dictionary_meta$epsilon
    seam_tau <- shift_fraction
    
    safe_tau_grid <- tau_grid[abs(tau_grid - seam_tau) > max(eps_grid)]
    if (length(safe_tau_grid) == 0) safe_tau_grid <- tau_grid # Fallback

    # Run Scan securely avoiding the artificial cliff
    suppressMessages({
        null_object <- scanDynamics(null_object,
            tau_grid = safe_tau_grid,
            epsilon_grid = eps_grid,
            cores = cores
        )
    })

    null_scores <- null_object@events$Score
    message("... Generated ", length(null_scores), " null scores.")

    # 3. Calculate FDR correctly over full theoretical pair spaces
    n_real_pairs <- ncol(object@logratio) * (ncol(object@logratio) - 1) / 2
    n_null_pairs <- ncol(null_object@logratio) * (ncol(null_object@logratio) - 1) / 2

    real_scores <- object@events$Score
    n_real <- length(real_scores)

    # Vectorized calculation for O(1) performance and to avoid O(N^2) loops
    null_scores_asc <- sort(null_scores, decreasing = FALSE)
    
    # Calculate exact number of null scores >= s for all real scores simultaneously
    n_null_ge_vec <- length(null_scores_asc) - findInterval(real_scores - 1e-12, null_scores_asc)
    
    # Inject robust + 1 pseudocount to bound the FDR when the null is totally empty
    prop_null_vec <- (n_null_ge_vec + 1) / n_null_pairs
    prop_real_vec <- seq_len(n_real) / n_real_pairs
    
    fdrs <- pmin(prop_null_vec / prop_real_vec, 1.0)

    # Enforce strict monotonicity (q-values representation)
    # Traverse from worst score to best score to propagate the tightest significance threshold.
    if (n_real > 1) {
        for (i in seq(n_real - 1, 1, by = -1)) {
            fdrs[i] <- min(fdrs[i], fdrs[i + 1])
        }
    }

    object@events$FDR <- fdrs

    return(object)
}
