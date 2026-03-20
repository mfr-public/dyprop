#' @include AllClasses.R
NULL

#' Estimate False Discovery Rate (FDR)
#'
#' Estimates the statistical significance of detected events using an empirical Null Model.
#' The method permutes the pseudotime vector to destroy biological signal, generates
#' a distribution of 'null' scan scores, and calculates the False Discovery Rate (FDR)
#' for the real events.
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

    # 1. Create Null Object
    # We permute pseudotime exactly once (as per docs "We permute... once")
    # This destroys the temporal structure for all genes simultaneously.
    # To get stable statistics, scanning a random subset of 10,000 pairs
    # from this single permutation is usually sufficient if N_pairs is large.

    null_object <- object
    set.seed(42) # Reproducibility
    null_object@pseudotime <- sample(object@pseudotime)

    # 2. Run Scan on Null Object
    # ISSUE: scanDynamics scans ALL pairs. For 10,000 genes, that's 50M pairs.
    # Doing that just for the null is expensive if we only need 10,000 null scores.
    # Ideally, scanDynamics would accept a 'subset' or 'n_pairs' argument.
    # Current scanDynamics implementation scans *everything* in @logratio.

    # Optimization: Subset the inputs to scanDynamics to limit computation.
    # If we only need ~10,000 null scores, we can pick a subset of genes
    # such that n_genes choose 2 ~= 10,000.
    # n*(n-1)/2 = 10000 => n^2 ~= 20000 => n ~= 141 genes.

    n_features <- ncol(object@logratio)

    if (n_features * (n_features - 1) / 2 > n_permutations) {
        # Subset features
        n_subset <- ceiling(sqrt(2 * n_permutations))
        # Ensure we don't exceed actual features
        n_subset <- min(n_subset, n_features)

        # Pick random genes
        keep_idx <- sample(seq_len(n_features), n_subset)

        # Update null object to have only these genes
        null_object@logratio <- object@logratio[, keep_idx, drop = FALSE]
        # Also counts just in case (though scanDynamics uses logratio)
        if (nrow(object@counts) > 0) null_object@counts <- object@counts[, keep_idx, drop = FALSE]

        message("... Subsetting to ", n_subset, " random genes to generate approx ", n_permutations, " null scores.")
    }

    # Run Scan
    # We suppress messages to avoid noise
    suppressMessages({
        null_object <- scanDynamics(null_object,
            tau_grid = object@dictionary_meta$tau,
            epsilon_grid = object@dictionary_meta$epsilon,
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

    fdrs <- numeric(n_real)

    null_scores <- sort(null_scores, decreasing = TRUE)

    for (i in seq_len(n_real)) {
        s <- real_scores[i]

        n_null_ge <- sum(null_scores >= s)

        prop_null <- n_null_ge / n_null_pairs
        prop_real <- i / n_real_pairs

        fdr <- prop_null / prop_real
        fdrs[i] <- min(fdr, 1.0) # Cap at 1
    }

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
