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

    # 3. Calculate FDR
    real_scores <- object@events$Score
    n_real <- length(real_scores)
    n_null <- length(null_scores)

    fdrs <- numeric(n_real)

    # For each real score, how many nulls are greater?
    # Optimization: Sort both and loop? Or just simple empirical calc.
    # FDR(s) = (N_null_ge / N_null) / (N_real_ge / N_real)

    # Since real_scores are already sorted descending in scanDynamics:
    # rank i implies i scores are >= real_scores[i]

    # Sort nulls descending
    null_scores <- sort(null_scores, decreasing = TRUE)

    for (i in seq_len(n_real)) {
        s <- real_scores[i]

        # Count nulls >= s
        # Since nulls sorted desc, find first index where null < s
        # Or typical approach: sum(nulls >= s)
        n_null_ge <- sum(null_scores >= s)

        prop_null <- n_null_ge / n_null
        prop_real <- i / n_real # i is count of real scores >= s (since sorted)

        fdr <- prop_null / prop_real
        fdrs[i] <- min(fdr, 1.0) # Cap at 1
    }

    # Enforce monotonicity? (Usually q-values are monotonic)
    # cummin? No, usually we do cummax from bottom up or similar.
    # Standard BH: sort p-values. Here we have scores.
    # Higher score -> Lower FDR.
    # So if you have a high score with high FDR, but a lower score has lower FDR, that's weird.
    # Usually we enforce that higher scores have lower (or equal) FDR.
    # We can use cummin() on the FDRs if we traverse from top score to bottom score? NO.
    # If score A > score B, FDR(A) should be <= FDR(B).
    # So as we go down the list (scores decrease), FDR should increase.
    # So we take cummax() ? No.
    # We want to ensure: FDR is increasing function of rank.
    # We can take cummax() of FDRs as we go down the list?
    # Yes: If a better score had worse FDR, we'd be confused.
    # Actually, standard approach (Storey) produces q-values which are monotonic.
    # Let's apply cummax to ensure monotonicity.

    # fdrs <- cummax(fdrs) # Wait, if top hit has FDR 0.05, next hit has FDR 0.01, that's impossible mathematically if prop_real increases faster.
    # But due to noise it might happen.
    # Let's leave raw empirical FDR for now, or just cummax?
    # q-value = min(FDR(t)) for all t <= s.
    # So for a given score s, the q-value is the minimum FDR found for any score <= s? No.
    # q-value is min FDR for any threshold t <= s (more stringent).
    # Here more stringent means HIGHER score.
    # So q[i] = min(FDR[1:i]).
    # Yes, valid FDR is monotonic.
    # Since we are iterating from Highest Score (Stringent) to Lowest,
    # The FDR usually rises.
    # If we find a lower FDR later, it effectively improves the previous ones? No.
    # The q-value of a feature is the expected proportion of false positives incurred
    # when calling that feature significant.
    # So we should validly use:
    # object@events$FDR <- cummin(fdrs(reversed)) ... standard q-value logic is tricky.
    # Let's stick to the raw ratio definition from the paper docs for now:
    # "Empirical p-values... calculated relative to this null distribution."

    object@events$FDR <- fdrs

    return(object)
}
