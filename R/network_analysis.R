#' @include AllClasses.R
NULL

#' Define Dynamic Regimes
#'
#' Partitions the pseudotime trajectory into "Pre-transition" and "Post-transition"
#' stable regimes based on the detected tipping point (\eqn{\tau}) and sharpness (\eqn{\epsilon}).
#'
#' @param pseudotime Numeric vector of pseudotime coordinates.
#' @param tau Numeric. The estimated tipping point.
#' @param epsilon Numeric. The estimated inverse sharpness (transition width).
#'
#' @return A list with two components:
#' \item{pre}{Integer vector of indices for the Pre-transition regime.}
#' \item{post}{Integer vector of indices for the Post-transition regime.}
#'
#' @export
get_regimes <- function(pseudotime, tau, epsilon) {
    boundary_width <- 2 * epsilon
    idx_pre <- which(pseudotime < (tau - boundary_width))
    idx_post <- which(pseudotime > (tau + boundary_width))
    return(list(pre = idx_pre, post = idx_post))
}

#' Calculate Proportionality Rewiring
#'
#' Computes the change in the static proportionality matrix (\eqn{\rho}) between
#' the Pre- and Post-transition regimes (\eqn{\Delta \mathbf{P}}).
#'
#' @param object A \code{dyprop} object.
#' @param geneA Character. Name of the first gene (or switch anchor).
#' @param geneB Character. Name of the second gene (or switch partner).
#' @param tau Numeric. Tipping point (optional; uses object@events if available).
#' @param epsilon Numeric. Sharpness (optional).
#'
#' @return A list containing:
#' \item{Rho_Pre}{The proportionality matrix for the Pre regime.}
#' \item{Rho_Post}{The proportionality matrix for the Post regime.}
#' \item{Delta_P}{The rewiring matrix (\eqn{Rho_{Post} - Rho_{Pre}}).}
#'
#' @export
calcRewiring <- function(object, geneA, geneB, tau = NULL, epsilon = NULL) {
    # If tau/epsilon not provided, try to find in events
    if (is.null(tau) || is.null(epsilon)) {
        ev <- object@events
        idx <- which(ev$GeneA == geneA & ev$GeneB == geneB)
        if (length(idx) == 0) {
            # Try swap
            idx <- which(ev$GeneA == geneB & ev$GeneB == geneA)
        }
        if (length(idx) > 0) {
            tau <- ev$Tau[idx[1]]
            epsilon <- ev$Epsilon[idx[1]]
        } else {
            stop("Event for given genes not found in @events. Please provide tau and epsilon manually.")
        }
    }

    regimes <- get_regimes(object@pseudotime, tau, epsilon)

    # We need counts or clr to recalculate Rho
    # object@logratio has the CLR data
    X_clr <- object@logratio

    if (length(regimes$pre) < 3 || length(regimes$post) < 3) {
        warning("Regimes contain too few samples for reliable Rho calculation.")
        return(list(Rho_Pre = NA, Rho_Post = NA, Delta_P = NA))
    }

    # Calculate Rho on subsets
    # propr::propr expects counts usually, but we have CLR.
    # Note: propr(counts) does CLR internally.
    # If we have CLR, we can calculate Rho directly as:
    # Rho = 1 - var(a - b) / (var(a) + var(b))
    # Or use weighted variance if needed, but here we assume static within regime.

    # Let's use the subset of CLR to calculate Rho
    # We need a helper for Rho from CLR
    # Calculate Rho on subsets using CLR data directly
    X_pre <- X_clr[regimes$pre, , drop = FALSE]
    X_post <- X_clr[regimes$post, , drop = FALSE]

    # Implementation: Calculate for the specific gene against ALL others to find hubs.

    # Vectorized Rho for specific gene vs ALL
    calc_rho_vec <- function(target_gene_idx, X, lambda = 1e-5) {
        tgt <- X[, target_gene_idx]
        var_tgt <- var(tgt)

        # Cols vars
        vars_all <- apply(X, 2, var)

        # VLRs
        # var(A - B) = var(A) + var(B) - 2cov(A,B)
        covs <- cov(X, tgt) # Vector
        vlrs <- vars_all + var_tgt - 2 * covs

        rhos <- 1 - vlrs / (vars_all + var_tgt + lambda)
        return(rhos)
    }

    idxA <- which(colnames(X_clr) == geneA)
    # idxB <- which(colnames(X_clr) == geneB) # Unused for now as we calculate A vs All

    rho_pre_A <- calc_rho_vec(idxA, X_pre)
    rho_post_A <- calc_rho_vec(idxA, X_post)

    delta_p_A <- rho_post_A - rho_pre_A

    # Store results?
    # Return Delta P vector for Gene A (as a Hub candidate)
    return(list(
        Gene = geneA,
        Rho_Pre = rho_pre_A,
        Rho_Post = rho_post_A,
        Delta_P = delta_p_A
    ))
}

#' Detect Haywire Hubs
#'
#' Identifies master regulators by calculating Differential Degree Centrality (\eqn{\Delta k}).
#'
#' @param delta_p_list List of Delta P vectors (output from calcRewiring).
#' @return Data frame of genes ranked by \eqn{\sum |\Delta P|}.
#' @export
detectHubs <- function(object, gene_of_interest, tau, epsilon) {
    # Wrapper to run calcRewiring and aggregate
    res <- calcRewiring(object, gene_of_interest, gene_of_interest, tau, epsilon) # geneB ignored if we look at A vs All

    delta_p <- res$Delta_P
    delta_k <- sum(abs(delta_p), na.rm = TRUE)

    # We might want to return the top rewiring partners
    top_partners <- head(sort(abs(delta_p), decreasing = TRUE), 10)

    return(list(
        Gene = gene_of_interest,
        Delta_K = delta_k,
        Top_Partners = top_partners
    ))
}
