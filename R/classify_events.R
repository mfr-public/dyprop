#' @include AllClasses.R
NULL

#' Classify Dynamic Events
#'
#' Applies a multi-evidence decision tree to distinguish between stoichiometric switches
#' (Type 1) and network decoupling (Type 2). This method performs "Lazy Evaluation,"
#' calculating the computationally expensive Phi and Rho metrics only for the
#' significant candidates identified in the scanning phase.
#'
#' @param object A \code{dyprop} object with populated \code{@events} slot.
#' @param lambda Numeric. Regularization constant for Phi calculation. Default 1e-5.
#'
#' @return A \code{dyprop} object with updated \code{@events$Event_Class} and
#'  populated \code{@Phi_Array} and \code{@Rho_Array} (sparse/selective population).
#'
#' @export
setGeneric("classifyEvents", function(object, lambda = 1e-5) standardGeneric("classifyEvents"))

#' @rdname classifyEvents
setMethod("classifyEvents", "dyprop", function(object, lambda = 1e-5) {
    if (!"FDR" %in% colnames(object@events)) {
        warning("FDR column missing. Please execute estimateFDR() before classification.")
        return(object)
    }

    # Strict Baseline Constraint: Only compute expensive metrics for valid events
    initial_n <- nrow(object@events)
    object@events <- object@events[!is.na(object@events$FDR) & object@events$FDR < 0.05, ]
    events <- object@events
    n_events <- nrow(events)

    if (n_events == 0) {
        message(sprintf(">>> Safety Filter: No significant events retained (FDR < 0.05). Reduced %d to 0.", initial_n))
        return(object)
    }

    message(sprintf(">>> Classifying %d significant candidates (FDR < 0.05)...", n_events))

    # Extract Data
    X_clr <- object@logratio # Samples x Genes
    pseudotime <- object@pseudotime

    h_opt <- 0.05 * (max(pseudotime) - min(pseudotime))

    # Identify numerical indices corresponding to strings
    geneA_idx <- match(events$GeneA, colnames(X_clr))
    geneB_idx <- match(events$GeneB, colnames(X_clr))

    # Safely migrate data structure to C++
    X_t <- as.matrix(t(X_clr))

    results <- calc_metrics(X_t, geneA_idx, geneB_idx, pseudotime, h_opt, lambda)

    Phi_traces <- results$Phi
    Rho_traces <- results$Rho

    classes <- character(n_events)

    for (i in seq_len(n_events)) {
        phi_t <- Phi_traces[i, ]
        rho_t <- Rho_traces[i, ]

        max_phi <- max(phi_t)
        min_rho <- min(rho_t)

        TH_PHI <- 0.2
        TH_RHO <- 0.5

        if (max_phi < TH_PHI) {
            classes[i] <- "Homeostasis"
        } else {
            if (min_rho < TH_RHO) {
                classes[i] <- "Decoupling"
            } else {
                classes[i] <- "Switch"
            }
        }

        # Lazy Evaluation Storage: Write trace safely to environment cache
        pair_key <- paste(events$GeneA[i], events$GeneB[i], sep = "_")
        object@metric_cache[[pair_key]] <- list(Phi = phi_t, Rho = rho_t)
    }

    object@events$Event_Class <- classes

    message(">>> Classification Complete.")
    table(classes)

    return(object)
})
