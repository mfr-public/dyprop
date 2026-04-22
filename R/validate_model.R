#' @include AllClasses.R
NULL

#' Hierarchical Validation using GLMM
#'
#' Fits a series of Generalized Linear Mixed Models (GLMM) via glmmTMB natively 
#' mapping to the Topological Vector Paradigm. Natively executes dynamic Likelihood
#' Ratio Tests routing splines, linear models, and variance break constraints.
#'
#' @param object A \code{dyprop} object with classified events.
#' @param max_candidates Integer. Top candidates to process conditionally by class. Default 1000.
#' @param p_value_cutoff Numeric. Final BH-adjusted significance threshold. Default 0.05.
#'
#' @return A \code{dyprop} object with populated \code{@events$P_Value_FDR}.
#'
#' @importFrom glmmTMB glmmTMB
#' @importFrom splines ns
#' @export
setGeneric("validateGLMM", function(object, max_candidates = 1000, p_value_cutoff = 0.05, cores = max(1L, parallel::detectCores() - 1L)) standardGeneric("validateGLMM"))

#' @rdname validateGLMM
setMethod("validateGLMM", "dyprop", function(object, max_candidates = 1000, p_value_cutoff = 0.05, cores = max(1L, parallel::detectCores() - 1L)) {
    if (nrow(object@events) == 0) {
        warning("No events to validate.")
        return(object)
    }

    # Extract Topological High-Confidence Candidates (Bypassing Smooth Drift for Spline destruction)
    candidates <- object@events[!is.na(object@events$Confidence) & object@events$Confidence >= 0.5 & object@events$Event_Class %in% c("Phase_Transition", "Decoupling"), ]

    if (nrow(candidates) == 0) {
        message("No high-confidence transition events found for validation.")
        return(object)
    }

    if (!is.null(max_candidates) && nrow(candidates) > max_candidates) {
        message(sprintf(">>> Subsetting %d candidates to strictly Top %d by Confidence to structurally preserve Memory bandwidth...", nrow(candidates), max_candidates))
        candidates <- candidates[order(candidates$Confidence, decreasing = TRUE), ]
        candidates <- candidates[1:max_candidates, ]
    }

    # Data setup
    X_clr <- object@logratio
    meta <- object@design 

    if (is.null(meta) || nrow(meta) == 0 || !"patient_id" %in% colnames(meta)) {
        warning("'patient_id' not found in @design. Falling back to dummy ID for structural matrix safety.")
        meta <- data.frame(patient_id = rep("P1", nrow(X_clr))) 
    }

    pseudotime <- object@pseudotime

    if (is.null(cores) || cores < 1) cores <- 1L
    message(sprintf(">>> Validating %d candidates with Dynamic GLMM routes on %d cores...", nrow(candidates), cores))

    pairs_list <- lapply(seq_len(nrow(candidates)), function(i) {
        list(
            geneA = candidates$GeneA[i], 
            geneB = candidates$GeneB[i],
            event_class = candidates$Event_Class[i],
            tau_event = candidates$Tau_Grid[i],
            tau_dec = candidates$Tau_Decouple[i]
        )
    })

    fit_worker <- function(pair) {
        geneA <- pair$geneA
        geneB <- pair$geneB
        event_class <- pair$event_class
        tau_event <- pair$tau_event
        tau_dec <- pair$tau_dec

        y_vec <- X_clr[, geneA] - X_clr[, geneB]
        df_fit <- data.frame(y = y_vec, pseudotime = pseudotime)
        df_fit <- cbind(df_fit, meta)

        res <- tryCatch(
            {
                if (event_class == "Phase_Transition") {
                    # Execute Targeted Structural Hinge (Chow Break) precisely at C++ tau parameter
                    df_fit$post_tau <- df_fit$pseudotime > tau_event
                    fit <- glmmTMB::glmmTMB(y ~ pseudotime * post_tau + (1|patient_id), data = df_fit, family = glmmTMB::t_family())
                    null_fit <- glmmTMB::glmmTMB(y ~ pseudotime + (1|patient_id), data = df_fit, family = glmmTMB::t_family())
                } else if (event_class == "Decoupling") {
                    # Execute Targeted Variance Fracture
                    df_fit$post_break <- df_fit$pseudotime > tau_dec
                    fit <- glmmTMB::glmmTMB(y ~ pseudotime + (1|patient_id), dispformula = ~ post_break, data = df_fit, family = glmmTMB::t_family())
                    null_fit <- glmmTMB::glmmTMB(y ~ pseudotime + (1|patient_id), dispformula = ~ 1, data = df_fit, family = glmmTMB::t_family())
                } else {
                    return(c(NA_real_))
                }

                lrt_res <- anova(fit, null_fit)
                p_val <- lrt_res$`Pr(>Chisq)`[2]
                c(p_val)
            },
            error = function(e) {
                c(NA_real_)
            }
        )
        return(res)
    }

    # Parallel Processing using Safe Blocks
    chunk_size <- 500
    n_pairs <- length(pairs_list)
    n_chunks <- ceiling(n_pairs / chunk_size)
    p_vals_list <- list()

    for (ch in seq_len(n_chunks)) {
        start_idx <- (ch - 1) * chunk_size + 1
        end_idx <- min(ch * chunk_size, n_pairs)
        chunk_pairs <- pairs_list[start_idx:end_idx]

        message(sprintf("... Processing GLMM Chunk %d of %d (Pairs %d to %d)...", ch, n_chunks, start_idx, end_idx))

        if (requireNamespace("pbmcapply", quietly = TRUE) && .Platform$OS.type != "windows") {
            chunk_pvals <- pbmcapply::pbmclapply(chunk_pairs, fit_worker, mc.cores = cores)
        } else {
            chunk_pvals <- parallel::mclapply(chunk_pairs, fit_worker, mc.cores = cores)
        }

        p_vals_list <- c(p_vals_list, chunk_pvals)
        gc()
    }

    p_vals_vec <- unlist(p_vals_list)

    object@events$P_Value <- rep(NA_real_, nrow(object@events))

    for (i in seq_len(nrow(candidates))) {
        geneA <- candidates$GeneA[i]
        geneB <- candidates$GeneB[i]
        idx <- which(object@events$GeneA == geneA & object@events$GeneB == geneB)
        if (length(idx) > 0) {
            object@events$P_Value[idx[1]] <- p_vals_vec[i]
        }
    }

    # Rigorous Multiple-Testing Correction (FDR) sequentially enforced on GLMM P-Values
    object@events$P_Value_FDR <- rep(NA_real_, nrow(object@events))
    
    valid_idx <- !is.na(object@events$P_Value)
    if (sum(valid_idx) > 0) {
        object@events$P_Value_FDR[valid_idx] <- p.adjust(object@events$P_Value[valid_idx], method = "BH")
    }

    # Filter to final mathematically secure structure (retaining Smooth_Drift which skips verification)
    keep_idx <- (valid_idx & object@events$P_Value_FDR < p_value_cutoff) | (object@events$Event_Class == "Smooth_Drift")
    object@events <- object@events[keep_idx, ]

    message(">>> Validation Complete. Remaining significant candidates: ", nrow(object@events))
    return(object)
})
