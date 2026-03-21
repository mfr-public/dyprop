#' @include AllClasses.R
NULL

#' Hierarchical Validation using GLMM
#'
#' Fits a Generalized Linear Mixed Model (GLMM) to significant gene pairs to
#' distinguish between true phase transitions and random condition/patient effects.
#' Uses glmmTMB to model heteroscedasticity and heavy-tailed residuals (t-distribution).
#'
#' @param object A \code{dyprop} object with classified events.
#' @param formula_str Character. The GLMM formula. Default:
#'  "y ~ splines::ns(pseudotime, df=3) + (1|patient_id)"
#' @param family The glmmTMB family. Default: \code{glmmTMB::t_family()}.
#'
#' @return A \code{dyprop} object with populated \code{@glmm_fits} slot.
#'
#' @importFrom glmmTMB glmmTMB t_family
#' @importFrom splines ns
#' @export
setGeneric("validateGLMM", function(object, formula_str = NULL, max_candidates = 10000, fdr_cutoff = 0.05) standardGeneric("validateGLMM"))

#' @rdname validateGLMM
setMethod("validateGLMM", "dyprop", function(object, formula_str = NULL, max_candidates = 10000, fdr_cutoff = 0.05) {
    if (nrow(object@events) == 0) {
        warning("No events to validate.")
        return(object)
    }

    # 1. Option B: Aggressive Statistical Gating via FDR
    object@events <- object@events[!is.na(object@events$FDR) & object@events$FDR <= fdr_cutoff, ]

    # Select candidates
    candidates <- object@events[object@events$Event_Class %in% c("Switch", "Decoupling"), ]

    if (nrow(candidates) == 0) {
        message("No significant transition events found for validation.")
        return(object)
    }

    # 2. Option A: Dynamic Top-K Resource Cap
    if (!is.null(max_candidates) && nrow(candidates) > max_candidates) {
        message(sprintf(">>> Subsetting %d candidates to strictly Top %d by scan Score to structurally preserve Memory bandwidth...", nrow(candidates), max_candidates))
        candidates <- candidates[order(candidates$Score, decreasing = TRUE), ]
        candidates <- candidates[1:max_candidates, ]
    }

    if (nrow(candidates) == 0) {
        message("No significant transition events found for validation.")
        return(object)
    }

    # Data setup
    X_clr <- object@logratio
    meta <- object@design # Assuming this contains patient_id, condition, etc.

    # Check metadata
    if (is.null(meta) || nrow(meta) == 0) {
        warning("@design slot is empty. Will fit generalized linear models without mixed/patient effects.")
        meta <- data.frame(dummy = rep(1, nrow(X_clr))) # Dummy structural metadata
    }

    pseudotime <- object@pseudotime

    if (is.null(formula_str)) {
        # Default formula
        # We construct a dataframe for fitting: y, pseudotime, patient_id
        # Check for patient_id col
        if (!"patient_id" %in% colnames(meta)) {
            # Fallback or error
            warning("'patient_id' not found in @design. Using simple LM.")
            formula_str <- "y ~ splines::ns(pseudotime, df=3)"
        } else {
            formula_str <- "y ~ splines::ns(pseudotime, df=3) + (1|patient_id)"
        }
    }

    cores <- getOption("mc.cores", parallel::detectCores() - 1L)
    if (cores < 1) cores <- 1L
    message(sprintf(">>> Validating %d candidates with GLMM on %d cores...", nrow(candidates), cores))

    # Identify pairs
    pairs_list <- lapply(seq_len(nrow(candidates)), function(i) {
        list(geneA = candidates$GeneA[i], geneB = candidates$GeneB[i])
    })

    fit_worker <- function(pair) {
        geneA <- pair$geneA
        geneB <- pair$geneB

        y_vec <- X_clr[, geneA] - X_clr[, geneB]

        df_fit <- data.frame(y = y_vec, pseudotime = pseudotime)
        df_fit <- cbind(df_fit, meta)

        res <- tryCatch(
            {
                fit <- glmmTMB::glmmTMB(as.formula(formula_str),
                    data = df_fit, family = glmmTMB::t_family(), dispformula = ~pseudotime
                )

                fit_summary <- summary(fit)
                cond_table <- fit_summary$coefficients$cond
                spline_rows <- grep("ns\\(pseudotime", rownames(cond_table))

                if (length(spline_rows) > 0) {
                    p_val <- min(cond_table[spline_rows, 4])
                    max_beta <- max(abs(cond_table[spline_rows, 1]))
                    c(p_val, max_beta)
                } else {
                    c(NA_real_, NA_real_)
                }
            },
            error = function(e) {
                # Silently catch non-convergent mathematical models into NA to protect IPC bandwidth
                c(NA_real_, NA_real_)
            }
        )
        return(res)
    }

    # Parallel map engineered to exclusively return primitive NUMERICS, structurally wrapped in rigid File-Descriptor Chunk bounds
    chunk_size <- 2000
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
        gc() # Safe pointer collection across child OS nodes natively
    }

    # Structurally unpack the float matrix [P_Value, Max_Beta] without caching massive C++ models
    p_vals_mat <- do.call(rbind, p_vals_list)

    # Initialize the target variable array
    object@events$P_Value <- rep(NA_real_, nrow(object@events))
    object@events$Max_Beta <- rep(NA_real_, nrow(object@events))

    for (i in seq_len(nrow(candidates))) {
        geneA <- candidates$GeneA[i]
        geneB <- candidates$GeneB[i]
        idx <- which(object@events$GeneA == geneA & object@events$GeneB == geneB)
        if (length(idx) > 0) {
            object@events$P_Value[idx[1]] <- p_vals_mat[i, 1]
            object@events$Max_Beta[idx[1]] <- p_vals_mat[i, 2]
        }
    }

    return(object)
})
