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
setGeneric("validateGLMM", function(object, formula_str = NULL) standardGeneric("validateGLMM"))

#' @rdname validateGLMM
setMethod("validateGLMM", "dyprop", function(object, formula_str = NULL) {
    if (nrow(object@events) == 0) {
        warning("No events to validate.")
        return(object)
    }

    # Select candidates (e.g., only Switch and Decoupling, ignore Homeostasis)
    candidates <- object@events[object@events$Event_Class %in% c("Switch", "Decoupling"), ]

    if (nrow(candidates) == 0) {
        message("No significant transition events found for validation.")
        return(object)
    }

    message(">>> Validating ", nrow(candidates), " candidates with GLMM...")

    # Data setup
    X_clr <- object@logratio
    meta <- object@design # Assuming this contains patient_id, condition, etc.

    # Check metadata
    if (is.null(meta) || nrow(meta) == 0) {
        warning("@design slot is empty. Cannot fit mixed models without patient metadata.")
        return(object)
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

    fits <- list()

    for (i in seq_len(nrow(candidates))) {
        geneA <- candidates$GeneA[i]
        geneB <- candidates$GeneB[i]

        # Calc Ratio Y
        y_vec <- X_clr[, geneA] - X_clr[, geneB]

        # Build DF
        df_fit <- data.frame(
            y = y_vec,
            pseudotime = pseudotime
        )
        df_fit <- cbind(df_fit, meta) # Add metadata columns

        # Fit Model
        tryCatch(
            {
                fit <- glmmTMB::glmmTMB(as.formula(formula_str),
                    data = df_fit,
                    family = glmmTMB::t_family(),
                    dispformula = ~pseudotime
                )

                # Store result (maybe just summary or p-values to save space)
                # For now, store summary
                fits[[paste(geneA, geneB, sep = "_")]] <- summary(fit)
            },
            error = function(e) {
                warning(paste("GLMM fit failed for", geneA, geneB, ":", e$message))
            }
        )
    }

    object@glmm_fits <- fits
    return(object)
})
