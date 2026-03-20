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

        df_fit <- data.frame(
            y = y_vec,
            pseudotime = pseudotime
        )
        df_fit <- cbind(df_fit, meta)

        res <- tryCatch(
            {
                fit <- glmmTMB::glmmTMB(as.formula(formula_str),
                    data = df_fit,
                    family = glmmTMB::t_family(),
                    dispformula = ~pseudotime
                )
                summary(fit)
            },
            error = function(e) {
                warning(paste("GLMM fit failed for", geneA, geneB, ":", e$message))
                NULL
            }
        )
        return(res)
    }

    # Parallel map
    fits_list <- parallel::mclapply(pairs_list, fit_worker, mc.cores = cores)

    # Associate names
    names(fits_list) <- sapply(pairs_list, function(p) paste(p$geneA, p$geneB, sep = "_"))

    # Pre-clean NULL errors
    fits_list <- fits_list[!vapply(fits_list, is.null, logical(1))]

    object@glmm_fits <- fits_list
    return(object)
})
