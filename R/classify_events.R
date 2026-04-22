#' @include AllClasses.R
NULL

#' Classify Dynamic Events
#'
#' Applies a deterministic, multi-evidence decision tree to logically map the
#' geometrical outputs of Phase 1 (R2_Linear, R2_Sigmoid, Var_Delta) into rigid
#' biological topologies (Phase_Transition, Smooth_Drift, Decoupling).
#'
#' @param object A \code{dyprop} object with populated \code{@events} slot.
#'
#' @return A \code{dyprop} object with updated \code{@events$Event_Class} and 
#' \code{@events$Confidence} columns.
#'
#' @export
setGeneric("classifyEvents", function(object) standardGeneric("classifyEvents"))

#' @rdname classifyEvents
setMethod("classifyEvents", "dyprop", function(object) {
    if (nrow(object@events) == 0) {
        warning("No events found. Please execute scanDynamics() first.")
        return(object)
    }

    message(">>> Phase 2: Deterministic Topological Classification...")

    events <- object@events

    # Extract Empirical Null Bound arrays safely
    if (length(object@fdr_cutoff) == 0) {
        warning("No empirical FDR cutoffs found. Run estimateFDR() securely first. Utilizing extremely conservative default bounds.", call. = FALSE)
        r2_cut <- 0.5
        var_cut <- 5.0
        lin_cut <- 0.5
    } else {
        r2_cut <- object@fdr_cutoff$R2_Sigmoid
        var_cut <- object@fdr_cutoff$Var_Delta
        lin_cut <- object@fdr_cutoff$R2_Linear
    }

    # Execute deterministic categorization using the Topological Vector Paradigm natively scaled to FDR
    classifications <- dplyr::case_when(
        events$Var_Delta >= var_cut ~ "Decoupling",
        events$R2_Sigmoid >= r2_cut & events$R2_Sigmoid > (1.1 * events$R2_Linear) & events$Epsilon_Grid < 0.15 ~ "Phase_Transition",
        events$R2_Linear >= lin_cut | (events$R2_Sigmoid >= r2_cut & events$Epsilon_Grid >= 0.15) ~ "Smooth_Drift",
        TRUE ~ "Homeostasis"
    )

    # Compute universal 0-1 Confidence score mapping
    conf_scores <- dplyr::case_when(
        classifications == "Decoupling" ~ 1.0 - (1.0 / pmax(events$Var_Delta, 1.0)),
        classifications == "Phase_Transition" ~ events$R2_Sigmoid,
        classifications == "Smooth_Drift" ~ events$R2_Linear,
        TRUE ~ 0.0
    )

    events$Event_Class <- classifications
    events$Confidence <- conf_scores
    
    # Filter pure homeostasis (noise) natively here 
    events <- events[events$Event_Class != "Homeostasis", ]

    # Sort strictly by Confidence so Top Tier mathematically rigid hits are front-loaded
    events <- events[order(events$Confidence, decreasing = TRUE), ]

    object@events <- events

    message(">>> Classification Complete.")
    print(table(object@events$Event_Class))

    return(object)
})
