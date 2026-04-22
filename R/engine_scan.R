#' @include AllClasses.R
NULL

#' Scan for Dynamic Proportionality Events
#'
#' Performs a vectorized template match of log-ratio trajectories against a dictionary
#' of sigmoid boundary functions. This method utilizes a high-performance C++ backend
#' (OpenMP threaded) to scan millions of pairs efficiently.
#'
#' @param object A \code{dyprop} object containing the CLR-transformed data.
#' @param tau_grid Numeric vector. The grid of tipping points (fraction of pseudotime) to scan.
#'  Defaults to seq(0.1, 0.9, length.out=20).
#' @param epsilon_grid Numeric vector. The grid of transition sharpness values.
#'  Defaults to c(0.05, 0.1, 0.15, 0.2).
#' @param cores Integer. Number of CPU cores to use for parallel processing.
#'  Defaults to 1 (or allow OpenMP to auto-detect).
#'
#' @return A \code{dyprop} object with the \code{events} slot populated.
#'  The \code{Phi_Array} and \code{Rho_Array} slots are left empty (Lazy Evaluation)
#'  until \code{classifyEvents} is called.
#'
#' @importFrom Rcpp evalCpp
#' @importFrom RcppParallel RcppParallelLibs
#' @useDynLib dyprop, .registration = TRUE
#' @export
scanDynamics <- function(object,
                         tau_grid = seq(0.1, 0.9, length.out = 20),
                         epsilon_grid = c(0.05, 0.1, 0.15, 0.2),
                         cores = 1L,
                         min_score = 0.5,
                         min_var_delta = 5.0,
                         method = c("pearson", "bicor")) {
  # 1. Input Validation
  method_idx <- match.arg(method)
  method_int <- ifelse(method_idx == "bicor", 1L, 0L)
  if (!inherits(object, "dyprop")) stop("Input must be a 'dyprop' object.")
  if (length(object@pseudotime) == 0) stop("Pseudotime slot is empty.")

  # 2. Extract Data
  # dyprop inherits from propr, so we use @logratio if available, or compute ratios
  # Note: The C++ engine expects the Ratio Matrix Y (Pairs x Time)
  # For efficiency, we scan *all* pairs in the propr object (which are usually pre-filtered)

  # However, propr stores a Matrix (Samples x Genes).
  # We need to construct the Ratio Matrix Y (Pairs x Samples)
  # This can be huge, so we process it carefully or rely on the propr @logratio
  # structure if it stores pairwise ratios (it usually stores CLR).

  # Strategy: We assume 'object' is a propr object where we want to scan specific pairs
  # OR we generate all pairs on the fly.
  # For the vector scan, we often calculate Y = clr[,A] - clr[,B]

  # To avoid memory explosion in R, we pass the CLR matrix to C++ and let it
  # generate ratios on the fly?
  # *Correction*: Our C++ engine 'scan_sequences' takes 'arma::mat Y'.
  # So we must compute Y in R first.

  message(">>> Preparing Ratio Matrix for ", nrow(object@counts), " samples...")

  # Get CLR data
  X_clr <- object@logratio
  n_genes <- ncol(X_clr)
  gene_names <- colnames(X_clr)
  if (is.null(gene_names) || length(gene_names) == 0) {
    gene_names <- paste0("Gene_", seq_len(n_genes))
  }

  n_pairs <- as.numeric(n_genes) * (n_genes - 1) / 2

  # Nyquist safety check
  min_eps <- min(epsilon_grid)
  step_tau <- if (length(tau_grid) > 1) tau_grid[2] - tau_grid[1] else 0
  if (step_tau > min_eps) {
    stop(sprintf("Nyquist Theorem violated: Grid step size tau (%f) must be <= minimum epsilon (%f).", step_tau, min_eps))
  }

  message(">>> Memory-Safe Scanning ", n_pairs, " gene pairs across ", length(tau_grid) * length(epsilon_grid) * 2, " topological archetypes...")

  # Construct X_t (Genes x Time) for C++
  X_t <- as.matrix(t(X_clr))

  # 4. Call Engine
  Sys.setenv("OMP_NUM_THREADS" = as.character(cores))

  results_list <- scan_sequences(X_t,
    time_vec = object@pseudotime,
    tau_grid = tau_grid,
    epsilon_grid = epsilon_grid,
    min_score = min_score,
    min_var_delta = min_var_delta,
    method = method_int
  )

  # 5. Format Results
  if (length(results_list$GeneA_Idx) == 0) {
    object@events <- data.frame(
      GeneA = character(0),
      GeneB = character(0),
      Tau_Grid = numeric(0),
      Epsilon_Grid = numeric(0),
      R2_Linear = numeric(0),
      R2_Sigmoid = numeric(0),
      Tau_Decouple = numeric(0),
      Var_Delta = numeric(0),
      Model_Type = character(0),
      stringsAsFactors = FALSE
    )
    return(object)
  }

  events_df <- data.frame(
    GeneA = gene_names[results_list$GeneA_Idx],
    GeneB = gene_names[results_list$GeneB_Idx],
    Tau_Grid = results_list$Tau,
    Epsilon_Grid = results_list$Epsilon,
    R2_Linear = results_list$R2_Linear,
    R2_Sigmoid = results_list$R2_Sigmoid,
    Tau_Decouple = results_list$Tau_Decouple,
    Var_Delta = results_list$Var_Delta,
    Model_Type = ifelse(results_list$Model_Type == 1, "Logistic", "Gaussian"),
    stringsAsFactors = FALSE
  )

  # Update Object
  object@events <- events_df
  object@dictionary_meta <- list(tau = tau_grid, epsilon = epsilon_grid)

  message(">>> Scan Complete. Topological candidates extracted: ", nrow(events_df))

  return(object)
}
