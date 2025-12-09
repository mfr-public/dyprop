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
#' @useDynLib dyprop, .registration = TRUE
#' @export
scanDynamics <- function(object, 
                         tau_grid = seq(0.1, 0.9, length.out = 20),
                         epsilon_grid = c(0.05, 0.1, 0.15, 0.2),
                         cores = 1L) {
  
  # 1. Input Validation
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
  
  # Generate All Unique Pairs (Combinatorics)
  # WARNING: For 10,000 genes, this is 50 million pairs. 
  # We strongly recommend the user pre-filters the object using propr::subset
  idx <- combn(n_genes, 2)
  n_pairs <- ncol(idx)
  
  message(">>> Scanning ", n_pairs, " gene pairs across ", length(tau_grid) * length(epsilon_grid), " archetypes...")
  
  # Construct Y (Pairs x Time) - efficiently
  # We transpose X_clr to (Genes x Samples) for faster column access
  X_t <- t(X_clr)
  
  # 3. RUN C++ SCAN
  # We'll need to update the C++ signature to accept X_t (Genes x Time) and the index pairs
  # to avoid allocating the massive Y matrix in R RAM.
  # But given our current C++ 'scan_kernels.cpp' accepts Y, let's construct it.
  # Ideally, for huge data, we would push the "Pair Construction" to C++.
  
  # For now (v0.1), let's construct Y in chunks or assume the user has subsetted.
  Y <- X_clr[, idx[1, ]] - X_clr[, idx[2, ]] # (Samples x Pairs)
  Y <- t(Y) # (Pairs x Samples) - matches C++ expectation
  
  # 4. Call Engine
  # Ensure OpenMP uses the requested cores
  Sys.setenv("OMP_NUM_THREADS" = as.character(cores))
  
  results_list <- scan_sequences(Y, 
                                 time_vec = object@pseudotime, 
                                 tau_grid = tau_grid, 
                                 epsilon_grid = epsilon_grid)
  
  # 5. Format Results
  events_df <- data.frame(
    GeneA = gene_names[idx[1, ]],
    GeneB = gene_names[idx[2, ]],
    Tau_Grid = results_list$tau,
    Epsilon_Grid = results_list$epsilon,
    Score = results_list$score,
    Event_Class = NA, # To be filled by classifyEvents
    FDR = NA
  )
  
  # Sort by Score (Best matches first)
  events_df <- events_df[order(events_df$Score, decreasing = TRUE), ]
  
  # Update Object
  object@events <- events_df
  object@dictionary_meta <- list(tau = tau_grid, epsilon = epsilon_grid)
  
  message(">>> Scan Complete. Top match score: ", round(max(events_df$Score), 2))
  
  return(object)
}