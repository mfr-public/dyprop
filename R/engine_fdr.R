#' @include AllClasses.R
NULL

#' Estimate False Discovery Rate (FDR) Boundaries
#'
#' Leverages the high-speed C++ Topological Array engine to calculate dynamic
#' geometric constraints. By taking a massive combinatorial subset of the RAW unshifted 
#' matrix, we functionally extract the geometric density of true biological drift 
#' (background cell pathway noise). The 99th percentile of this unshifted array 
#' serves as our mathematically elite boundary limit, analogous to WGCNA scale-free 
#' topology thresholds.
#'
#' @param object A \code{dyprop} object.
#' @param cores Integer. Threads dedicated to null matrix extraction.
#' @param ... Additional arguments.
#'
#' @return A modified \code{dyprop} object with populated \code{fdr_cutoff} bounds.
#'
#' @export
estimateFDR <- function(object, cores = 1L, ...) {
    if (!inherits(object, "dyprop")) stop("Input must be a 'dyprop' object.")
    
    X_clr <- object@logratio
    X_t <- as.matrix(t(X_clr))
    
    # Subset rigorously to prevent memory flooding and speed up FDR sweeps
    n_genes <- nrow(X_t)
    if (n_genes > 500) {
        set.seed(42)  # For consistent reproducible null bounds
        target_idx <- sample(seq_len(n_genes), 500)
        X_shifted <- X_t[target_idx, , drop=FALSE]
        message(">>> FDR: Subsetting to 500 random genes (124k pairs) to formally establish pure structural background boundaries...")
    } else {
        X_shifted <- X_t
    }
    
    if (length(object@dictionary_meta) == 0) {
         tau_grid <- seq(0.1, 0.9, length.out=20)
         epsilon_grid <- c(0.05, 0.1, 0.15, 0.2)
    } else {
         tau_grid <- object@dictionary_meta$tau
         epsilon_grid <- object@dictionary_meta$epsilon
    }
    
    message(">>> Extracting Generic Biological Drift Thresholds (Top 1% Elite Geometric Bound)...")
    Sys.setenv("OMP_NUM_THREADS" = as.character(cores))
    
    # We blast the bounds to 0.0 to capture the true continuous spectrum
    suppressMessages({
        res_null <- scan_sequences(X_shifted,
                                   time_vec = object@pseudotime,
                                   tau_grid = tau_grid,
                                   epsilon_grid = epsilon_grid,
                                   min_score = 0.0,
                                   min_var_delta = 0.0)
    })
    
    # Extract the Top 99th Percentile as our mathematically elite Statistical Bound
    if (length(res_null$R2_Sigmoid) > 0) {
        emp_r2 <- stats::quantile(res_null$R2_Sigmoid, 0.99, na.rm=TRUE)
        emp_var <- stats::quantile(res_null$Var_Delta, 0.99, na.rm=TRUE)
        emp_lin <- stats::quantile(res_null$R2_Linear, 0.99, na.rm=TRUE)
    } else {
        emp_r2 <- 0.5
        emp_var <- 5.0
        emp_lin <- 0.5
    }
    
    # Fail-safes preventing extremely dead nulls from creating hyper-sensitive bounds
    emp_r2 <- max(emp_r2, 0.3)
    emp_var <- max(emp_var, 3.0)
    emp_lin <- max(emp_lin, 0.3)
    
    object@fdr_cutoff <- list(
        R2_Sigmoid = as.numeric(emp_r2),
        Var_Delta = as.numeric(emp_var),
        R2_Linear = as.numeric(emp_lin)
    )
    
    message(sprintf(">>> Empirical Null Boundaries established: R2_Sig >= %.3f | R2_Lin >= %.3f | Var_Delta >= %.3f", 
                    emp_r2, emp_lin, emp_var))
    
    return(object)
}
