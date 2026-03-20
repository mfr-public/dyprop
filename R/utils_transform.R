#' @include AllClasses.R
NULL

#' Prepare Compositional Data
#'
#' A wrapper function to perform the pre-processing pipeline:
#' 1. Technical Dropout Correction (kNN pooling)
#' 2. Biological Zero Imputation (Bayesian Replacement)
#' 3. CLR Transformation
#'
#' @param counts Matrix. Raw counts (Samples x Genes).
#' @param pseudotime Numeric vector.
#' @param design Matrix. Design matrix for validation (optional).
#' @param k Integer. Number of neighbors for kNN. Default 5.
#'
#' @return A \code{dyprop} object ready for scanning.
#'
#' @importFrom propr propr
#' @export
prepareComposition <- function(counts, pseudotime, design = matrix(nrow = 0, ncol = 0), k = 5) {
    message(">>> Phase 0: Pre-processing Data...")
    message("DEBUG: Input counts dim: ", paste(dim(counts), collapse = "x"))
    message("DEBUG: Input counts range: ", paste(range(as.matrix(counts)), collapse = " - "))


    # 1. kNN Pooling (Technical Dropout Correction)
    # Objective: Smooth technical zeros using manifold proximity (GLM-PCA)
    # k = 5 nearest neighbors

    counts_pooled <- counts # Default fallback

    if (k > 0 && requireNamespace("glmpca", quietly = TRUE) && requireNamespace("FNN", quietly = TRUE)) {
        message(">>> Step A: kNN Smoothing via GLM-PCA...")

        # A. Dimension Reduction (GLM-PCA)
        # We use Poisson family to handle counts directly without log-pseudocount bias
        L_dim <- min(ncol(counts), 50) # Latent dimensions
        counts_pooled <- tryCatch(
            {
                gpca <- glmpca::glmpca(counts, L = L_dim, fam = "poi")
                X_latent <- gpca$factors

                # B. Find Nearest Neighbors
                knn_res <- FNN::get.knn(X_latent, k = k)
                knn_idx <- knn_res$nn.index

                # C. Pool Counts (Raw Averaging)
                counts_mat <- as.matrix(counts)
                cp_temp <- matrix(0, nrow = nrow(counts), ncol = ncol(counts))
                rownames(cp_temp) <- rownames(counts)
                colnames(cp_temp) <- colnames(counts)

                for (i in seq_len(nrow(counts))) {
                    # Neighbors + Self
                    idx <- c(i, knn_idx[i, ])
                    # Average raw counts to preserve scale
                    cp_temp[i, ] <- colMeans(counts_mat[idx, , drop = FALSE])
                }
                message(sprintf("... pooled %d neighbors.", k))
                cp_temp
            },
            error = function(e) {
                warning("GLM-PCA failed: ", e$message, ". Skipping kNN pooling.")
                return(counts)
            }
        )
    } else {
        if (k > 0) {
            if (!requireNamespace("glmpca", quietly = TRUE)) message("Package 'glmpca' not found. Skipping kNN smoothing.")
            if (!requireNamespace("FNN", quietly = TRUE)) message("Package 'FNN' not found. Skipping kNN smoothing.")
        } else {
            message("... Skipping kNN pooling (k=0).")
        }
        counts_pooled <- counts
    }

    # 2. Imputation (Biological Zeros)
    # Check if zCompositions is available
    if (k > 0 && requireNamespace("zCompositions", quietly = TRUE)) {
        # Check for zeros
        if (any(counts_pooled == 0)) {
            message(">>> Step B: Imputing zeros with zCompositions::cmultRepl (GBM)...")
            # cmultRepl expects Samples x Genes
            counts_imp <- tryCatch(
                {
                    zCompositions::cmultRepl(counts_pooled, label = 0, method = "GBM", output = "p-counts", suppress.print = TRUE)
                },
                error = function(e) {
                    warning("GBM Imputation failed: ", e$message, ". Falling back to +1 pseudo-count.")
                    return(counts_pooled + 1)
                }
            )
        } else {
            counts_imp <- counts_pooled
        }
    } else {
        if (k > 0) {
            warning("zCompositions not installed. Skipping imputation (simple +1 pseudo-count used).")
        } else {
            message("... Skipping GBM Imputation (k=0). Using +1 pseudo-count.")
        }
        counts_imp <- counts_pooled + 1
    }

    # 3. CLR Transform (Manual Implementation)
    message("... Calculating CLR transform")
    message("DEBUG: counts_imp[1:5,1:5]:")
    print(counts_imp[1:5, 1:5])
    message("DEBUG: row sums: ", paste(head(rowSums(counts_imp)), collapse = " "))


    # helper for CLR
    clr_fun <- function(x) {
        log_x <- log(x)
        return(log_x - mean(log_x))
    }

    # Apply CLR row-wise (Samples)
    # counts_imp is Samples x Genes
    if (any(counts_imp <= 0)) {
        warning("Zeros/Negatives found in CLR input. Adding small constant.")
        counts_imp[counts_imp <= 0] <- 1e-6
    }
    X_clr <- t(apply(counts_imp, 1, clr_fun))

    # 4. Create dyprop object
    new("dyprop",
        counts = as.data.frame(counts),
        logratio = as.data.frame(X_clr),
        matrix = t(X_clr), # Store CLR (Genes x Samples) for scanning
        pairs = numeric(), # Empty pairs vector (numeric)
        pseudotime = as.numeric(pseudotime),
        design = as.data.frame(design),
        Phi_Array = array(0, dim = c(0, 0, 0)),
        Rho_Array = array(0, dim = c(0, 0, 0))
    )
}
