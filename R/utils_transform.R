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
prepareComposition <- function(counts, pseudotime, design = matrix(nrow = 0, ncol = 0), k = 5, impute_zeros = TRUE) {
    message(">>> Phase 0: Pre-processing Data...")
    message("DEBUG: Input counts dim: ", paste(dim(counts), collapse = "x"))
    message("DEBUG: Input counts range: ", paste(range(as.matrix(counts)), collapse = " - "))


    # 1. kNN Pooling (Technical Dropout Correction)
    # Objective: Smooth technical zeros using manifold proximity (GLM-PCA)
    # k = 5 nearest neighbors

    counts_pooled <- counts # Default fallback

    if (k > 0 && requireNamespace("glmpca", quietly = TRUE) && requireNamespace("FNN", quietly = TRUE)) {
        message(">>> Step A: kNN Smoothing via Randomized SVD (irlba)...")

        # A. Dimension Reduction (Randomized SVD via irlba)
        # Scalable Log-normalization matching Seurat standards for highly sparse matrices
        L_dim <- min(ncol(counts), 50) # Latent dimensions
        counts_pooled <- tryCatch(
            {
                if (!requireNamespace("irlba", quietly = TRUE)) {
                    stop("Package 'irlba' is required for Fast SVD kNN but not found.")
                }
                # Log-normalization
                norm_counts <- log1p(as.matrix(counts))

                # Fast SVD
                svd_res <- irlba::irlba(norm_counts, nv = L_dim)
                X_latent <- svd_res$u %*% diag(svd_res$d)

                # B. Find Nearest Neighbors
                knn_res <- FNN::get.knn(X_latent, k = k)
                knn_idx <- knn_res$nn.index

                # C. Pool Counts (Sparse Matrix Multiplication)
                if (requireNamespace("Matrix", quietly = TRUE)) {
                    n_samples <- nrow(counts)
                    i_idx <- rep(seq_len(n_samples), times = k + 1)
                    j_idx <- c(seq_len(n_samples), as.vector(knn_idx))

                    W <- Matrix::sparseMatrix(
                        i = i_idx,
                        j = j_idx,
                        x = 1 / (k + 1),
                        dims = c(n_samples, n_samples)
                    )
                    counts_mat <- as.matrix(counts)
                    cp_temp <- as.matrix(W %*% counts_mat)
                    rownames(cp_temp) <- rownames(counts)
                    colnames(cp_temp) <- colnames(counts)
                } else {
                    # Fallback
                    counts_mat <- as.matrix(counts)
                    cp_temp <- matrix(0, nrow = nrow(counts), ncol = ncol(counts))
                    rownames(cp_temp) <- rownames(counts)
                    colnames(cp_temp) <- colnames(counts)
                    for (i in seq_len(nrow(counts))) {
                        idx <- c(i, knn_idx[i, ])
                        cp_temp[i, ] <- colMeans(counts_mat[idx, , drop = FALSE])
                    }
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
            if (!requireNamespace("irlba", quietly = TRUE)) message("Package 'irlba' not found. Skipping kNN smoothing.")
            if (!requireNamespace("FNN", quietly = TRUE)) message("Package 'FNN' not found. Skipping kNN smoothing.")
        } else {
            message("... Skipping kNN pooling (k=0).")
        }
        counts_pooled <- counts
    }

    # 2. Imputation (Biological Zeros)
    # Check if zCompositions is available
    if (impute_zeros && k > 0 && requireNamespace("zCompositions", quietly = TRUE)) {
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

    # 3. CLR Transform (Ecosystem Implementation)
    message("... Calculating CLR transform via propr::propr")

    if (any(counts_imp <= 0)) {
        warning("Zeros/Negatives found in CLR input. Adding small constant 1e-6.")
        counts_imp[counts_imp <= 0] <- 1e-6
    }

    # Use propr framework to handle CoDa structure securely
    # Note: winsorization of the difference trajectories happens natively in C++ layer
    pr <- propr::propr(counts_imp, metric = "rho", p = 0)
    X_clr <- pr@logratio

    # 4. Create dyprop object
    new("dyprop",
        counts = as.data.frame(counts_imp),
        logratio = as.data.frame(X_clr),
        matrix = pr@matrix,
        pairs = pr@pairs,
        pseudotime = as.numeric(pseudotime),
        design = as.data.frame(design),
        metric_cache = new.env(parent = emptyenv()),
        events = data.frame()
    )
}
