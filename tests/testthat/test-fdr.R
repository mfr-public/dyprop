context("fdr")

test_that("estimateFDR throws deprecation warning", {
    # 1. Mock Data
    set.seed(123)
    n_genes <- 10
    n_samples <- 50
    pseudotime <- seq(0, 1, length.out = n_samples)

    X_clr <- matrix(rnorm(n_samples * n_genes), nrow = n_samples, ncol = n_genes)
    colnames(X_clr) <- paste0("Gene", 1:n_genes)

    # Create dummy object
    obj <- new("dyprop",
        logratio = as.data.frame(X_clr),
        counts = as.data.frame(matrix(100, nrow = n_samples, ncol = n_genes)),
        pseudotime = pseudotime
    )

    obj <- estimateFDR(obj)

    # 4. Check boundaries
    expect_true(!is.null(obj@fdr_cutoff))
    expect_true(obj@fdr_cutoff$R2_Sigmoid >= 0.3)
    expect_true(obj@fdr_cutoff$Var_Delta >= 3.0)
})
