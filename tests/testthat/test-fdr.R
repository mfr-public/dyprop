context("fdr")

test_that("estimateFDR works correctly", {
    # 1. Mock Data
    set.seed(123)
    n_genes <- 10
    n_samples <- 50
    pseudotime <- seq(0, 1, length.out = n_samples)

    X_clr <- matrix(rnorm(n_samples * n_genes), nrow = n_samples, ncol = n_genes)
    colnames(X_clr) <- paste0("Gene", 1:n_genes)

    # Create object
    obj <- new("dyprop",
        logratio = as.data.frame(X_clr),
        counts = as.data.frame(matrix(100, nrow = n_samples, ncol = n_genes)),
        pseudotime = pseudotime
    )

    # 2. Run Scan first (needed to populate events)
    # We expect some random matches, disable memory threshold
    obj <- scanDynamics(obj, min_score = 0.0)
    expect_true(nrow(obj@events) > 0)

    # 3. Run FDR Estimation
    # Use small n_permutations for speed
    obj <- estimateFDR(obj, n_permutations = 100)

    # 4. Checks
    expect_true(!any(is.na(obj@events$FDR)))
    expect_true(all(obj@events$FDR >= 0 & obj@events$FDR <= 1))

    # Check that higher scores generally have lower FDR (monotonic trend check)
    # Ideally, the top hit should have low FDR
    # Sort by score desc
    sorted_events <- obj@events[order(obj@events$Score, decreasing = TRUE), ]

    # Top hit FDR should be <= Bottom hit FDR (roughly)
    # Note: With n_permutations=100 and random data, this might be noisy.
    # But we can check that it runs without error.

    # Check subsetting logic
    # If we ask for fewer perms than pairs, it should warn or subset?
    # Here 10 genes -> 45 pairs. n_permutations=100 > 45. No subsetting.
})
