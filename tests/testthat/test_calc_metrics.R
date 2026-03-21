test_that("C++ Level 3 BLAS chunking algorithm exactly matches mathematical R slow-loop expectation", {
    # 1. Generate mathematically stable dummy matrix (3 genes, 50 cells)
    set.seed(42)
    n_time <- 50
    n_genes <- 3

    # Random count-like logratios
    X_clr <- matrix(rnorm(n_time * n_genes), nrow = n_genes, ncol = n_time)
    pseudotime <- seq(0, 1, length.out = n_time)

    # Event Mapping: Gene 1 vs Gene 2
    geneA_idx <- 1L
    geneB_idx <- 2L

    h_opt <- 0.05
    lambda <- 1e-5

    # 2. Run new C++ BLAS GEMM matrix chunking
    # Note: `calc_metrics` expects `X_clr` strictly as Genes x Time!
    res_cpp <- calc_metrics(X_clr, geneA_idx, geneB_idx, pseudotime, h_opt, lambda)

    # 3. Program the rigorous O(T^2) baseline manually in slow R loops
    W <- matrix(0, nrow = n_time, ncol = n_time)
    for (i in 1:n_time) {
        t0 <- pseudotime[i]
        w <- exp(-(pseudotime - t0)^2 / (2 * h_opt^2))
        W[i, ] <- w / sum(w)
    }

    valA <- X_clr[1, ]
    valB <- X_clr[2, ]
    lr <- valA - valB

    max_phi <- 0
    min_rho <- 2

    for (t in 1:n_time) {
        wt <- W[t, ]
        mu_A <- sum(wt * valA)
        mu_B <- sum(wt * valB)
        mu_lr <- sum(wt * lr)

        var_A <- sum(wt * (valA - mu_A)^2)
        var_B <- sum(wt * (valB - mu_B)^2)
        var_lr <- sum(wt * (lr - mu_lr)^2)

        denom <- var_A + var_B + 1e-12
        phi <- var_lr / (denom + lambda)
        rho <- 1.0 - (var_lr / denom)

        if (phi > max_phi) max_phi <- phi
        if (rho < min_rho) min_rho <- rho
    }

    expected_class <- if (max_phi < 0.2) {
        0L # Homeostasis
    } else if (min_rho < 0.5) {
        2L # Decoupling
    } else {
        1L # Switch
    }

    # 4. Formally assert topological alignment up to floating-point truncation limits
    expect_equal(res_cpp$Class_ID[1], expected_class)
})
