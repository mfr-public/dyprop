context("network")

test_that("Network Reconstruction functions work", {
    # 1. Mock Data
    set.seed(123)
    n_genes <- 10
    n_samples <- 50
    pseudotime <- seq(0, 1, length.out = n_samples)

    # Create synthetic switch
    # Gene 1 and 2 switch at tau=0.5
    # CLR matrix: random noise for others
    X_clr <- matrix(rnorm(n = n_samples * n_genes), nrow = n_samples, ncol = n_genes)
    colnames(X_clr) <- paste0("Gene", 1:n_genes)

    # Imprint correlation change
    tau <- 0.5
    epsilon <- 0.05

    # Regime A (t < tau): Correlated
    idx_A <- which(pseudotime < tau - 2 * epsilon)
    X_clr[idx_A, 1] <- rnorm(length(idx_A))
    X_clr[idx_A, 2] <- X_clr[idx_A, 1] + rnorm(length(idx_A), 0, 0.1) # High rho

    # Regime B (t > tau): Uncorrelated (Decoupling)
    idx_B <- which(pseudotime > tau + 2 * epsilon)
    X_clr[idx_B, 1] <- rnorm(length(idx_B))
    X_clr[idx_B, 2] <- rnorm(length(idx_B)) # Low rho

    # Mock object
    obj <- new("dyprop",
        logratio = as.data.frame(X_clr),
        counts = as.data.frame(matrix(100, nrow = n_samples, ncol = n_genes)), # Dummy counts
        pseudotime = pseudotime
    )

    # Mock Event
    obj@events <- data.frame(
        GeneA = "Gene1",
        GeneB = "Gene2",
        Tau = tau,
        Epsilon = epsilon,
        Event_Class = "Decoupling",
        stringsAsFactors = FALSE
    )

    # 2. Test get_regimes
    regimes <- dyprop::get_regimes(pseudotime, tau, epsilon)
    expect_true(length(regimes$pre) > 0)
    expect_true(length(regimes$post) > 0)
    # Check disjointness / boundary buffer
    boundary_width <- 2 * epsilon
    # Pre indices should satisfy t < tau - width
    expect_true(all(pseudotime[regimes$pre] < (tau - boundary_width)))

    # 3. Test calcRewiring
    res <- dyprop::calcRewiring(obj, "Gene1", "Gene2") # Should find tau automatically

    expect_is(res, "list")
    expect_true(!any(is.na(res$Delta_P)))

    # Check logic: Pre should have high rho for Gene1 vs Gene2 (which is index 2 if Gene1 is target)
    # Rho definition used in code: 1 - VLR / (varA+varB)
    # Pre: VLR small -> Rho high.
    # Post: VLR large -> Rho low/zero.
    # So Delta P = Post - Pre should be Negative for Gene1-Gene2 pair.

    idx_target <- which(colnames(X_clr) == "Gene2")
    rho_pre <- res$Rho_Pre[idx_target]
    rho_post <- res$Rho_Post[idx_target]

    # High correlation in Pre means high Rho.
    # Uncorrelated in Post means low Rho.
    # Expect Rho_Pre > Rho_Post
    expect_gt(rho_pre, rho_post)

    # 4. Test detectHubs
    hubs <- dyprop::detectHubs(obj, "Gene1", tau, epsilon)
    expect_is(hubs, "list")
    expect_true(hubs$Delta_K > 0)
})
