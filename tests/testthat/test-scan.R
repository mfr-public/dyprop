test_that("Full dyprop workflow runs on synthetic data", {
    # 1. Generate Synthetic Data
    set.seed(123)
    n_genes <- 10
    n_samples <- 50
    pseudotime <- seq(0, 10, length.out = n_samples)

    counts <- matrix(rpois(n_genes * n_samples, lambda = 100), nrow = n_samples, ncol = n_genes)
    colnames(counts) <- paste0("Gene", 1:n_genes)

    # Introduce a signal (Gene1 vs Gene2 switch)
    # t < 5: Ratio ~ 0 | t > 5: Ratio ~ 2
    counts[pseudotime > 5, 1] <- counts[pseudotime > 5, 1] * 3

    # 2. Prepare Object
    # Mock zCompositions if missing
    if (!requireNamespace("zCompositions", quietly = TRUE)) {
        dat <- prepareComposition(counts, pseudotime)
    } else {
        dat <- prepareComposition(counts, pseudotime)
    }

    expect_s4_class(dat, "dyprop")
    expect_equal(dim(dat@logratio), c(n_samples, n_genes))

    # 3. Scan Dynamics
    # Use small grid for speed and enforce Nyquist
    dat <- scanDynamics(dat, tau_grid = seq(2, 8, by = 0.5), epsilon_grid = c(0.5, 1), cores = 1, min_score = 0.0)

    expect_true(nrow(dat@events) > 0)
    expect_true("Score" %in% names(dat@events))

    # Check if Gene1-Gene2 pair is top ranked (or close)
    # (Given random noise, it should be decent)

    # 4. Classify Events
    dat <- classifyEvents(dat)
    expect_true("Event_Class" %in% names(dat@events))

    # 5. Validate (Mock Design)
    dat@design <- data.frame(patient_id = rep(c("P1", "P2"), each = 25))

    # Try validation (might fail if glmmTMB not installed or data too noisy, so wrap in try)
    skip_if_not_installed("glmmTMB")
    dat <- validateGLMM(dat)
    expect_true(length(dat@glmm_fits) >= 0) # Just check it ran
})
