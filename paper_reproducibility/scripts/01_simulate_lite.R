# ==============================================================================
# Dynamic Proportionality - "Lite" Benchmark for Rapid Verification
# Topology: Rooted Bifurcation (Stem -> A / Stem -> B)
# Scale: 500 Cells, 1,000 Genes (Scaled down for CI/Testing)
# ==============================================================================

library(splatter)
library(scater)
library(SingleCellExperiment)
library(dplyr)

set.seed(2025)

# 1. CONFIGURATION (Scaled Down) ---------------------------------
N_CELLS <- 500
N_GENES <- 1000

# PAIR CONFIGURATION
N_PAIRS <- 20
N_IMPLANTS <- N_PAIRS * 3 # 60 Genes involved

# 2. GENERATE BIFURCATING BACKBONE ---------------------------------------------
message(">>> Generating Bifurcating Backbone (N=", N_CELLS, ")...")

params <- newSplatParams()
params <- setParams(params, list(
    batchCells = N_CELLS,
    nGenes = N_GENES,
    group.prob = c(0.2, 0.4, 0.4),
    path.from = c(0, 1, 1),
    path.nSteps = c(10, 20, 20), # Shorter paths
    lib.loc = 11,
    lib.scale = 0.5,
    dropout.type = "experiment",
    dropout.mid = 0,
    dropout.shape = -1,
    bcv.common = 0.2
))

sim <- splatSimulate(params, method = "path", verbose = T, sparsify = F)

# 3. EXTRACT LINEAGE A -----------------------------------
cell_subset <- colData(sim)$Group %in% c("Path1", "Path2")
sim_A <- sim[, cell_subset]

counts_A <- counts(sim_A)
meta_A <- colData(sim_A)

meta_A$GlobalTime <- NA
is_root <- meta_A$Group == "Path1"
is_branch <- meta_A$Group == "Path2"

steps_root <- meta_A$Step[is_root]
if (length(steps_root) > 0) meta_A$GlobalTime[is_root] <- (steps_root / max(steps_root)) * 0.33

steps_branch <- meta_A$Step[is_branch]
if (length(steps_branch) > 0) meta_A$GlobalTime[is_branch] <- 0.33 + ((steps_branch / max(steps_branch)) * 0.67)

pseudotime <- meta_A$GlobalTime
lib_factors <- colSums(counts_A) / mean(colSums(counts_A))

message(">>> Lineage A Extracted. N=", ncol(sim_A), " Cells.")

# 4. IMPLANTATION ------------------------------------------
message(">>> Implanting ", N_IMPLANTS, " Ground Truth Pairs...")

truth_table <- data.frame(
    PairID = paste0("Pair_", 1:N_IMPLANTS),
    GeneA = seq(1, by = 2, length.out = N_IMPLANTS),
    GeneB = seq(2, by = 2, length.out = N_IMPLANTS),
    Type = rep(c("Null", "Switch", "Decouple"), each = N_PAIRS),
    Tau = runif(N_IMPLANTS, 0.4, 0.8),
    Epsilon = runif(N_IMPLANTS, 0.03, 0.15),
    EffectSize = runif(N_IMPLANTS, 3, 8)
)

counts_implanted <- counts_A

logistic <- function(t, start, end, tau, epsilon) {
    start + (end - start) / (1 + exp(-(t - tau) / epsilon))
}

for (i in 1:nrow(truth_table)) {
    row <- truth_table[i, ]
    t <- pseudotime

    base_A <- runif(1, 20, 60) * lib_factors * rlnorm(ncol(sim_A), 0, 0.3)
    val_A <- base_A
    val_B <- numeric(ncol(sim_A))

    if (row$Type == "Null") {
        val_B <- val_A * 1.0 * rlnorm(ncol(sim_A), 0, 0.15)
    } else if (row$Type == "Switch") {
        trend <- logistic(t, 1, row$EffectSize, row$Tau, row$Epsilon)
        val_B <- val_A * trend * rlnorm(ncol(sim_A), 0, 0.15)
    } else if (row$Type == "Decouple") {
        is_after <- t > row$Tau
        val_B_coupled <- val_A * 1.0 * rlnorm(ncol(sim_A), 0, 0.15)
        val_B_decoupled <- runif(ncol(sim_A), 5, 50) * lib_factors
        val_B <- ifelse(is_after, val_B_decoupled, val_B_coupled)
    }

    counts_implanted[row$GeneA, ] <- rpois(ncol(sim_A), val_A)
    counts_implanted[row$GeneB, ] <- rpois(ncol(sim_A), val_B)
}

counts(sim_A) <- counts_implanted
colData(sim_A)$Pseudotime <- pseudotime
metadata(sim_A)$GroundTruth <- truth_table

# 6. SAVE ARTIFACTS --------------------------------------------------------
dir.create("paper_reproducibility/output", recursive = TRUE, showWarnings = FALSE)
dir.create("paper_reproducibility/figures", recursive = TRUE, showWarnings = FALSE)

saveRDS(sim_A, "paper_reproducibility/output/sim_ground_truth_lite.rds")

message(">>> Lite Simulation Saved.")
