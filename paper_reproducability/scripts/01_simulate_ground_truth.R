# ==============================================================================
# Dynamic Proportionality - "Nature Methods" Scale Benchmark
# Topology: Rooted Bifurcation (Stem -> A / Stem -> B)
# Scale: 5,000 Cells, 10,000 Genes
# ==============================================================================

library(splatter)
library(scater)
library(SingleCellExperiment)
library(dplyr)

set.seed(2025)

# 1. CONFIGURATION (The "Gold Standard" Scale) ---------------------------------
N_CELLS <- 5000
N_GENES <- 10000 

# PAIR CONFIGURATION
# We will implant 500 of each type into "Lineage A" to test detection
N_PAIRS <- 500 
N_IMPLANTS <- N_PAIRS * 3  # 1500 Genes involved in pairs

# 2. GENERATE BIFURCATING BACKBONE ---------------------------------------------
# We use 'splatter' paths to create a Y-shape
message(">>> Generating Bifurcating Backbone (N=", N_CELLS, ")...")

params <- newSplatParams()
params <- setParams(params, list(
    #SCALE
    batchCells = N_CELLS,		# tradeSeq x2
        nGenes = N_GENES,
    
    # TOPOLOGY: Y-Shape
    #method = "paths",				# same as tradeSeq
    group.prob = c(0.2, 0.4, 0.4), # 20% Root, 40% Path A, 40% Path B
    path.from = c(0, 1, 1),        # Path 1 is Root. Path 2&3 come from 1.
    path.nSteps = c(50, 100, 100), # Root is shorter, Branches are long
    
    # NOISE MODEL (Crucial for beating Correlation)
    lib.loc = 11,				# same a tradeSeq (standard)
    lib.scale = 0.5,             # High Depth Variance (tradeseq 0.2) forces artifacts
    dropout.type = "experiment", # Realistic Zeros (potentially harder than tradeseq)
    dropout.mid = 0,
    dropout.shape = -1,
    # BIOLOGY
    bcv.common = 0.2             # Noisier than tradeSeq (0.1)
    
))

sim <- splatSimulate(params, method = "path", verbose = T, sparsify = F)

# 3. EXTRACT LINEAGE A (The Analysis Target) -----------------------------------
# We focus the benchmark on the trajectory: Root -> Lineage A
# This gives us a clean 0 -> 1 pseudotime vector for the math.

# Splatter stores group info in 'Group'
# Path 1 = Group1, Path 2 = Group2, Path 3 = Group3
cell_subset <- colData(sim)$Group %in% c("Path1", "Path2")
sim_A <- sim[, cell_subset]

# Extract Pseudotime for Lineage A
# Splatter 'Step' needs to be concatenated from Root to Branch
counts_A <- counts(sim_A)
meta_A <- colData(sim_A)

# Normalize Pseudotime (0 to 1)
# Root (Path1) goes 0->0.33, Branch (Path2) goes 0.33->1.0
meta_A$GlobalTime <- NA
is_root <- meta_A$Group == "Path1"
is_branch <- meta_A$Group == "Path2"

# Normalize Root Steps
steps_root <- meta_A$Step[is_root]
meta_A$GlobalTime[is_root] <- (steps_root / max(steps_root)) * 0.33

# Normalize Branch Steps
steps_branch <- meta_A$Step[is_branch]
meta_A$GlobalTime[is_branch] <- 0.33 + ((steps_branch / max(steps_branch)) * 0.67)

pseudotime <- meta_A$GlobalTime
lib_factors <- colSums(counts_A) / mean(colSums(counts_A))

message(">>> Lineage A Extracted. N=", ncol(sim_A), " Cells.")

# 4. IMPLANTATION: ROBUST & SCALABLE ------------------------------------------
message(">>> Implanting ", N_IMPLANTS, " Ground Truth Pairs...")

truth_table <- data.frame(
    PairID = paste0("Pair_", 1:N_IMPLANTS),
    GeneA = seq(1, by=2, length.out=N_IMPLANTS),
    GeneB = seq(2, by=2, length.out=N_IMPLANTS),
    Type = rep(c("Null", "Switch", "Decouple"), each=N_PAIRS),
    # Randomize Parameters
    Tau = runif(N_IMPLANTS, 0.4, 0.8),      # Events happen in the Branch (after 0.33)
    Epsilon = runif(N_IMPLANTS, 0.03, 0.15),
    EffectSize = runif(N_IMPLANTS, 3, 8)
)

counts_implanted <- counts_A

logistic <- function(t, start, end, tau, epsilon) {
    start + (end - start) / (1 + exp(-(t - tau) / epsilon))
}

for(i in 1:nrow(truth_table)) {
    row <- truth_table[i, ]
    t <- pseudotime
    
    # Base Gene A (Constitutive background)
    base_A <- runif(1, 20, 60) * lib_factors * rlnorm(ncol(sim_A), 0, 0.3)
    val_A <- base_A
    val_B <- numeric(ncol(sim_A))
    
    if(row$Type == "Null") {
        # Fixed Ratio (Stable)
        val_B <- val_A * 1.0 * rlnorm(ncol(sim_A), 0, 0.15)
        
    } else if(row$Type == "Switch") {
        # Changing Mean
        trend <- logistic(t, 1, row$EffectSize, row$Tau, row$Epsilon)
        val_B <- val_A * trend * rlnorm(ncol(sim_A), 0, 0.15)
        
    } else if(row$Type == "Decouple") {
        # Changing Variance
        is_after <- t > row$Tau
        val_B_coupled <- val_A * 1.0 * rlnorm(ncol(sim_A), 0, 0.15)
        val_B_decoupled <- runif(ncol(sim_A), 5, 50) * lib_factors # Pure Noise
        val_B <- ifelse(is_after, val_B_decoupled, val_B_coupled)
    }
    
    # Apply Poisson Noise
    counts_implanted[row$GeneA, ] <- rpois(ncol(sim_A), val_A)
    counts_implanted[row$GeneB, ] <- rpois(ncol(sim_A), val_B)
    
    if(i %% 100 == 0) message("   Processed ", i, "/", N_IMPLANTS)
}

# Save Result
counts(sim_A) <- counts_implanted
colData(sim_A)$Pseudotime <- pseudotime
metadata(sim_A)$GroundTruth <- truth_table

# 5. OUTPUT SUMMARY -----------------------------------------------------------
message(">>> Simulation Ready.")
message("Object: sim_A")
message("Cells: ", ncol(sim_A), " (Root + Branch A)")
message("Genes: ", nrow(sim_A))
message("Ground Truth: metadata(sim_A)$GroundTruth")

# Optional: Visual Check of one Switch Pair
 df <- data.frame(t=pseudotime, A=counts_implanted[1,], B=counts_implanted[2,])
 plot(df$t, log2(df$A/df$B))
 
 # Check a Switch Pair (Pair 501)
 # Indices: (501 * 2) - 1 = 1001
 df_switch <- data.frame(t=pseudotime, A=counts_implanted[1001,], B=counts_implanted[1002,])
 plot(df_switch$t, log2((df_switch$A+1)/(df_switch$B+1)), main="The Switch (Pair 501)")
 abline(v=0.5, col="red", lty=2) # Expected Tipping Point approx here

 # Check a Decoupling Pair (Pair 1001)
 # Indices: (1001 * 2) - 1 = 2001
 df_decouple <- data.frame(t=pseudotime, A=counts_implanted[2001,], B=counts_implanted[2002,])
 plot(df_decouple$t, log2((df_decouple$A+1)/(df_decouple$B+1)), main="The Decoupling (Pair 1001)")
 
 # --- 6. SAVE ARTIFACTS --------------------------------------------------------
dir.create("paper_reproducibility/output", recursive = TRUE, showWarnings = FALSE)
dir.create("paper_reproducibility/figures", recursive = TRUE, showWarnings = FALSE)

# Save the SingleCellExperiment Object
saveRDS(sim_A, "paper_reproducibility/output/sim_ground_truth.rds")

# Save the Plots
pdf("paper_reproducibility/figures/fig_simulation_type1_switch.pdf")
plot(df_switch$t, log2((df_switch$A+1)/(df_switch$B+1)), 
     main="Type 1: Stoichiometric Switch", xlab="Pseudotime", ylab="Log-Ratio")
abline(v=metadata(sim_A)$GroundTruth$Tau[501], col="red", lty=2)
dev.off()

pdf("paper_reproducibility/figures/fig_simulation_type2_decoupling.pdf")
plot(df_decouple$t, log2((df_decouple$A+1)/(df_decouple$B+1)), 
     main="Type 2: Decoupling Event", xlab="Pseudotime", ylab="Log-Ratio")
dev.off()

message(">>> Simulation Artifacts Saved.")
 
 