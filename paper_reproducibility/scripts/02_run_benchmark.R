library(dyprop)
library(SingleCellExperiment)
library(dplyr)
library(data.table)

# 1. Load Ground Truth Data
# ==============================================================================
# rds_path <- "paper_reproducibility/output/sim_ground_truth.rds"
rds_path <- "paper_reproducibility/output/sim_ground_truth.rds"
# rds_path <- "paper_reproducibility/output/sim_ground_truth_lite.rds"
if (!file.exists(rds_path)) stop("Simulation file not found. Run 01_simulate_ground_truth.R first.")

message(">>> Loading simulation data...")
sim <- readRDS(rds_path)
counts_mat <- t(counts(sim)) # dyprop expects Samples x Genes
pseudotime <- colData(sim)$Pseudotime
truth_table <- metadata(sim)$GroundTruth

message("   Cells: ", ncol(sim))
message("   Genes: ", nrow(sim))
message("   Injected Events: ", nrow(truth_table))

# 2. Run DyProp Pipeline
# ==============================================================================
message(">>> Step 1: Preprocessing...")
# Use basic pooling if glmpca is slow for 10k genes, but for paper we use best settings.
# We'll stick to default.
dobj <- prepareComposition(counts_mat, pseudotime, k = 0)



message(">>> Step 2: Scanning...")
# Scan with standard grid
# For 10k genes, we have 50M pairs.

start_time <- Sys.time()
dobj <- scanDynamics(dobj, cores = 4)
end_time <- Sys.time()
message("   Scan Time: ", round(difftime(end_time, start_time, units = "mins"), 2), " mins")

message(">>> Step 3: FDR Estimation...")
dobj <- estimateFDR(dobj, n_permutations = 5000, cores = 4)

message(">>> Step 4: Filtering & Classifying...")
# Optimization: Only classify potential hits (FDR < 0.05)
# Cap at 2000 events to ensure reproduction script finishes quickly for this initial test
dobj@events <- dobj@events[dobj@events$FDR < 0.05, ]
if (nrow(dobj@events) > 2000) {
    message("... Capping at top 2000 events for speed.")
    dobj@events <- dobj@events[1:2000, ]
}
dobj <- classifyEvents(dobj)

# 3. Calculate Metrics (Precision / Recall)
# ==============================================================================
message(">>> Calculating Benchmarks...")

detected <- dobj@events
# Filter significant
sig_events <- detected %>% filter(FDR < 0.05)

# We need to match detected pairs (GeneA, GeneB) to Truth Table.
# Truth Table has PairID, GeneA (index), GeneB (index), Type.
# Detected has Gene Name "GeneX".

# Map detected names to indices
gene_map <- setNames(seq_len(nrow(sim)), rownames(sim))
detected$IdxA <- gene_map[detected$GeneA]
detected$IdxB <- gene_map[detected$GeneB]

# Construct Pair Key "min_max" to handle undirected matching
make_key <- function(idx1, idx2) {
    pmin <- pmin(idx1, idx2)
    pmax <- pmax(idx1, idx2)
    paste(pmin, pmax, sep = "_")
}

detected$Key <- make_key(detected$IdxA, detected$IdxB)
truth_table$Key <- make_key(truth_table$GeneA, truth_table$GeneB)

# Merge
res <- detected %>%
    inner_join(truth_table, by = "Key") %>%
    select(Key, Score, FDR, Event_Class, Type_Truth = Type, Tau_Truth = Tau)

# Evaluate Classification Accuracy
# Truth: "Switch", "Decouple", "Null"
# Predicted: "Switch", "Decoupling"

# Define Positives (exclude Null)
truth_positives <- truth_table %>% filter(Type != "Null")
n_positives <- nrow(truth_positives)

# 3a. Global Detection (Detecting Non-Null)
tp_global <- res %>% filter(Type_Truth != "Null" & FDR < 0.05)
fp_global <- detected %>% filter(FDR < 0.05 & !Key %in% truth_positives$Key)

precision <- nrow(tp_global) / (nrow(tp_global) + nrow(fp_global))
recall <- nrow(tp_global) / n_positives
f1 <- 2 * (precision * recall) / (precision + recall)

message("--- Global Results (FDR < 0.05) ---")
message("Precision: ", round(precision, 3))
message("Recall:    ", round(recall, 3))
message("F1 Score:  ", round(f1, 3))

# 3b. Classification Accuracy (Switch vs Decouple)
# Subset to correctly detected positives
res_correct <- tp_global

# Confusion Matrix
table(Predicted = res_correct$Event_Class, Truth = res_correct$Type_Truth)

# Save Metrics
metrics <- data.frame(
    Metric = c("Precision", "Recall", "F1", "ScanTimeMins"),
    Value = c(precision, recall, f1, as.numeric(difftime(end_time, start_time, units = "mins")))
)

write.csv(metrics, "paper_reproducibility/output/benchmark_metrics.csv", row.names = FALSE)
message(">>> Benchmark metrics saved to paper_reproducibility/output/benchmark_metrics.csv")
