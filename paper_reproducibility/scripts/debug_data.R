library(dyprop)
library(SingleCellExperiment)

rds_path <- "paper_reproducibility/output/sim_ground_truth_lite.rds"
sim <- readRDS(rds_path)
counts_mat <- t(counts(sim))
pseudotime <- colData(sim)$Pseudotime

message("Counts Range: ", min(counts_mat), " - ", max(counts_mat))
message("Counts Variance: ", mean(apply(counts_mat, 2, var)))
message("Zero fraction: ", mean(counts_mat == 0))

# Try manual prep
dobj <- prepareComposition(counts_mat, pseudotime, k = 5)
message("CLR Range: ", min(dobj@matrix), " - ", max(dobj@matrix))
message("CLR NAs: ", sum(is.na(dobj@matrix)))

# Try mini scan
dobj <- scanDynamics(dobj, cores = 1)
if (!is.null(dobj@events)) {
    message("Events found: ", nrow(dobj@events))
    if (nrow(dobj@events) > 0) print(head(dobj@events))
} else {
    message("No events slot populated.")
}
