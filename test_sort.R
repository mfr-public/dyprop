library(dplyr)
library(SingleCellExperiment)

suppressWarnings({
    dobj <- readRDS("paper_reproducibility/output/sim_ground_truth_classified.rds")
    sim <- readRDS("paper_reproducibility/output/sim_ground_truth.rds")

    truth_table <- metadata(sim)$GroundTruth
    truth_positives <- truth_table %>% dplyr::filter(Type != "Null")

    candidates <- dobj@events[!is.na(dobj@events$FDR) & dobj@events$FDR <= 0.05, ]
    candidates <- candidates[candidates$Event_Class %in% c("Switch", "Decoupling"), ]

    # Sort mathematically
    candidates_sorted <- candidates[order(candidates$FDR, -candidates$Tau), ]
    top_5k <- head(candidates_sorted, 5000)

    # Sort 'randomly' due to missing Score parameter (Previous Logic)
    candidates_random <- candidates[order(candidates$Score, decreasing=TRUE), ]
    top_5k_random <- head(candidates_random, 5000)

    gene_map <- setNames(seq_len(nrow(sim)), rownames(sim))
    make_key <- function(idx1, idx2) paste(pmin(idx1, idx2), pmax(idx1, idx2), sep="_")

    top_5k$Key <- make_key(gene_map[top_5k$GeneA], gene_map[top_5k$GeneB])
    top_5k_random$Key <- make_key(gene_map[top_5k_random$GeneA], gene_map[top_5k_random$GeneB])
    truth_positives$Key <- make_key(truth_positives$GeneA, truth_positives$GeneB)

    cat("\n=== DIAGNOSTIC 1: SUBSETTING LOGIC ===\n")
    cat("If the sorting logic was destructive, the 'Sorted' parameter will have 0 true positives natively captured.\n")
    cat("Base Total True Dynamic Events Injected: 1000\n")  # Out of 1500 total, 500 are likely null
    cat("True Positives captured in Top 5000 (Random/NULL Score): ", sum(top_5k_random$Key %in% truth_positives$Key), "\n")
    cat("True Positives captured in Top 5000 (Sorted by FDR/Tau): ", sum(top_5k$Key %in% truth_positives$Key), "\n")

    cat("\n=== DIAGNOSTIC 2: GLMM DOUBLE-FDR PENALTY ===\n")
    if (file.exists("paper_reproducibility/output/benchmark_tractability.csv")) {
        res <- read.csv("paper_reproducibility/output/benchmark_tractability.csv")
        print(res)
    }
})
