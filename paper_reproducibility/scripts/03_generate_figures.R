library(ggplot2)

# Load Metrics
metrics <- read.csv("paper_reproducibility/output/benchmark_metrics.csv")
print(metrics)

# Create Bar Plot
pdf("paper_reproducibility/figures/fig_benchmark_metrics.pdf", width = 6, height = 4)
p <- ggplot(metrics[metrics$Metric != "ScanTimeMins", ], aes(x = Metric, y = Value, fill = Metric)) +
    geom_bar(stat = "identity") +
    ylim(0, 1) +
    theme_minimal() +
    labs(title = "DyProp Performance (Lite Benchmark)", y = "Score") +
    geom_text(aes(label = round(Value, 2)), vjust = -0.5)

print(p)
dev.off()

message(">>> Figure saved to paper_reproducibility/figures/fig_benchmark_metrics.pdf")
