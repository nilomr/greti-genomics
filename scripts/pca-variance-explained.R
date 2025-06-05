# Import configuration file
config <- config::get()
box::use(R / plot[titheme, base_breaks_x, base_breaks_y])

# Parameters (matching the RF script)
num.eigenvectors <- 15
num.threads <- 44

# Input files
all.gds <- file.path(config$path$data, "all_recode.gds")

# Load all the data
all.maf <- SNPRelate::snpgdsOpen(all.gds, allow.duplicate = TRUE)

# Get all sample IDs
all.id <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(all.maf, "sample.id"))

cat("Performing PCA on all samples to calculate variance explained...\n")

# Perform PCA on all samples
pca_full <- SNPRelate::snpgdsPCA(
    all.maf,
    sample.id = all.id,
    num.thread = num.threads,
    verbose = TRUE
)

# Calculate variance explained by each PC - FIXED
eigenvalues <- pca_full$eigenval
# Use varprop which gives the correct proportion of variance
variance_explained <- pca_full$varprop * 100  # Convert to percentage
cumulative_variance <- cumsum(variance_explained)

# Alternative calculation using TraceXTX (should give same result)
# total_variance <- pca_full$TraceXTX
# variance_explained_alt <- eigenvalues / total_variance * 100

# Create results dataframe
pca_results <- data.frame(
    PC = 1:length(eigenvalues),
    eigenvalue = eigenvalues,
    variance_explained = variance_explained,
    cumulative_variance = cumulative_variance
)

# Report variance explained by the first num.eigenvectors components
cat("\n=== VARIANCE EXPLAINED SUMMARY ===\n")
cat(sprintf("Number of eigenvectors used in RF: %d\n", num.eigenvectors))
cat(sprintf("Variance explained by PC1-PC%d: %.2f%%\n", 
    num.eigenvectors, cumulative_variance[num.eigenvectors]))
cat(sprintf("Individual PC contributions (PC1-PC%d):\n", num.eigenvectors))

for (i in 1:num.eigenvectors) {
    cat(sprintf("  PC%d: %.2f%%\n", i, variance_explained[i]))
}

# Create scree plot - FIXED
scree_plot <- ggplot2::ggplot(
    pca_results[1:min(50, nrow(pca_results)), ], 
    ggplot2::aes(x = PC, y = variance_explained)
) +
    ggplot2::geom_col(fill = "#5d8566", alpha = 0.7) +
    ggplot2::geom_vline(
        xintercept = num.eigenvectors, 
        linetype = "dashed", 
        color = "#c29007",
        linewidth = 1
    ) +
    ggplot2::annotate(
        "text",
        x = num.eigenvectors + 2,
        y = max(variance_explained[1:min(50, length(variance_explained))]) * 0.8,
        label = paste0("Used in RF\n(", num.eigenvectors, " PCs)"),
        color = "#c29007",
        size = 3
    ) +
    ggplot2::labs(
        title = "Scree plot: Variance explained by principal components",
        x = "Principal Component",
        y = "Variance Explained (%)"
    ) +
    base_breaks_x(seq(0, 50, 5), expand = ggplot2::expansion(mult = 0.02)) +
    base_breaks_y(seq(0, ceiling(max(variance_explained, na.rm = TRUE)), 
                     max(1, ceiling(max(variance_explained, na.rm = TRUE)/10))), 
                  expand = ggplot2::expansion(mult = c(0, 0.05))) +
    titheme() +
    ggplot2::theme(aspect.ratio = 0.6)

# Create cumulative variance plot
cumulative_plot <- ggplot2::ggplot(
    pca_results[1:min(50, nrow(pca_results)), ], 
    ggplot2::aes(x = PC, y = cumulative_variance)
) +
    ggplot2::geom_line(color = "#5d8566", linewidth = 1) +
    ggplot2::geom_point(color = "#5d8566", size = 1) +
    ggplot2::geom_vline(
        xintercept = num.eigenvectors, 
        linetype = "dashed", 
        color = "#c29007",
        linewidth = 1
    ) +
    ggplot2::geom_hline(
        yintercept = cumulative_variance[num.eigenvectors],
        linetype = "dotted",
        color = "#c29007",
        linewidth = 0.8
    ) +
    ggplot2::annotate(
        "text",
        x = num.eigenvectors + 5,
        y = cumulative_variance[num.eigenvectors] - 5,
        label = sprintf("%.1f%% variance\nexplained", cumulative_variance[num.eigenvectors]),
        color = "#c29007",
        size = 3
    ) +
    ggplot2::labs(
        title = "Cumulative variance explained by principal components",
        x = "Principal Component",
        y = "Cumulative Variance Explained (%)"
    ) +
    base_breaks_x(seq(0, 50, 5), expand = ggplot2::expansion(mult = 0.02)) +
    base_breaks_y(seq(0, 100, 10), expand = ggplot2::expansion(mult = c(0, 0.02))) +
    titheme() +
    ggplot2::theme(aspect.ratio = 0.6)

# Save plots
ggplot2::ggsave(
    file.path(config$path$figures, "pca_scree_plot.jpg"),
    plot = scree_plot,
    width = 10, height = 6, bg = "white"
)

ggplot2::ggsave(
    file.path(config$path$figures, "pca_cumulative_variance.jpg"),
    plot = cumulative_plot,
    width = 10, height = 6, bg = "white"
)

# Save results for paper reporting
pca_summary <- list(
    num_eigenvectors_used = num.eigenvectors,
    variance_explained_by_selected = cumulative_variance[num.eigenvectors],
    individual_pc_contributions = variance_explained[1:num.eigenvectors],
    total_pcs = length(eigenvalues),
    eigenvalues = eigenvalues[1:num.eigenvectors]
)

# Save detailed results
write.csv(pca_results, 
    file.path(config$path$data, "pca_variance_results.csv"), 
    row.names = FALSE)

saveRDS(pca_summary, 
    file.path(config$path$data, "pca_summary_for_paper.rds"))

# Print summary for paper
cat("\n=== FOR PAPER REPORTING ===\n")
cat(sprintf("The first %d principal components explained %.1f%% of the total genetic variance.\n",
    num.eigenvectors, cumulative_variance[num.eigenvectors]))

# Close GDS file
SNPRelate::snpgdsClose(all.maf)

# cat("\nAnalysis complete. Results saved to:\n")
# cat(sprintf("- %s\n", file.path(config$path$data, "pca_variance_results.csv")))
# cat(sprintf("- %s\n", file.path(config$path$data, "pca_summary_for_paper.rds")))
# cat(sprintf("- %s\n", file.path(config$path$figures, "pca_scree_plot.jpg")))
# cat(sprintf("- %s\n", file.path(config$path$figures, "pca_cumulative_variance.jpg")))
