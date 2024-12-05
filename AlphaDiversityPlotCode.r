# Load necessary libraries
library(phyloseq)
library(ggplot2)

# Read in your phyloseq object
physeq <- readRDS("/path/to/your/phyloseq/object/")

# Normalization
# One approach is to subsample to the minimum number of reads in a sample
# Another approach is to subsample to a chosen number of reads to capture most of the samples, without normalizing to a very small number of reads
# For example, if you have a sample with 10 reads, subsampling your other samples with hundreds to thousands of reads would skew diversity analyses

# In this example, we subsample to 400 reads per sample using rarefy_even_depth
# Samples with less than 400 reads are discarded
physeq_400 <- rarefy_even_depth(physeq, sample.size = 400, rngseed = 42)

# Transform the sample counts using the 'ceiling' function to ensure that all numbers are integers.
ps.diversity <- transform_sample_counts(physeq_400, ceiling)

# Perform diversity calculations using estimate_richness function
rich <- estimate_richness(ps.diversity, measures = c("Observed", "Shannon", "Simpson", "InvSimpson", "Chao1"))

# Add the grouping variable to the richness data
# Here "mat_color" is used
rich$mat_color <- sample_data(ps.diversity)$mat_color

# Plot alpha diversity using plot_richness and facet by 'mat_color'
all.div <- plot_richness(ps.diversity,  
                         measures = c("Observed", "Shannon", "Simpson", "InvSimpson", "Chao1"),  
                         x = "mat_color",  
                         color = "mat_color") +
  ggplot2::geom_boxplot() +  
  ggplot2::scale_color_manual(values = c("#AA7942", "#ED7D31", "#0070C0")) +
  labs(x = "Sampling Site", color = "Sampling Site")

# Display the diversity plot
print(all.div)

# Save the plot as a PDF
ggsave("alpha_diversity_plot.pdf", plot = all.div, 
       path = "/path/to/where/you/want/to/save/plots/", 
       width = 12, height = 10)

# Perform Kruskal-Wallis test (or another test) for each diversity measure
kruskal_results <- list()

# Loop over each alpha diversity measure (excluding 'mat_color')
for (measure in colnames(rich)[-ncol(rich)]) {  
  # Perform Kruskal-Wallis test comparing each diversity measure across 'mat_color' (sampling sites)
  kruskal_results[[measure]] <- kruskal.test(rich[[measure]] ~ rich$mat_color)
}

# Print the Kruskal-Wallis test results (p-values) for each metric
cat("Kruskal-Wallis Test Results for Each Metric:\n")
for (measure in names(kruskal_results)) {
  cat("\nKruskal-Wallis Test for", measure, "\n")
  print(kruskal_results[[measure]])
}

# Extract p-values and format them
p_values <- sapply(kruskal_results, function(x) x$p.value)
cat("\nKruskal-Wallis p-values for each diversity metric:\n")
print(p_values)