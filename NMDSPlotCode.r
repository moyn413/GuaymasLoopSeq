# Load necessary libraries
library(phyloseq)
library(ggplot2)
library(vegan)

# Read in your phyloseq object
physeq <- readRDS("/path/to/your/phyloseq/object/")

# Extract sample data
sample_data <- as(sample_data(physeq), "data.frame")

# Discard low abundance samples
total_reads_per_sample <- colSums(otu_table(physeq))
physeq_filtered <- prune_samples(total_reads_per_sample >= 400, physeq)

# Calculate Bray-Curtis distance matrix
dist_matrix <- distance(physeq_filtered, method = "bray")

# Perform NMDS
nmds <- metaMDS(dist_matrix, k = 2, trymax = 100)

# Extract NMDS scores and combine with sample data
nmds_scores <- as.data.frame(scores(nmds))
nmds_scores$SampleID <- rownames(nmds_scores)
nmds_scores <- merge(nmds_scores, sample_data, by.x = "SampleID", by.y = "row.names")

# Define colors for each sampling site
site_colors <- c("bare" = "#AA7942", "white" = "#0070C0", "orange" = "#ED7D31")
site_colors <- site_colors[match(nmds_scores$mat_color, names(site_colors))]

# Map last two digits of SampleID to depth range
# In our data set the last two digits of the SampleID corresponded with depth range of the sample
# For example, SampleID E4872-06-01 is from 0-3 cm deep
depth_map <- c(
  "01" = "0-3cm", "02" = "3-6cm", "03" = "6-9cm", "04" = "9-12cm", 
  "05" = "12-15cm", "06" = "15-18cm", "07" = "18-21cm", 
  "08" = "21-26cm", "09" = "26-31cm", "10" = "31-36cm", "11" = "36-41cm"
)
nmds_scores$Depth <- substr(nmds_scores$SampleID, nchar(nmds_scores$SampleID)-1, nchar(nmds_scores$SampleID))
nmds_scores$DepthLabel <- depth_map[nmds_scores$Depth]

# Calculate NMDS stress
nmds_stress <- round(nmds$stress, 4)

# Plotting
plt <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = mat_color)) + 
  geom_point(size = 3) +  # Plot the points
  geom_text(aes(label = DepthLabel), vjust = -1, size = 5) +  # Add depth labels above the points
  scale_color_manual(values = site_colors) +  # Use the defined colors
  theme_minimal() + 
  labs(x = "NMDS1", y = "NMDS2", color = "Sampling Site") + 
  annotate("text", x = Inf, y = Inf, label = paste("Stress =", nmds_stress), 
           hjust = 1.1, vjust = 1.1, size = 5, color = "black", fontface = "bold") +  # Add stress value to the plot
  theme(
    legend.title = element_text(size = 14, face = "bold"),  # Legend title style
    legend.text = element_text(size = 12),  # Legend text size
    legend.position = "right",  # Legend position
    axis.text = element_text(size = 14),   # Axis font size
    axis.title = element_text(size = 16),  # Axis title font size
    plot.margin = margin(10, 10, 10, 10)   # Set the margins around the plot
  ) + 
  coord_cartesian(clip = 'off')

# Save the plot as a PDF
ggsave("nmds_plot.pdf", plot = plt, 
       path = "/path/to/where/you/want/to/save/plots/", 
       width = 12, height = 10)
# OPTIONAL: Create a dendrogram to investigate clusters on the NMDS plot

# Perform hierarchical clustering on Bray-Curtis distance matrix
hc <- hclust(dist_matrix)

# Convert the hclust object to a dendrogram
dend <- as.dendrogram(hc)

# Plot the dendrogram
dendro_plot <- dend %>%
  set("branches_lwd", 2) %>%  # Adjust branch thickness
  set("labels_cex", 0.8) %>%  # Adjust label size
  set("labels_col", "black") %>%  # Set label color
  set("hang", 0.1)  # Adjust label position for better readability
plot(dendro_plot, 
     main = "Dendrogram based on Bray-Curtis Dissimilarity", 
     xlab = "Sample", 
     ylab = "Bray-Curtis Dissimilarity")

# Save the plot as a PDF
ggsave("dendrogram_plot.pdf", 
       plot = dendro_plot, 
       width = 12, height = 10)
             
