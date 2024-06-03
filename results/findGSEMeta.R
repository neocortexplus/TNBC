library(GEOquery)
library(limma)
library(umap)
library(tidyr)
library(dplyr)
library(ggplot2)
library(Biobase)
library(readxl)

# Load the environments and data
env_GSE76250 <- new.env()
env_GSE38959 <- new.env()

load("~/RCode/TNBC/results/GSE76250_tT2.RData", envir = env_GSE76250)
load("~/RCode/TNBC/results/GSE38959_tT2.RData", envir = env_GSE38959)

GSE76250_tT2 <- env_GSE76250$tT2
GSE38959_tT2 <- env_GSE38959$tT2


p_value_threshold <- 0.05

# Apply criteria for upregulated and downregulated
upregulated_GSE76250 <- GSE76250_tT2[GSE76250_tT2$logFC >= 2 & GSE76250_tT2$adj.P.Val < p_value_threshold, ]
downregulated_GSE76250 <- GSE76250_tT2[GSE76250_tT2$logFC <= -2 & GSE76250_tT2$adj.P.Val < p_value_threshold, ]

upregulated_GSE38959 <- GSE38959_tT2[GSE38959_tT2$logFC >= 4 & GSE38959_tT2$adj.P.Val < p_value_threshold, ]
downregulated_GSE38959 <- GSE38959_tT2[GSE38959_tT2$logFC  <=-4 & GSE38959_tT2$adj.P.Val < p_value_threshold, ]
dim(downregulated_GSE38959)
# Find common genes
common_upregulated_genes <- intersect(upregulated_GSE76250$Gene.Symbol, upregulated_GSE38959$Gene.Symbol)
common_downregulated_genes <- intersect(downregulated_GSE76250$Gene.Symbol, downregulated_GSE38959$Gene.Symbol)

# Print common genes
# cat("Common Upregulated Genes:\n",sep="")
cat(common_upregulated_genes,sep="\n")
# cat("Common Downregulated Genes:\n",sep="")
cat(common_downregulated_genes,sep="\n")





# Gene semantic similarity measurement #####

g1 <- c("84842", "2524", "10590", "3070", "91746")
g2 <- c("84289", "6045", "56999", "9869")

DOSE::geneSim(g1[1], g2[1], measure="Wang", combine="BMA")

# Gene cluster semantic similarity measurement ####

g3 <- c("57491", "6296", "51438", "5504", "27319", "1643")
clusters <- list(a=g1, b=g2, c=g3)
DOSE::mclusterSim(clusters, measure="Wang", combine="BMA")

# Gene semantic similarity measurement between Genes #####

library(DOSE)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(reshape2)
library(viridis)  # For better color palettes


gene_symbols <- c("BRCA1", "TP53")  
entrez_ids <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
gene_sim_result <- geneSim(entrez_ids$ENTREZID[1], entrez_ids$ENTREZID[2], measure="Wang", combine="BMA")

# Gene cluster semantic similarity measurement between Genes####

g1 <- c("BRCA1", "TP53", "MDM2")  # Example genes in cluster a
g2 <- c("CDK2", "CCNA2", "CDK4")   # Example genes in cluster b
g3 <- c("BCL2", "CASP8", "FAS")    # Example genes in cluster c

clusters <- list(a=g1, b=g2, c=g3)

convert_to_entrez <- function(gene_list) {
  entrez <- tryCatch({
    bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  }, error = function(e) data.frame())
  return(entrez$ENTREZID)
}

clusters_entrez <- lapply(clusters, convert_to_entrez)

similarity_matrix <- mclusterSim(clusters_entrez, measure="Wang", combine="BMA")

heatmap(as.matrix(similarity_matrix), symm = TRUE, 
        main = "Semantic Similarity between Gene Clusters", 
        xlab = "Clusters", ylab = "Clusters", 
        col = colorRampPalette(c("blue", "white", "red"))(256))


melted_similarity_matrix <- melt(as.matrix(similarity_matrix))

# Name the columns appropriately
colnames(melted_similarity_matrix) <- c("X1", "X2", "value")

# Calculate the rank for each value (lower rank for higher values)
melted_similarity_matrix$rank <- rank(-melted_similarity_matrix$value, ties.method = "min")

# Create the heatmap using ggplot2
heatmap_plot <- ggplot(melted_similarity_matrix, aes(x = X1, y = X2, fill = value)) +
  geom_tile() +  # This adds the tiles
  geom_text(aes(label = rank), color = "black", size = 3, check_overlap = TRUE) +  # Add rank numbers
  scale_fill_viridis(option = "C", direction = -1, begin = 0, end = 1, limits = c(min(melted_similarity_matrix$value, na.rm = TRUE), max(melted_similarity_matrix$value, na.rm = TRUE))) +  # Color scale
  labs(title = "Semantic Similarity between Gene Clusters",
       x = "Clusters",
       y = "Clusters",
       fill = "Similarity") +  # Labels including legend title
  theme_minimal() +  # Minimal theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Improve x-axis label readability
        plot.title = element_text(hjust = 0.5),  # Center title
        legend.position = "right")  # Position of the legend

# Print the plot
print(heatmap_plot)

