library(GEOquery)
library(limma)
library(umap)
library(tidyr)
library(dplyr)
library(ggplot2)
library(Biobase)
library(readxl)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)

# Load the environments and data
env_GSE76250 <- new.env()
env_GSE38959 <- new.env()

load("~/RCode/TNBC/results/GSE76250_tT2.RData", envir = env_GSE76250)
load("~/RCode/TNBC/results/GSE38959_tT2.RData", envir = env_GSE38959)

GSE76250_tT2 <- env_GSE76250$tT2
GSE38959_tT2 <- env_GSE38959$tT2


p_value_threshold <- 0.05

all_regulated_GSE76250 <-GSE76250_tT2[abs(GSE76250_tT2$logFC >= 1) & GSE76250_tT2$adj.P.Val < p_value_threshold, ]

# Apply criteria for upregulated and downregulated
upregulated_GSE76250 <- GSE76250_tT2[GSE76250_tT2$logFC >= 1 & GSE76250_tT2$adj.P.Val < p_value_threshold, ]
downregulated_GSE76250 <- GSE76250_tT2[GSE76250_tT2$logFC <= -1 & GSE76250_tT2$adj.P.Val < p_value_threshold, ]

all_regulated_GSE38959 <- GSE38959_tT2[abs(GSE38959_tT2$logFC) >= 2 & GSE38959_tT2$adj.P.Val < p_value_threshold, ]

upregulated_GSE38959 <- GSE38959_tT2[GSE38959_tT2$logFC >= 1 & GSE38959_tT2$adj.P.Val < p_value_threshold, ]
downregulated_GSE38959 <- GSE38959_tT2[GSE38959_tT2$logFC  <=-1 & GSE38959_tT2$adj.P.Val < p_value_threshold, ]
dim(downregulated_GSE38959)

# Find common genes
common_all_regulated_genes <- as.data.frame(intersect(all_regulated_GSE38959,all_regulated_GSE76250))
common_upregulated_genes <- as.data.frame(intersect(upregulated_GSE76250$Gene.Symbol, upregulated_GSE38959$Gene.Symbol))
common_downregulated_genes <- as.data.frame(intersect(downregulated_GSE76250$Gene.Symbol, downregulated_GSE38959$Gene.Symbol))


# Print common genes
# cat("Common Upregulated Genes:\n",sep="")
cat(common_upregulated_genes,sep="\n")
# cat("Common Downregulated Genes:\n",sep="")
cat(common_downregulated_genes,sep="\n")



# GO enrichment analysis ####


head(upregulated_GSE38959)
head(upregulated_GSE38959$Gene.ID)
gene_vector <- setNames(GSE38959_tT2$logFC, GSE38959_tT2$Gene.ID)

ego <- enrichGO(gene          = upregulated_GSE38959$Gene.ID,
                universe      = names(gene_vector),
                OrgDb         = org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
dim(ego)
head(ego,100)
x <- as.data.frame(ego)
dim(x)


library(enrichplot)
barplot(ego, showCategory=20) 
mutate(ego, qscore = -log(p.adjust, base=10)) %>% 
  barplot(x="qscore")




# decrease text size by 20%
par(cex=0.2)  

goplot(ego)


gene.df <- bitr(upregulated_GSE38959$Gene.ID, fromType = "ENTREZID",
                toType = c("ENSEMBL", "SYMBOL"),
                OrgDb = org.Hs.eg.db)

ego2 <- enrichGO(gene         = gene.df$ENSEMBL,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'ENSEMBL',
                 ont           = "CC",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)



dim(ego2)
head(ego2,100)
x <- as.data.frame(ego2)
dim(x)


dotplot(ego, showCategory=30) + ggtitle("dotplot for ORA")
dotplot(ego2, showCategory=30) + ggtitle("dotplot for GSEA")


sorted_gene_vector <- sort(gene_vector, decreasing = TRUE)


edox <- setReadable(ego, 'org.Hs.eg.db', 'ENTREZID')
p1 <- cnetplot(edox, foldChange=gene_vector)
p2 <- cnetplot(edox, categorySize="pvalue", foldChange=sorted_gene_vector)
p3 <- cnetplot(edox, foldChange=geneList, circular = TRUE, colorEdge = TRUE) 
cowplot::plot_grid(p1, p2, p3, ncol=3, labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2))

p1 <- cnetplot(edox, node_label="category", 
               cex_label_category = 1.2) 
p2 <- cnetplot(edox, node_label="gene", 
               cex_label_gene = 0.8) 
p3 <- cnetplot(edox, node_label="all") 

p4 <- cnetplot(edox, node_label="none", 
               color_category='firebrick', 
               color_gene='steelblue') 
cowplot::plot_grid(p1, p2, p3, p4, ncol=2, labels=LETTERS[1:4])


p1 <- heatplot(edox, showCategory=5)
p2 <- heatplot(edox, foldChange=geneList, showCategory=5)
cowplot::plot_grid(p1, p2, ncol=1, labels=LETTERS[1:2])


edox2 <- pairwise_termsim(edox)
p1 <- treeplot(edox2)
p2 <- treeplot(edox2, hclust_method = "average")
aplot::plot_list(p1, p2, tag_levels='A')


edo <- pairwise_termsim(ego)
p1 <- emapplot(edo)
p2 <- emapplot(edo, cex_category=1.5)
p3 <- emapplot(edo, layout="kk")
p4 <- emapplot(edo, cex_category=1.5,layout="kk") 
cowplot::plot_grid(p1, p2, p3, p4, ncol=2, labels=LETTERS[1:4])

upsetplot(ego)

ridgeplot(ego)

# Now use the sorted vector in the gseGO function
ego3 <- gseGO(geneList     = sorted_gene_vector,
              OrgDb        = org.Hs.eg.db,
              ont          = "CC",
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)


dim(ego3)
head(ego3,100)
x <- as.data.frame(ego3)
dim(x)


kk <- enrichKEGG(gene         = as.character(upregulated_GSE38959$Gene.ID),
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
head(kk)



kk2 <- gseKEGG(geneList     = sorted_gene_vector,
               organism     = 'hsa',
               minGSSize    = 120,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
head(kk2)



mkk <- enrichMKEGG(gene = as.character(upregulated_GSE38959$Gene.ID),
                   organism = 'hsa',
                   pvalueCutoff = 1,
                   qvalueCutoff = 1)
head(mkk)                   



mkk2 <- gseMKEGG(geneList = sorted_gene_vector,
                 organism = 'hsa',
                 pvalueCutoff = 1)
head(mkk2)

browseKEGG(kk, 'hsa04110')

library("pathview")
hsa04110 <- pathview(gene.data  = sorted_gene_vector,
                     pathway.id = "hsa04110",
                     species    = "hsa",
                     limit      = list(gene=max(abs(sorted_gene_vector)), cpd=1))



w1 <- enrichWP(as.character(upregulated_GSE38959$Gene.ID), organism = "Homo sapiens") 

head(w1)

w2 <- gseWP(sorted_gene_vector, organism = "Homo sapiens")
head(w2)

library(ReactomePA)

x <- enrichPathway(gene=upregulated_GSE38959$Gene.ID, pvalueCutoff = 0.05, readable=TRUE)
head(x)


y <- gsePathway(sorted_gene_vector, 
                pvalueCutoff = 0.2,
                pAdjustMethod = "BH", 
                verbose = FALSE)
head(y)


viewPathway("E2F mediated regulation of DNA replication", 
            readable = TRUE, 
            foldChange = sorted_gene_vector)



x <- enrichDO(gene          = upregulated_GSE38959$Gene.ID,
              ont           = "DO",
              pvalueCutoff  = 0.05,
              pAdjustMethod = "BH",
              universe      = names(sorted_gene_vector),
              minGSSize     = 5,
              maxGSSize     = 500,
              qvalueCutoff  = 0.01,
              readable      = FALSE)
head(x,16)
dim(x)


ncg <- enrichNCG(upregulated_GSE38959$Gene.ID) 
head(ncg)


dgn <- enrichDGN(upregulated_GSE38959$Gene.ID) 
head(dgn)

y <- gseDO(sorted_gene_vector,
           minGSSize     = 120,
           pvalueCutoff  = 0.2,
           pAdjustMethod = "BH",
           verbose       = FALSE)
head(y, 3)


ncg <- gseNCG(sorted_gene_vector,
              pvalueCutoff  = 0.5,
              pAdjustMethod = "BH",
              verbose       = FALSE)
ncg <- setReadable(ncg, 'org.Hs.eg.db')
head(ncg, 3) 



dgn <- gseDGN(sorted_gene_vector,
              pvalueCutoff  = 0.2,
              pAdjustMethod = "BH",
              verbose       = FALSE)
dgn <- setReadable(dgn, 'org.Hs.eg.db')
head(dgn, 3) 

library(AnnotationHub)
library(MeSHDbi)

ah <- AnnotationHub(localHub=TRUE)
hsa <- query(ah, c("MeSHDb", "Homo sapiens"))
file_hsa <- hsa[[1]]
db <- MeSHDbi::MeSHDb(file_hsa)

library(meshes)
x <- enrichMeSH(upregulated_GSE38959$Gene.ID, MeSHDb = db, database='gendoo', category = 'C')
head(x)

y <- gseMeSH(sorted_gene_vector, MeSHDb = db, database = 'gene2pubmed', category = "C")
head(y)

x1 <- enricher(gene, TERM2GENE = cells)
head(x1)

cell_marker_data <- readxl::read_excel('~/RCode/TNBC/Cell_marker_Human.xlsx')

## instead of `cellName`, users can use other features (e.g. `cancerType`)
cells <- cell_marker_data %>%
  dplyr::select(cell_name, GeneID) %>%
  dplyr::mutate(GeneID = strsplit(as.character(GeneID), ',\\s*')) %>%
  tidyr::unnest(GeneID)

x <- enricher(upregulated_GSE38959$Gene.ID, TERM2GENE = cells)
head(x)

y <- GSEA(sorted_gene_vector, TERM2GENE = cells)
head(y)


library(msigdbr)
msigdbr_show_species()


m_t2g <- msigdbr(species = "Homo sapiens", category = "C6") %>% 
  dplyr::select(gs_name, entrez_gene)
head(m_t2g)

em <- enricher(upregulated_GSE38959$Gene.ID, TERM2GENE=m_t2g)
head(em)

C3_t2g <- msigdbr(species = "Homo sapiens", category = "C3") %>% 
  dplyr::select(gs_name, entrez_gene)
head(C3_t2g)

em2 <- GSEA(sorted_gene_vector, TERM2GENE = C3_t2g)
head(em2)



m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)
head(m_t2g)

em <- enricher(upregulated_GSE38959$Gene.ID, TERM2GENE=m_t2g)
head(em)

em2 <- GSEA(sorted_gene_vector, TERM2GENE = m_t2g)
head(em2)









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

g1 <- c("BRCA1", "TP53", "MDM2") 
g2 <- c("CDK2", "CCNA2", "CDK4")   
g3 <- c("BCL2", "CASP8", "FAS")   

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


# Load necessary libraries
library(ggplot2)
library(reshape2)
library(viridis)  # For better color palettes

# Assuming similarity_matrix is already defined
# Transform the matrix to a data frame suitable for ggplot2
melted_similarity_matrix <- melt(as.matrix(similarity_matrix))

# Name the columns appropriately
colnames(melted_similarity_matrix) <- c("X1", "X2", "value")

# Create the heatmap using ggplot2
heatmap_plot <- ggplot(melted_similarity_matrix, aes(x = X1, y = X2, fill = value)) +
  geom_tile() +  # This adds the tiles
  geom_text(aes(label = sprintf("%.2f", value)), color = "black", size = 3, check_overlap = TRUE) +  # Add value labels formatted to two decimal places
  scale_fill_viridis(option = "C", direction = -1, begin = 0, end = 1, limits = c(0, 1)) +  # Color scale with values between 0 and 1
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



# GO enrichment analysis #############

library(clusterProfiler)
library(org.Hs.eg.db)

data(geneList, package="DOSE")
gene <- names(geneList)[abs(geneList) > 2]

# Entrez gene ID
head(gene)

ggo <- groupGO(gene     = gene,
               OrgDb    = org.Hs.eg.db,
               ont      = "CC",
               level    =3,
               readable = TRUE)

head(ggo)



ego <- enrichGO(gene          = gene,
                universe      = names(geneList),
                OrgDb         = org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego)


gene.df <- bitr(gene, fromType = "ENTREZID",
                toType = c("ENSEMBL", "SYMBOL"),
                OrgDb = org.Hs.eg.db)

#  readable=TRUE / FALSE

ego2 <- enrichGO(gene         = gene.df$ENSEMBL,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'ENSEMBL',
                 ont           = "CC",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05
                 )
head(ego2, 10)                


ego3 <- gseGO(geneList     = geneList,
              OrgDb        = org.Hs.eg.db,
              ont          = "CC",
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)





goplot(ego2)



# KEGG enrichment analysis #####

library(clusterProfiler)
search_kegg_organism('ece', by='kegg_code')

ecoli <- search_kegg_organism('Escherichia coli', by='scientific_name')
dim(ecoli)


head(ecoli)


data(geneList, package="DOSE")
gene <- names(geneList)[abs(geneList) > 2]

kk <- enrichKEGG(gene         = gene,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
head(kk)


kk2 <- gseKEGG(geneList     = geneList,
               organism     = 'hsa',
               minGSSize    = 120,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
head(kk2)


mkk <- enrichMKEGG(gene = gene,
                   organism = 'hsa',
                   pvalueCutoff = 1,
                   qvalueCutoff = 1)
head(mkk)                   



mkk2 <- gseMKEGG(geneList = geneList,
                 organism = 'hsa',
                 pvalueCutoff = 1)
head(mkk2)


library(DOSE)
data(geneList)
de <- names(geneList)[abs(geneList) > 2]

edo <- enrichDGN(de)
library(enrichplot)
barplot(edo, showCategory=20) 
mutate(edo, qscore = -log(p.adjust, base=10)) %>% 
  barplot(x="qscore",showCategory=20)



edo2 <- gseDO(geneList)
head(edo2,10)
dotplot(edo, showCategory=10) + ggtitle("dotplot for ORA")
dotplot(edo2, showCategory=10) + ggtitle("dotplot for GSEA")

## convert gene ID to Symbol
edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
p1 <- cnetplot(edox, foldChange=geneList)
## categorySize can be scaled by 'pvalue' or 'geneNum'
p2 <- cnetplot(edox, categorySize="pvalue", foldChange=geneList)
p3 <- cnetplot(edox, foldChange=geneList, circular = TRUE, colorEdge = TRUE) 
cowplot::plot_grid(p1, ncol=1, labels=LETTERS[1])

p1_labeled <- cowplot::plot_grid(p1, labels=LETTERS[1])
print(p1_labeled)

p2_labeled <- cowplot::plot_grid(p2, labels=LETTERS[2])
print(p2_labeled)


p3_labeled <- cowplot::plot_grid(p3, labels=LETTERS[3])
print(p3_labeled)


library(ggplot2)
p3 <- cnetplot(edox, foldChange=geneList, circular = TRUE, colorEdge = TRUE) 

# Assuming p3 is a ggplot object
p3 <- p3 + theme(
  plot.title = element_text(size = 3),    # Smaller title text
  axis.text.x = element_text(angle = 12, hjust = 1)  # Tilt x-axis labels
)

# Redo the plot grid with adjusted plot
p3_labeled <- cowplot::plot_grid(p3, labels=LETTERS[3])
print(p3_labeled)



library(clusterProfiler)
library(ggplot2)

# Generate the plot with adjustments
p3 <- cnetplot(edox, foldChange=geneList, circular = TRUE, colorEdge = TRUE) +
  theme(
    text = element_text(size = 2),  # Smaller text size for all text elements
    plot.title = element_text(size = 5)  # Adjust title size separately if needed
  )


# Open a larger plot window in RStudio before plotting
dev.new(width = 10, height = 10)

# Plot and print
print(p3)


p1 <- cnetplot(edox, node_label="category", 
               cex_label_category = 1.2) 
p2 <- cnetplot(edox, node_label="gene", 
               cex_label_gene = 0.8) 
p3 <- cnetplot(edox, node_label="all") 
p4 <- cnetplot(edox, node_label="none", 
               color_category='firebrick', 
               color_gene='steelblue') 
cowplot::plot_grid(p1, p2, p3, p4, ncol=2, labels=LETTERS[1:4])




p1 <- heatplot(edox, showCategory=5)
p2 <- heatplot(edox, foldChange=geneList, showCategory=5)
cowplot::plot_grid(p1, p2, ncol=1, labels=LETTERS[1:2])


edox2 <- pairwise_termsim(edox)
p1 <- treeplot(edox2)
p2 <- treeplot(edox2, hclust_method = "average")
aplot::plot_list(p1, p2, tag_levels='A')


library(clusterProfiler)  # Ensure the necessary package is loaded for treeplot
library(gridExtra)

# Adjust plot parameters before plotting
par(mar = c(5, 8, 4, 2) + 0.1)  # Adjust margins: bottom, left, top, right

# Generate plots with customized settings
p1 <- treeplot(edox2, cex = 0.01, margin = 1)  # Adjust text size and margins
print(p1)
p2 <- treeplot(edox2, hclust_method = "complete", cex = 0.01, margin = 1)  # Try complete method with adjusted margin
print(p2)

grid.arrange(p1, p2, ncol = 2)


library(cowplot)

# Using plot_grid to combine plots
plot_grid(p1, p2, labels = c("A", "B"), ncol = 2)

edo <- pairwise_termsim(edo)
p1 <- emapplot(edo, layout="fr", cex_category=0.05, cex_label_category=0.4, cex_line=0.1)
p1 <- emapplot(edo, layout="graphopt",cex_line=0.7,cex_label_category=0.6)
print(p1)

p2 <- emapplot(edo, cex_category=0.05)
p3 <- emapplot(edo, layout="kk")
p4 <- emapplot(edo, cex_category=0.5,layout="kk") 
cowplot::plot_grid(p1, p2, p3, p4, ncol=2, labels=LETTERS[1:4])



