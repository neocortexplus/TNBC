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
env_meta <- new.env()


load("~/RCode/TNBC/results/GSE76250_tT2.RData", envir = env_GSE76250)
load("~/RCode/TNBC/results/GSE38959_tT2.RData", envir = env_GSE38959)
load("~/RCode/TNBC/meta_tT2.RData",envir = env_meta)


GSE76250_tT2 <- env_GSE76250$tT2
GSE38959_tT2 <- env_GSE38959$tT2
meta <- env_meta$tT2



# Convert Gene Symbols to Entrez IDs, including handling of multiple mappings
gene_conversion <- tryCatch({
  bitr(GSE76250_tT2$Gene.Symbol, fromType = "SYMBOL",
       toType = "ENTREZID", OrgDb = org.Hs.eg.db, drop = FALSE)
}, error = function(e) {
  message("Error during conversion: ", e)
  return(data.frame(SYMBOL = character(), ENTREZID = character()))  # return empty dataframe on error
})

# Assuming taking the first mapping if multiple mappings exist
gene_conversion <- gene_conversion[!duplicated(gene_conversion$SYMBOL), ]

GSE76250_tT2 <- merge(GSE76250_tT2, gene_conversion, by.x = "Gene.Symbol", by.y = "SYMBOL", all.x = TRUE)

# Rename the column for clarity
colnames(GSE76250_tT2)[colnames(GSE76250_tT2) == "ENTREZID"] <- "Gene.ID"

head(GSE76250_tT2)

GSE76250_tT2$ID <- GSE76250_tT2$Gene.ID
head(GSE76250_tT2)

GSE76250_tT2$ID <- as.integer(GSE76250_tT2$ID)
GSE76250_tT2 <- GSE76250_tT2[!is.na(GSE76250_tT2$Gene.ID), ]

summary(GSE76250_tT2);summary(GSE38959_tT2)


p_value_threshold <- 0.05

all_regulated_GSE76250 <- GSE76250_tT2[abs(GSE76250_tT2$logFC) >= 1 & GSE76250_tT2$adj.P.Val < p_value_threshold, ]
upregulated_GSE76250 <- GSE76250_tT2[GSE76250_tT2$logFC >= 1 & GSE76250_tT2$adj.P.Val < p_value_threshold, ]
downregulated_GSE76250 <- GSE76250_tT2[GSE76250_tT2$logFC <= -1 & GSE76250_tT2$adj.P.Val < p_value_threshold, ]
dim(all_regulated_GSE76250);dim(upregulated_GSE76250);dim(downregulated_GSE76250)

all_regulated_GSE38959 <- GSE38959_tT2[abs(GSE38959_tT2$logFC) >= 2 & GSE38959_tT2$adj.P.Val < p_value_threshold, ]
upregulated_GSE38959 <- GSE38959_tT2[GSE38959_tT2$logFC >= 2 & GSE38959_tT2$adj.P.Val < p_value_threshold, ]
downregulated_GSE38959 <- GSE38959_tT2[GSE38959_tT2$logFC  <=-2 & GSE38959_tT2$adj.P.Val < p_value_threshold, ]
dim(all_regulated_GSE38959);dim(upregulated_GSE38959);dim(downregulated_GSE38959)



# common_all_regulated_genes <- inner_join(all_regulated_GSE38959, all_regulated_GSE76250, by = "Gene.Symbol")
# dim(common_all_regulated_genes)

common_all_regulated_genes <- as.data.frame(intersect(all_regulated_GSE38959$Gene.Symbol , all_regulated_GSE76250$Gene.Symbol))
common_upregulated_genes <- as.data.frame(intersect(upregulated_GSE76250$Gene.Symbol, upregulated_GSE38959$Gene.Symbol))
common_downregulated_genes <- as.data.frame(intersect(downregulated_GSE76250$Gene.Symbol, downregulated_GSE38959$Gene.Symbol))
dim(common_all_regulated_genes);dim(common_upregulated_genes);dim(common_downregulated_genes)

# META

meta_GSE76250_GSE38959 <- meta[abs(meta$logFC) >= 1 & meta$adj.P.Val < p_value_threshold, ]
upregulated_meta_GSE76250_GSE38959 <- meta[meta$logFC >= 1 & meta$adj.P.Val < p_value_threshold, ]
downregulated_meta_GSE76250_GSE38959 <- meta[meta$logFC <= -1 & meta$adj.P.Val < p_value_threshold, ]
dim(meta_GSE76250_GSE38959);dim(upregulated_meta_GSE76250_GSE38959);dim(downregulated_meta_GSE76250_GSE38959)


# cat("Common Upregulated Genes:\n",sep="")
cat(common_upregulated_genes,sep="\n")
# cat("Common Downregulated Genes:\n",sep="")
cat(common_downregulated_genes,sep="\n")



# GO - Pathway Analysis ####


gene_vector1 <- setNames(GSE76250_tT2$logFC, GSE76250_tT2$Gene.ID)
gene_vector2 <- setNames(GSE38959_tT2$logFC,GSE38959_tT2$Gene.ID)
gene_vector3 <- setNames(meta$logFC,meta$Gene.ID)


ego1 <- enrichGO(gene          = all_regulated_GSE76250$Gene.ID,
                universe      = names(gene_vector1),
                OrgDb         = org.Hs.eg.db,
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

ego2 <- enrichGO(gene          = all_regulated_GSE38959$Gene.ID,
                 universe      = names(gene_vector2),
                 OrgDb         = org.Hs.eg.db,
                 ont           = "ALL",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.05,
                 readable      = TRUE)

ego3 <- enrichGO(gene          = meta_GSE76250_GSE38959$Gene.ID,
                 universe      = names(gene_vector3),
                 OrgDb         = org.Hs.eg.db,
                 ont           = "ALL",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.05,
                 readable      = TRUE)


dim(ego1);dim(ego2);dim(ego3)
head(ego1,100)
head(ego2,100)

# x <- as.data.frame(ego)
# dim(x)


top_terms_combined <- x %>%
  filter(p.adjust < 0.05) %>%   # Filter for adjusted p-values less than 0.05
  group_by(ONTOLOGY) %>%        # Group data by the ONTOLOGY column
  arrange(p.adjust) %>%         # Arrange data within each group by p.adjust
  slice_head(n = 10) %>%        # Select the top 10 entries for each group
  ungroup()                     # Remove the grouping



gene.df1 <- bitr(all_regulated_GSE76250$Gene.ID, fromType = "ENTREZID",
                toType = c("ENSEMBL", "SYMBOL"),
                OrgDb = org.Hs.eg.db)

gene.df2 <- bitr(all_regulated_GSE38959$Gene.ID, fromType = "ENTREZID",
                 toType = c("ENSEMBL", "SYMBOL"),
                 OrgDb = org.Hs.eg.db)

gene.df3 <- bitr(meta_GSE76250_GSE38959$Gene.ID, fromType = "ENTREZID",
                 toType = c("ENSEMBL", "SYMBOL"),
                 OrgDb = org.Hs.eg.db)


enrich1 <- enrichGO(gene         = gene.df1$ENSEMBL,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'ENSEMBL',
                 ont           = "ALL",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)

enrich2 <- enrichGO(gene         = gene.df2$ENSEMBL,
                    OrgDb         = org.Hs.eg.db,
                    keyType       = 'ENSEMBL',
                    ont           = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05)


enrich3 <- enrichGO(gene         = gene.df3$ENSEMBL,
                    OrgDb         = org.Hs.eg.db,
                    keyType       = 'ENSEMBL',
                    ont           = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05)


dim(enrich1);dim(enrich2);dim(enrich3)
head(enrich1,100)
head(enrich2,100)

# x <- as.data.frame(ego2)
# dim(x)




sorted_gene_vector1 <- sort(gene_vector1, decreasing = TRUE)
sorted_gene_vector2 <- sort(gene_vector2, decreasing = TRUE)
sorted_gene_vector3 <- sort(gene_vector3, decreasing = TRUE)






gseGO1 <- gseGO(geneList     = sorted_gene_vector1,
              OrgDb        = org.Hs.eg.db,
              ont          = "ALL",
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)

gseGO2 <- gseGO(geneList     = sorted_gene_vector2,
                OrgDb        = org.Hs.eg.db,
                ont          = "ALL",
                minGSSize    = 100,
                maxGSSize    = 500,
                pvalueCutoff = 0.05,
                verbose      = FALSE)


dim(gseGO1);dim(gseGO2)
head(ego3,100)
x <- as.data.frame(ego3)
dim(x)

gseGO3 <- gseGO(geneList     = sorted_gene_vector3,
                OrgDb        = org.Hs.eg.db,
                ont          = "ALL",
                minGSSize    = 100,
                maxGSSize    = 500,
                pvalueCutoff = 0.05,
                verbose      = FALSE)


kk1 <- enrichKEGG(gene         = as.character(all_regulated_GSE76250$Gene.ID),
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)


kk2 <- enrichKEGG(gene         = as.character(all_regulated_GSE38959$Gene.ID),
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)

kk_m1 <- enrichKEGG(gene         = as.character(meta_GSE76250_GSE38959$Gene.ID),
                  organism     = 'hsa',
                  pvalueCutoff = 0.05)

head(kk1,100);head(kk2,100);head(kk_m1)



kk3 <- gseKEGG(geneList     = sorted_gene_vector1,
               organism     = 'hsa',
               minGSSize    = 120,
               pvalueCutoff = 0.05,
               verbose      = FALSE)

kk4 <- gseKEGG(geneList     = sorted_gene_vector2,
               organism     = 'hsa',
               minGSSize    = 120,
               pvalueCutoff = 0.05,
               verbose      = FALSE)


kkm1 <- gseKEGG(geneList     = sorted_gene_vector3,
               organism     = 'hsa',
               minGSSize    = 120,
               pvalueCutoff = 0.05,
               verbose      = FALSE)

head(kk3);head(kk4)



mkk1 <- enrichMKEGG(gene = as.character(all_regulated_GSE76250$Gene.ID),
                   organism = 'hsa',
                   pvalueCutoff = 1,
                   qvalueCutoff = 1)
mkk2 <- enrichMKEGG(gene = as.character(all_regulated_GSE38959$Gene.ID),
                   organism = 'hsa',
                   pvalueCutoff = 1,
                   qvalueCutoff = 1)
mkk3 <- enrichMKEGG(gene = as.character(meta_GSE76250_GSE38959$Gene.ID),
                    organism = 'hsa',
                    pvalueCutoff = 1,
                    qvalueCutoff = 1)


dim(mkk1);dim(mkk2)
head(mkk1);head(mkk2)                   



mkk3 <- gseMKEGG(geneList = sorted_gene_vector1,
                 organism = 'hsa',
                 pvalueCutoff = 1)

mkk4 <- gseMKEGG(geneList = sorted_gene_vector2,
                 organism = 'hsa',
                 pvalueCutoff = 1)


mkk5 <- gseMKEGG(geneList = sorted_gene_vector3,
                 organism = 'hsa',
                 pvalueCutoff = 1)


head(mkk3);head(mkk4)


browseKEGG(kk, 'hsa04110')



library("pathview")
hsa04110 <- pathview(gene.data  = sorted_gene_vector,
                     pathway.id = "hsa04110",
                     species    = "hsa",
                     limit      = list(gene=max(abs(sorted_gene_vector)), cpd=1))



w1 <- enrichWP(as.character(all_regulated_GSE38959$Gene.ID), organism = "Homo sapiens") 
w2 <- enrichWP(as.character(all_regulated_GSE76250$Gene.ID), organism = "Homo sapiens") 
w3 <- enrichWP(as.character(meta_GSE76250_GSE38959$Gene.ID), organism = "Homo sapiens") 

head(w1);head(w2)



w2 <- gseWP(sorted_gene_vector, organism = "Homo sapiens")
head(w2)

library(ReactomePA)

x <- enrichPathway(gene=all_regulated_GSE76250$Gene.ID, pvalueCutoff = 0.05, readable=TRUE)
head(x)


y <- gsePathway(sorted_gene_vector, 
                pvalueCutoff = 0.2,
                pAdjustMethod = "BH", 
                verbose = FALSE)
head(y)


viewPathway("E2F mediated regulation of DNA replication", 
            readable = TRUE, 
            foldChange = sorted_gene_vector)



x1 <- enrichDO(gene          = all_regulated_GSE76250$Gene.ID,
              ont           = "DO",
              pvalueCutoff  = 0.05,
              pAdjustMethod = "BH",
              universe      = names(sorted_gene_vector1),
              minGSSize     = 5,
              maxGSSize     = 500,
              qvalueCutoff  = 0.01,
              readable      = FALSE)

x2 <- enrichDO(gene          = all_regulated_GSE38959$Gene.ID,
              ont           = "DO",
              pvalueCutoff  = 0.05,
              pAdjustMethod = "BH",
              universe      = names(sorted_gene_vector2),
              minGSSize     = 5,
              maxGSSize     = 500,
              qvalueCutoff  = 0.01,
              readable      = FALSE)

x3 <- enrichDO(gene          = meta_GSE76250_GSE38959$Gene.ID,
               ont           = "DO",
               pvalueCutoff  = 0.05,
               pAdjustMethod = "BH",
               universe      = names(sorted_gene_vector2),
               minGSSize     = 5,
               maxGSSize     = 500,
               qvalueCutoff  = 0.01,
               readable      = FALSE)

head(x1,16)
dim(x1)


ncg1 <- enrichNCG(all_regulated_GSE76250$Gene.ID) 
ncg2 <- enrichNCG(all_regulated_GSE38959$Gene.ID) 
ncg3 <- enrichNCG(meta_GSE76250_GSE38959$Gene.ID) 
head(ncg2)


dgn1 <- enrichDGN(all_regulated_GSE76250$Gene.ID) 
dgn2 <- enrichDGN(all_regulated_GSE38959$Gene.ID) 
dgn3 <- enrichDGN(meta_GSE76250_GSE38959$Gene.ID) 

head(dgn1)

y1 <- gseDO(sorted_gene_vector1,
           minGSSize     = 120,
           pvalueCutoff  = 0.2,
           pAdjustMethod = "BH",
           verbose       = FALSE)


y2 <- gseDO(sorted_gene_vector2,
            minGSSize     = 120,
            pvalueCutoff  = 0.2,
            pAdjustMethod = "BH",
            verbose       = FALSE)

y3 <- gseDO(sorted_gene_vector3,
            minGSSize     = 120,
            pvalueCutoff  = 0.2,
            pAdjustMethod = "BH",
            verbose       = FALSE)

head(y, 3)


ncg1 <- gseNCG(sorted_gene_vector1,
              pvalueCutoff  = 0.5,
              pAdjustMethod = "BH",
              verbose       = FALSE)
ncg2 <- gseNCG(sorted_gene_vector2,
               pvalueCutoff  = 0.5,
               pAdjustMethod = "BH",
               verbose       = FALSE)

ncg3 <- gseNCG(sorted_gene_vector3,
               pvalueCutoff  = 0.5,
               pAdjustMethod = "BH",
               verbose       = FALSE)

ncg <- setReadable(ncg, 'org.Hs.eg.db')
head(ncg, 3) 



dgn1 <- gseDGN(sorted_gene_vector1,
              pvalueCutoff  = 0.2,
              pAdjustMethod = "BH",
              verbose       = FALSE)
dgn <- setReadable(dgn, 'org.Hs.eg.db')
head(dgn1, 3) 

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



