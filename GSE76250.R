# Version info: R 4.2.2, Biobase 2.58.0, GEOquery 2.66.0, limma 3.54.0
################################################################
#   Differential expression analysis with limma
library(GEOquery)
library(limma)
library(umap)
library(tidyr)
library(dplyr)
library(ggplot2)
library(Biobase)
library(ggrepel)  # Load ggrepel for better label management

gset <- getGEO("GSE76250", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL17586", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples


gsms <- paste0("00000000000000000000000000000000000000000000000000",
               "00000000000000000000000000000000000000000000000000",
               "00000000000000000000000000000000000000000000000000",
               "000000000000000111111111111111111111111111111111")
sml <- strsplit(gsms, split="")[[1]]




sml <- strsplit(gsms, split="")[[1]]

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

ex <- exprs(gset)

colnames(ex)

ex_df <- as.data.frame(ex)
feature_data <- fData(gset)
pheno_data <- pData(phenoData(gset))

# Assuming 'feature_data' is your dataframe and 'gene_assignment' is the column with the data


feature_data$Gene.Symbol <- sapply(strsplit(as.character(feature_data$gene_assignment), " // "), function(x) x[2])
feature_data$Gene.ID <- sapply(strsplit(as.character(feature_data$gene_assignment), " // "), function(x) x[5])

feature_data$Gene.ID1 <- sapply(strsplit(as.character(feature_data$gene_assignment), " // "), function(x) x[1])
feature_data$Gene.Symbol1 <- sapply(strsplit(as.character(feature_data$gene_assignment), " // "), function(x) x[2])

feature_data$Gene.ID <- sapply(strsplit(as.character(feature_data$Gene.ID), " /// "), function(x) x[1])


count_missing <- function(x) {
  sum(is.na(x) | x %in% c("--", "nan", "NA", "", " ", "---"))
}

# Apply this function to each of the specified columns

missing_counts <- sapply(feature_data[c("Gene.ID", "Gene.Symbol", "Gene.ID1", "Gene.Symbol1")], count_missing)

print(missing_counts)


feature_data <- feature_data %>%
  filter(!is.na(Gene.Symbol) & !Gene.Symbol %in% c("--", "nan", "NA", "", " ", "---"))


# feature_data <- feature_data[, c("ID", "Gene.ID", "Gene.symbol", "Platform_SPOTID")]
feature_data <- feature_data[, c("ID","Gene.ID","Gene.Symbol")]

combined_data <- merge(ex_df, feature_data, by.x = "row.names", by.y = "ID")
# Rename the 'Row.names' column to 'ID' in the resulting data frame for clarity
names(combined_data)[1] <- "ID"




aggregated_data <- combined_data %>%
  group_by(Gene.Symbol) %>%
  summarise(across(starts_with("GSM"), mean, na.rm = TRUE))


aggregated_data <- as.data.frame(aggregated_data)


feature_data <- feature_data[!duplicated(feature_data$Gene.Symbol), ]

aggregated_data_with_ids <- aggregated_data %>%
  left_join(feature_data, by = "Gene.Symbol")

# Check the first few rows of the new dataframe to confirm the merge
head(aggregated_data_with_ids)



fs_data <- aggregated_data_with_ids[, c("ID", "Gene.Symbol","Gene.ID")]
row.names(aggregated_data_with_ids) <- aggregated_data_with_ids$Gene.Symbol

aggregated_data_with_ids <- aggregated_data_with_ids %>%
  select(-ID, -Gene.Symbol)
aggregated_data_with_ids <- aggregated_data_with_ids %>%
  select(-Gene.ID)



expression_matrix <- as.matrix(aggregated_data_with_ids)

rownames(fs_data) <- fs_data$Gene.Symbol


new_gset <- ExpressionSet(
  assayData = expression_matrix,       # Replace 'expression_matrix' with your actual matrix of expression data
  phenoData = phenoData(gset),         # Use the 'phenoData' directly if it is already an AnnotatedDataFrame
  featureData = AnnotatedDataFrame(fs_data)  # Replace 'fs_data' with your feature data
)



# Inspecting the new ExpressionSet

gset <- new_gset





# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("TNBC","normal"))
levels(gs) <- groups
gset$group <- gs






design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)





gset <- gset[complete.cases(exprs(gset)), ] # skip missing values

fit <- lmFit(gset, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- paste(groups[1], groups[2], sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)

tT <- subset(tT, select=c("ID" ,"Gene.Symbol","adj.P.Val","P.Value","t","B","logFC"))
write.table(tT, file=stdout(), row.names=F, sep="\t")





# Visualize and quality control test results.
# Build histogram of P-values for all genes. Normal test
# assumption is that most genes are not differentially expressed.

tT2 <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)

save(tT2, file = "/home/aiusrdata/RCode/TNBC/results/GSE76250_tT2.RData")

hist(tT2$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
     ylab = "Number of genes", main = "P-adj value distribution")


upregulated_all <- tT2[tT2$logFC >= 1.5 & tT2$adj.P.Val <= 0.01, ]
downregulated_all <- tT2[tT2$logFC <= -1.5 & tT2$adj.P.Val <= 0.01, ]

dim(upregulated_all);dim(downregulated_all)

p_value_threshold <- 0.01
tmp1_tT2 <- tT2[abs(tT2$logFC) >= 1 & tT2$adj.P.Val < p_value_threshold, ]
upregulated_tmp1_tT2 <- tT2[tT2$logFC >= 1 & tT2$adj.P.Val < p_value_threshold, ]
downregulated_tmp1_tT2 <- tT2[tT2$logFC <= -1 & tT2$adj.P.Val < p_value_threshold, ]
dim(tmp1_tT2 );dim(upregulated_tmp1_tT2);dim(downregulated_tmp1_tT2 )

kk_temp1 <- enrichKEGG(gene         = as.character(tmp1_tT2 $Gene.ID),
                       organism     = 'hsa',
                       pvalueCutoff = 0.05)



head(kk_temp1@result$Description,10)




# summarize test results as "up", "down" or "not expressed"
dT <- decideTests(fit2, adjust.method="fdr", p.value=0.05, lfc=0)
dt_df <- as.data.frame(dT)
dim(dt_df)


upregulated_genes <- dt_df[dt_df$`normal-TNBC` == 1, , drop = FALSE]

downregulated_genes <- dt_df[dt_df$`normal-TNBC` == -1, , drop = FALSE]

x <- tT2[rownames(upregulated_genes),]

dim(upregulated_genes);dim(downregulated_genes)



# Venn diagram of results
vennDiagram(dT, circle.col=palette())

# create Q-Q plot for t-statistic
t.good <- which(!is.na(fit2$F)) # filter out bad probes
qqt(fit2$t[t.good], fit2$df.total[t.good], main="Moderated t statistic")

# volcano plot (log P-value vs log fold change)
colnames(fit2) # list contrast names
ct <- 1        # choose contrast of interest
# Please note that the code provided to generate graphs serves as a guidance to
# the users. It does not replicate the exact GEO2R web display due to multitude
# of graphical options.
# 
# The following will produce basic volcano plot using limma function:
volcanoplot(fit2, coef=ct, main=colnames(fit2)[ct], pch=20,
            highlight=length(which(dT[,ct]!=0)), names=rep('+', nrow(fit2)))



volcano_data <- data.frame(
  Log2FoldChange = fit2$coefficients[, ct],
  NegLog10PValue = -log10(fit2$p.value[, ct]),
  GeneStatus = factor(dT[, ct], levels = c(-1, 0, 1), labels = c("down", "not significant", "up"))
)

ggplot(volcano_data, aes(x = Log2FoldChange, y = NegLog10PValue, color = GeneStatus)) +
  geom_point(alpha = 1, size = 1) +  # Adjust point size and opacity
  scale_color_manual(values = c("down" = "blue", "up" = "red"),
                     labels = c("down", "up")) +  # Assign colors and labels to legend
  labs(title = "Volcano Plot: Normal vs TNBC",
       x = "Log2 Fold Change",
       y = "-Log10 P-value",
       color = "Padj < 0.05") +  # Title for the color legend
  theme_minimal() +
  theme(legend.title = element_text(size = 10),  # Adjust legend title text size
        legend.position = "bottom",  # Position the legend at the bottom
        legend.justification = "left",  # Justify legend to the left
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 10),  # Adjust text size within the legend
        axis.text = element_text(color = "black"),  # Ensure axis labels are visible
        axis.title = element_text(color = "black"))  # Ensure axis titles are visible


# Ensure the volcano_data dataframe is defined as previously described
volcano_data <- data.frame(
  Log2FoldChange = fit2$coefficients[, ct],
  NegLog10PValue = -log10(fit2$p.value[, ct]),
  GeneStatus = factor(dT[, ct], levels = c(-1, 0, 1), labels = c("down", "not significant", "up"))
)

# Create the 'Status' column based on the significance criteria
volcano_data <- volcano_data %>%
  mutate(Status = case_when(
    NegLog10PValue > -log10(0.05) & Log2FoldChange > 2  ~ "Upregulated Significant",
    NegLog10PValue > -log10(0.05) & Log2FoldChange < -2 ~ "Downregulated Significant",
    Log2FoldChange > 0 ~ "Upregulated Not Significant",
    TRUE ~ "Downregulated Not Significant"
  ))

ggplot(volcano_data, aes(x = Log2FoldChange, y = NegLog10PValue, color = Status)) +
  geom_point(alpha = 1, size = 1) +
  scale_color_manual(values = c(
    "Upregulated Significant" = "brown",
    "Downregulated Significant" = "green",
    "Upregulated Not Significant" = "red",
    "Downregulated Not Significant" = "blue"
  )) +
  labs(title = "Volcano Plot: Normal vs TNBC",
       x = "Log2 Fold Change",
       y = "-Log10 P-value",
       color = "Gene Status") +
  theme_minimal() +
  theme(legend.title = element_text(size = 10),
        legend.position = "bottom",  # Position the legend at the bottom
        legend.justification = "left",  # Justify legend to the left
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 10),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"))




# Ensure the volcano_data dataframe is defined as previously described
volcano_data <- data.frame(
  Log2FoldChange = fit2$coefficients[, ct],
  NegLog10PValue = -log10(fit2$p.value[, ct]),
  GeneStatus = factor(dT[, ct], levels = c(-1, 0, 1), labels = c("down", "not significant", "up"))
)

# Create the 'Status' column based on the significance criteria including a new "Non-significant" category for p-values > 0.05
volcano_data <- volcano_data %>%
  mutate(Status = case_when(
    NegLog10PValue <= -log10(0.05) ~ "Non-significant",
    NegLog10PValue > -log10(0.05) & Log2FoldChange > 2  ~ "Upregulated Significant",
    NegLog10PValue > -log10(0.05) & Log2FoldChange < -2 ~ "Downregulated Significant",
    Log2FoldChange >= 0 ~ "Upregulated Not Significant",
    TRUE ~ "Downregulated Not Significant"
  ))
volcano_data$GeneSymbol <- row.names(volcano_data)


# Plot with the updated color scheme
ggplot(volcano_data, aes(x = Log2FoldChange, y = NegLog10PValue, color = Status)) +
  geom_point(alpha = 1, size = 1) +
  scale_color_manual(values = c(
    "Non-significant" = "grey",
    "Upregulated Significant" = "brown",
    "Downregulated Significant" = "green",
    "Upregulated Not Significant" = "red",
    "Downregulated Not Significant" = "blue"
  )) +
  labs(title = "Volcano Plot: Normal vs TNBC",
       x = "Log2 Fold Change",
       y = "-Log10 P-value",
       color = "Gene Status") +
  theme_minimal() +
  theme(legend.title = element_text(size = 10),
        legend.position = "bottom",  # Position the legend at the bottom
        legend.justification = "left",  # Justify legend to the left
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 10),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"))


# Load necessary library
library(ggplot2)
library(dplyr)

# Assuming volcano_data is already loaded and structured as described

# Define criteria for labeling
label_data <- volcano_data %>%
  filter(Status %in% c("Upregulated Significant", "Downregulated Significant"))

# Plot with labels for significant points
ggplot(volcano_data, aes(x = Log2FoldChange, y = NegLog10PValue, color = Status)) +
  geom_point(alpha = 1, size = 1) +  # Basic point layer
  geom_text(data = label_data, aes(label = GeneSymbol), vjust = 1.5, hjust = 0.5, 
            check_overlap = TRUE, size = 3, position = position_jitter(width = 0.2, height = 0)) +  # Add jitter to labels for clarity
  scale_color_manual(values = c(
    "Non-significant" = "grey",
    "Upregulated Significant" = "brown",
    "Downregulated Significant" = "green",
    "Upregulated Not Significant" = "red",
    "Downregulated Not Significant" = "blue"
  )) +
  labs(title = "Volcano Plot: Normal vs TNBC",
       x = "Log2 Fold Change",
       y = "-Log10 P-value",
       color = "Gene Status") +
  theme_minimal() +
  theme(legend.title = element_text(size = 10),
        legend.position = "bottom",
        legend.justification = "left",
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 10),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"))



# all labels

ggplot(volcano_data, aes(x = Log2FoldChange, y = NegLog10PValue, color = Status)) +
  geom_point(alpha = 1, size = 1) +  # Basic point layer
  geom_text_repel(data = label_data, aes(label = GeneSymbol),
                  box.padding = 0.35, point.padding = 0.5, 
                  max.overlaps = Inf,  # Allow infinite attempts to place labels
                  size = 3) +  # Consider reducing size if labels still overlap
  scale_color_manual(values = c(
    "Non-significant" = "grey",
    "Upregulated Significant" = "brown",
    "Downregulated Significant" = "green",
    "Upregulated Not Significant" = "red",
    "Downregulated Not Significant" = "blue"
  )) +
  labs(title = "Volcano Plot: Normal vs TNBC",
       x = "Log2 Fold Change",
       y = "-Log10 P-value",
       color = "Gene Status") +
  theme_minimal() +
  theme(legend.title = element_text(size = 10),
        legend.position = "bottom",
        legend.justification = "left",
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 10),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"),
        panel.spacing = unit(1, "lines"))  # Increase spacing around the plot area

#ggsave("volcano_plot.png", width = 10, height = 8, dpi = 300)




# MD plot (log fold change vs mean log expression)
# highlight statistically significant (p-adj < 0.05) probes
plotMD(fit2, column=ct, status=dT[,ct], legend=F, pch=20, cex=1)
abline(h=0)

logFC <- fit2$coefficients[, "TNBC-normal"]  # Extract log fold changes for the "normal-TNBC" contrast
aveExpr <- fit2$Amean                      # Extract average expressions
geneNames <- fit2$genes$Gene.Symbol        # Extract gene symbols

# Create a data frame for ggplot
data <- data.frame(
  AveExpr = aveExpr,
  LogFC = logFC,
  GeneNames = geneNames,
  HighImpact = abs(logFC) >= 2,
  LogFCSign = ifelse(logFC > 0, "Positive", "Negative")  # Create a new column for logFC sign
)

# Plot using ggplot2 and ggrepel
p <- ggplot(data, aes(x = AveExpr, y = LogFC, label = GeneNames)) +
  geom_point(aes(color = LogFCSign), size = 1.5) +  # Color points based on logFC sign
  scale_color_manual(values = c("Positive" = "red", "Negative" = "blue")) +  # Assign colors
  geom_text_repel(data = subset(data, HighImpact),  # Only label high impact genes
                  aes(label = GeneNames), 
                  box.padding = 0.35, 
                  point.padding = 0.5,
                  size = 3,
                  max.overlaps = Inf) +
  theme_minimal() +
  labs(title = "Expression Plot", x = "Average Expression", y = "Log Fold Change") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(legend.position = "none")  # Remove the legend

# Print the plot
print(p)

# Create a new column for logFC sign with updated conditions
data$LogFCSign <- ifelse(logFC > 0 & logFC < 2, "Red",
                         ifelse(logFC < 0 & logFC > -2, "Blue",
                                ifelse(logFC >= 2, "Brown", "Green")))

# Update the plot using ggplot2 and ggrepel
p <- ggplot(data, aes(x = AveExpr, y = LogFC, label = GeneNames)) +
  geom_point(aes(color = LogFCSign), size = 1.5) +  # Color points based on the new logFC sign
  scale_color_manual(values = c("Red" = "red", "Blue" = "blue", "Brown" = "brown", "Green" = "green")) +  # Assign colors
  geom_text_repel(data = subset(data, HighImpact),  # Only label high impact genes
                  aes(label = GeneNames), 
                  box.padding = 0.35, 
                  point.padding = 0.5,
                  size = 3,
                  max.overlaps = Inf) +
  theme_minimal() +
  labs(title = "Expression Plot", x = "Average Expression", y = "Log Fold Change") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(legend.position = "none")  # Remove the legend

# Print the plot
print(p)


################################################################
# General expression data analysis
ex <- exprs(gset)

# box-and-whisker plot
dev.new(width=3+ncol(gset)/6, height=5)
ord <- order(gs)  # order samples by group
palette(c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02",
          "#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666"))
par(mar=c(7,4,2,1))
title <- paste ("GSE76250", "/", annotation(gset), sep ="")
boxplot(ex[,ord], boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=gs[ord])
legend("topleft", groups, fill=palette(), bty="n")
dev.off()

# expression value distribution
par(mar=c(4,4,2,1))
title <- paste ("GSE76250", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex, group=gs, main=title, legend ="topright")

# UMAP plot (dimensionality reduction)
ex <- na.omit(ex) # eliminate rows with NAs
ex <- ex[!duplicated(ex), ]  # remove duplicates
ump <- umap(t(ex), n_neighbors = 15, random_state = 123)
par(mar=c(3,3,2,6), xpd=TRUE)
plot(ump$layout, main="UMAP plot, nbrs=15", xlab="", ylab="", col=gs, pch=10, cex=1.5)
legend("topright", inset=c(-0.15,0), legend=levels(gs), pch=20,
       col=1:nlevels(gs), title="Group", pt.cex=1.5)
library("maptools")  # point labels without overlaps
pointLabel(ump$layout, labels = rownames(ump$layout), method="SANN", cex=0.6)

# mean-variance trend, helps to see if precision weights are needed
plotSA(fit2, main="Mean variance trend, GSE76250")




gene_data <- as.data.frame(t(ex))
pheno_data <- pData(phenoData(gset))

complete <- merge(gene_data, pheno_data, by.x = "row.names", by.y = "row.names")

data_to_plot <- complete %>%
  select(TOP2A, group)  

# Create the plot
ggplot(data_to_plot, aes(x = group, y = TOP2A, color = group)) +
  geom_boxplot() +
  labs(title = "Expression of TOP2A in TNBC and Normal Samples",
       x = "Sample Group",
       y = "Expression Level") +
  theme_minimal()  # Adds a minimalistic theme to the plot

# Print the plot
print(ggplot_object)




