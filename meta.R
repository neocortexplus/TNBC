
library(GEOquery)
library(limma)
library(umap)
library(tidyr)
library(dplyr)
library(ggplot2)
library(Biobase)
library(readxl)
library(sva)

###### GSE38959

# 47 samples 
################################################################
#   Differential expression analysis with limma

library(GEOquery)
library(limma)
library(umap)
library(tidyr)
library(dplyr)
library(ggplot2)
library(Biobase)
library(readxl)
# load series and platform data from GEO
# 
# gset <- tryCatch({
#   getGEO("GSE38959", GSEMatrix = TRUE, AnnotGPL = TRUE)
# }, error = function(e) e)


gset <- getGEO("GSE38959", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL4133", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- "11111111111111111111111111111100000000000000000"



sml <- strsplit(gsms, split="")[[1]]

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

ex <- exprs(gset)


ex_df <- as.data.frame(ex)
feature_data <- fData(gset)
pheno_data <- pData(phenoData(gset))

# Assuming 'feature_data' is your dataframe and 'gene_assignment' is the column with the data

feature_data <- feature_data %>% rename(Gene.Symbol = Gene.symbol)

count_missing <- function(x) {
  sum(is.na(x) | x %in% c("--", "nan", "NA", "", " ", "---"))
}


# Apply this function to each of the specified columns

missing_counts <- sapply(feature_data[c("Gene.ID", "Gene.Symbol")], count_missing)

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

GSE38959 <- new_gset

max(exprs(GSE38959))
min(exprs(GSE38959))

save(GSE38959, file = "/home/aiusrdata/RCode/TNBC/GSE38959.RData")
load("/home/aiusrdata/RCode/TNBC/GSE38959.RData")

summary(GSE38959)  # Gives a summary of the ExpressionSet
str(GSE38959)      # Shows the structure of the ExpressionSet




########## GSE76250


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

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

ex <- exprs(gset)



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
feature_data <- feature_data[, c("ID", "Gene.Symbol")]

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



fs_data <- aggregated_data_with_ids[, c("ID", "Gene.Symbol")]
row.names(aggregated_data_with_ids) <- aggregated_data_with_ids$Gene.Symbol

aggregated_data_with_ids <- aggregated_data_with_ids %>%
  select(-ID, -Gene.Symbol)



expression_matrix <- as.matrix(aggregated_data_with_ids)

rownames(fs_data) <- fs_data$Gene.Symbol


new_gset <- ExpressionSet(
  assayData = expression_matrix,       # Replace 'expression_matrix' with your actual matrix of expression data
  phenoData = phenoData(gset),         # Use the 'phenoData' directly if it is already an AnnotatedDataFrame
  featureData = AnnotatedDataFrame(fs_data)  # Replace 'fs_data' with your feature data
)

# Inspecting the new ExpressionSet

GSE76250 <- new_gset

max(exprs(GSE76250))
min(exprs(GSE76250))

save(GSE76250, file = "/home/aiusrdata/RCode/TNBC/GSE76250.RData")
load("/home/aiusrdata/RCode/TNBC/GSE76250.RData")

summary(GSE76250)  # Gives a summary of the ExpressionSet
str(GSE76250)      # Shows the structure of the ExpressionSet




######### meta 

load("/home/aiusrdata/RCode/TNBC/GSE38959.RData")
load("/home/aiusrdata/RCode/TNBC/GSE76250.RData")

GSE38959_ex <- exprs(GSE38959)
GSE76250_ex <- exprs(GSE76250)

GSE38959_fs <- fData(GSE38959)
GSE76250_fs <- fData(GSE76250)


GSE38959_pd <- pData(phenoData(GSE38959))
GSE76250_pd <- pData(phenoData(GSE76250))


common_genes <- intersect(rownames(GSE38959_ex), rownames(GSE76250_ex))

GSE38959_ex_common <- GSE38959_ex[common_genes, ]
GSE76250_ex_common <- GSE76250_ex[common_genes, ]

all <- cbind(GSE38959_ex_common, GSE76250_ex_common)

all_fs <- GSE38959_fs[common_genes,]

dim(all_fs);dim(GSE38959_ex_common) ; dim(GSE76250_ex_common);dim(all )


GSE38959_pd <- GSE38959_pd %>%
  rename(
    age = `age (y):ch1`,  # Correctly reference the age column if that's its name
    disease = `disease state:ch1`,
    gender = `gender:ch1`
  )

# Adjustments for GSE76250_pd
GSE76250_pd <- GSE76250_pd %>%
  rename(
    age = `age (yrs):ch1`,  # Ensure this is the correct column name in your dataframe
    disease = `tissue type:ch1`,
    gender = `gender:ch1`
  )


all_columns <- union(names(GSE38959_pd), names(GSE76250_pd))

# Add missing columns as NA
missing_columns_38959 <- setdiff(all_columns, names(GSE38959_pd))
GSE38959_full <- GSE38959_pd
for (col in missing_columns_38959) {
  GSE38959_full[[col]] <- NA
}

# Ensure the columns are in the same order
GSE38959_full <- GSE38959_full %>% select(all_of(all_columns))

# Repeat for GSE76250
missing_columns_76250 <- setdiff(all_columns, names(GSE76250_pd))
GSE76250_full <- GSE76250_pd
for (col in missing_columns_76250) {
  GSE76250_full[[col]] <- NA
}

# Ensure the columns are in the same order
GSE76250_full <- GSE76250_full %>% select(all_of(all_columns))

# Now rbind them together
all_pd <- rbind(GSE38959_full, GSE76250_full)


all_gset <- ExpressionSet(
  assayData = as.matrix(all),
  phenoData = AnnotatedDataFrame(all_pd),
  featureData = AnnotatedDataFrame(all_fs)  # Assuming all_fs is already prepared
)


ex <- exprs(all_gset)


t_all <- t(all)

dim(t_all)

batch <- c(rep(0, ncol(GSE38959_ex_common)), rep(1, ncol(GSE76250_ex_common)))
batch <-   factor(batch)

ex <- na.omit(ex)
ex <- ex[!duplicated(ex), ]


# Perform UMAP dimensionality reduction
ump <- umap(t_all, n_neighbors = 5, random_state = 123)

# Plotting
par(mar=c(3,3,2,6), xpd=TRUE)
plot(ump$layout, main="UMAP plot, nbrs=5", xlab="", ylab="", col=batch, pch=20, cex=1.5)

# Add legend
legend("topright", inset=c(-0.15,0), legend=c("GSE38959", "GSE76250"), pch=20,
       col=1:2, title="Group", pt.cex=1.5)

allc<- ComBat(all, batch)
dim(allc)


# Perform UMAP dimensionality reduction
ump <- umap(t(allc), n_neighbors = 5, random_state = 123)

# Plotting
par(mar=c(3,3,2,6), xpd=TRUE)
plot(ump$layout, main="UMAP plot, nbrs=5", xlab="", ylab="", col=batch, pch=20, cex=1.5)

# Add legend
legend("topright", inset=c(-0.15,0), legend=c("GSE38959", "GSE76250"), pch=20,
       col=1:2, title="Group", pt.cex=1.5)



