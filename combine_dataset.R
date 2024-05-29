
################################################################
#   Differential expression analysis with limma
library(GEOquery)
library(limma)
library(umap)
library(tidyr)
library(dplyr)
library(ggplot2)
library(Biobase)

# load series and platform data from GEO
# 
# gset <- tryCatch({
#   getGEO("GSE76250", GSEMatrix = TRUE, AnnotGPL = TRUE)
# }, error = function(e) e)


gset <- getGEO("GSE102088", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL17586", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples

gsms <- paste0(rep("0", ncol(exprs(gset))), collapse = "")

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

feature_data$Gene.ID2 <- sapply(strsplit(as.character(feature_data$gene_assignment), " /// "), function(x) strsplit(x[1], " // ")[[1]][1])
feature_data$Gene.Symbol2 <- sapply(strsplit(as.character(feature_data$gene_assignment), " /// "), function(x) strsplit(x[1], " // ")[[1]][2])



feature_data$Gene.Symbol <- sapply(strsplit(as.character(feature_data$gene_assignment), " // "), function(x) x[2])
feature_data$Gene.ID <- sapply(strsplit(as.character(feature_data$gene_assignment), " // "), function(x) x[5])

feature_data$Gene.ID1 <- sapply(strsplit(as.character(feature_data$gene_assignment), " // "), function(x) x[1])
feature_data$Gene.Symbol1 <- sapply(strsplit(as.character(feature_data$gene_assignment), " // "), function(x) x[2])

feature_data$Gene.ID <- sapply(strsplit(as.character(feature_data$Gene.ID), " /// "), function(x) x[1])



count_missing <- function(x) {
  sum(is.na(x) | x %in% c("--", "nan", "NA", "", " ", "---"))
}

missing_counts <- sapply(feature_data[c("Gene.ID","Gene.ID1","Gene.ID2", "Gene.Symbol", "Gene.Symbol1","Gene.Symbol2")], count_missing)

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

gset <- new_gset

GSE102088 <- gset

save(GSE102088, file = "/home/aiusrdata/RCode/TNBC/GSE102088.RData")

######################## GSE76124

gset <- getGEO("GSE76124", GSEMatrix =TRUE, getGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples

gsms <- paste0(rep("1", ncol(exprs(gset))), collapse = "")

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


feature_data$Gene.Symbol <- sub(" ///.*", "", feature_data$Gene.Symbol)
feature_data$Gene.ID <- sub(" ///.*", "", feature_data$ENTREZ_GENE_ID)


count_missing <- function(x) {
  sum(is.na(x) | x %in% c("--", "nan", "NA", "", " ", "---"))
}

missing_counts <- sapply(feature_data[c("Gene.ID", "Gene.Symbol", "ID")], count_missing)

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

gset <- new_gset


GSE76124 <- gset

save(GSE76124, file = "/home/aiusrdata/RCode/TNBC/GSE76124.RData")




############# combine 


load("/home/aiusrdata/RCode/TNBC/GSE102088.RData")
load("/home/aiusrdata/RCode/TNBC/GSE76124.RData")


GSE102088_ex <- exprs(GSE102088)
GSE76124_ex <- exprs(GSE76124)

GSE102088_fs <-fData(GSE102088)
GSE76124_fs <- fData(GSE76124)





common_features <- intersect(rownames(GSE102088_ex), rownames(GSE76124_ex))
all_df <- cbind(GSE102088_ex[common_features, ], GSE76124_ex[common_features, ])



all_fs <- subset(GSE102088_fs, Gene.Symbol %in% common_features)




GSE102088_pd <- pData(phenoData(GSE102088))
GSE76124_pd <- pData(phenoData(GSE76124))



GSE102088_pd <- GSE102088_pd %>%
  rename(
    age = `age:ch1`,  # Correctly reference the age column if that's its name
  )

GSE76124_pd <- GSE76124_pd %>%
  rename(
    age = `age (years):ch1`,  # Ensure this is the correct column name in your dataframe
    disease = `triple-negative status:ch1`,
    gender = `gender:ch1`
  )


all_columns <- union(names(GSE102088_pd), names(GSE76124_pd))

# Add missing columns as NA
missing_columns_GSE102088_pd <- setdiff(all_columns, names(GSE102088_pd))
GSE102088_full <- GSE102088_pd

for (col in missing_columns_GSE102088_pd) {
  GSE102088_full[[col]] <- NA
}

# Ensure the columns are in the same order
GSE102088_full <- GSE102088_full %>% select(all_of(all_columns))

# Repeat for GSE76250
missing_columns_GSE76124 <- setdiff(all_columns, names(GSE76124_pd))
GSE76124_full <- GSE76124_pd
for (col in missing_columns_GSE76124) {
  GSE76124_full[[col]] <- NA
}

# Ensure the columns are in the same order
GSE76124_full <- GSE76124_full %>% select(all_of(all_columns))

# Now rbind them together
all_pd <- rbind(GSE102088_full, GSE76124_full)


all_gset <- ExpressionSet(
  assayData = as.matrix(all_df),
  phenoData = AnnotatedDataFrame(all_pd),
  featureData = AnnotatedDataFrame(all_fs)  # Assuming all_fs is already prepared
)

