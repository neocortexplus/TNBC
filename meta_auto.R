# اصل اینه ها یادت باشه فردا

library(GEOquery)
library(limma)
library(umap)
library(tidyr)
library(dplyr)
library(ggplot2)
library(Biobase)
library(readxl)
library(sva)

#GSE65194 #########################

gset <- getGEO("GSE65194", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples

gsms <- paste0("11XX1111111X1111111XXXX1111XXXXXX111XXX1111X1111XX",
               "XXXX1XXX11X11XX111111111X1XXXXXXXX111X111111XXXXXX",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXX0X0XX00XX0XX0X00XX0X00XX")
sml <- strsplit(gsms, split="")[[1]]
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]



# log2 transformation
ex <- exprs(gset)
colnames(ex)
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

GSE65194 <- new_gset

max(exprs(GSE65194))
min(exprs(GSE65194))

save(GSE65194, file = "/home/aiusrdata/RCode/TNBC/GSE65194.RData")
load("/home/aiusrdata/RCode/TNBC/GSE65194.RData")

summary(GSE65194)  # Gives a summary of the ExpressionSet
str(GSE65194)      # Shows the structure of the ExpressionSet




#GSEGSE65212

#GSE31448  ####################

gset <- getGEO("GSE31448", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))




# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

ex <- exprs(gset)

gsms <- paste0(rep("1", ncol(ex)), collapse = "")

sml <- strsplit(gsms, split="")[[1]]


ex_df <- as.data.frame(ex)
feature_data <- fData(gset)
pheno_data <- pData(phenoData(gset))

feature_data$Gene.Symbol <- sub("///.*", "", feature_data$Gene.symbol)
# feature_data$ENTREZ_GENE_ID <- sub(" ///.*", "", feature_data$ENTREZ_GENE_ID)
feature_data$Gene.ID <- sub("///.*", "", feature_data$Gene.ID)

colnames(feature_data)

feature_data <- feature_data %>%
  rename(Gene.ID = ENTREZ_GENE_ID)

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

GSE31448 <- new_gset

max(exprs(GSE31448))
min(exprs(GSE31448))

save(GSE31448, file = "/home/aiusrdata/RCode/TNBC/GSE31448.RData")
load("/home/aiusrdata/RCode/TNBC/GSE31448.RData")

summary(GSE31448)  # Gives a summary of the ExpressionSet
str(GSE31448)      # Shows the structure of the ExpressionSet





#GSE45827 #########

gset <- getGEO("GSE45827", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples

gsms <- paste0("11X1XX11111XXXX1X11XXX1X11111111XX111XXX11XX1X11XX",
               "X11X11XX111111111XXXXXXXXXXXXXXXXXX00000000000XXXX",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXXX")
sml <- strsplit(gsms, split="")[[1]]

# filter out excluded samples (marked as "X")
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]




# log2 transformation
ex <- exprs(gset)
colnames(ex)
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

feature_data$Gene.Symbol <- sub("///.*", "", feature_data$Gene.Symbol)
feature_data$Gene.ID <- sub("///.*", "", feature_data$ID)



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

GSE45827 <- new_gset

max(exprs(GSE45827))
min(exprs(GSE45827))

save(GSE45827, file = "/home/aiusrdata/RCode/TNBC/GSE45827.RData")
load("/home/aiusrdata/RCode/TNBC/GSE45827")

summary(GSE45827)  # Gives a summary of the ExpressionSet
str(GSE45827)      # Shows the structure of the ExpressionSet




#GSE61724 ####


gset <- getGEO("GSE61724", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL6244", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples

gsms <- "XXXXXXXXX0X1XXXXXXXXXX111000XXX11XX1XX1XXXXXXX11X1XXX1X1XXXX1XX11XXX"
sml <- strsplit(gsms, split="")[[1]]

# filter out excluded samples (marked as "X")
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]






# log2 transformation
ex <- exprs(gset)
colnames(ex)
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

feature_data$Gene.Symbol <- sub("///.*", "", feature_data$Gene.Symbol)
feature_data$Gene.ID <- sub("///.*", "", feature_data$ID)



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

GSE61724 <- new_gset

max(exprs(GSE61724))
min(exprs(GSE61724))

save(GSE61724, file = "/home/aiusrdata/RCode/TNBC/GSE61724.RData")
load("/home/aiusrdata/RCode/TNBC/GSE61724.RData")

summary(GSE61724)  # Gives a summary of the ExpressionSet
str(GSE61724)      # Shows the structure of the ExpressionSet


#GSE36295  #####


gset <- getGEO("GSE36295", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL6244", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples

gsms <- "11111XXX000XXX0X0XXXXX00X0XXXXXXX0XXXXXXXX00XXXXXX"
sml <- strsplit(gsms, split="")[[1]]

# filter out excluded samples (marked as "X")
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]






# log2 transformation
ex <- exprs(gset)
colnames(ex)
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

feature_data$Gene.Symbol <- sub("///.*", "", feature_data$Gene.Symbol)
feature_data$Gene.ID <- sub("///.*", "", feature_data$ID)



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

GSE36295 <- new_gset

max(exprs(GSE36295))
min(exprs(GSE36295))

save(GSE36295, file = "/home/aiusrdata/RCode/TNBC/GSE36295.RData")
load("/home/aiusrdata/RCode/TNBC/GSE36295.RData")

summary(GSE36295)  # Gives a summary of the ExpressionSet
str(GSE36295)      # Shows the structure of the ExpressionSet



#GSE37751 ######


gset <- getGEO("GSE37751", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL6244", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples

gsms <- paste0("11111111111111111111111111111111111111111111111X0X",
               "XX0XXX0XX0XX000XXXXXXXXXXXXXXXXXX000XXXXXXX00X0XXX",
               "X0XXXXXX")
sml <- strsplit(gsms, split="")[[1]]

# filter out excluded samples (marked as "X")
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]






# log2 transformation
ex <- exprs(gset)
colnames(ex)
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

feature_data$Gene.Symbol <- sub("///.*", "", feature_data$Gene.Symbol)
feature_data$Gene.ID <- sub("///.*", "", feature_data$ID)



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

GSE37751 <- new_gset

max(exprs(GSE37751))
min(exprs(GSE37751))

save(GSE37751, file = "/home/aiusrdata/RCode/TNBC/GSE37751.RData")
load("/home/aiusrdata/RCode/TNBC/GSE37751.RData")

summary(GSE37751)  # Gives a summary of the ExpressionSet
str(GSE37751)      # Shows the structure of the ExpressionSet


#GSE38959  ###############################################################


# 47 samples 


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
colnames(ex)

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





#GSE76250 ###################################################################


gset <- getGEO("GSE76250", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL17586", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples

gsms <- paste0("11111111111111111111111111111111111111111111111111",
               "11111111111111111111111111111111111111111111111111",
               "11111111111111111111111111111111111111111111111111",
               "111111111111111000000000000000000000000000000000")



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

GSE76250 <- new_gset

max(exprs(GSE76250))
min(exprs(GSE76250))

save(GSE76250, file = "/home/aiusrdata/RCode/TNBC/GSE76250.RData")
load("/home/aiusrdata/RCode/TNBC/GSE76250.RData")

summary(GSE76250)  # Gives a summary of the ExpressionSet
str(GSE76250)      # Shows the structure of the ExpressionSet






#                     meta 


#GSE76124  ####################################################

# 198 TNBC samples 


gset <- getGEO("GSE76124", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))




# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

ex <- exprs(gset)

gsms <- paste0(rep("1", ncol(ex)), collapse = "")

sml <- strsplit(gsms, split="")[[1]]


ex_df <- as.data.frame(ex)
feature_data <- fData(gset)
pheno_data <- pData(phenoData(gset))

feature_data$Gene.Symbol <- sub("///.*", "", feature_data$Gene.symbol)
# feature_data$ENTREZ_GENE_ID <- sub(" ///.*", "", feature_data$ENTREZ_GENE_ID)
feature_data$Gene.ID <- sub("///.*", "", feature_data$Gene.ID)

colnames(feature_data)

feature_data <- feature_data %>%
  rename(Gene.ID = ENTREZ_GENE_ID)

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

GSE76124 <- new_gset

max(exprs(GSE76124))
min(exprs(GSE76124))

save(GSE76124, file = "/home/aiusrdata/RCode/TNBC/GSE76124.RData")
load("/home/aiusrdata/RCode/TNBC/GSE76124.RData")

summary(GSE76124)  # Gives a summary of the ExpressionSet
str(GSE76124)      # Shows the structure of the ExpressionSet



#GSE102088 #########################################


gset <- getGEO("GSE102088", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL17586", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))



# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

ex <- exprs(gset)



gsms <- paste0(rep("1", ncol(ex)), collapse = "")

sml <- strsplit(gsms, split="")[[1]]

ex_df <- as.data.frame(ex)
feature_data <- fData(gset)
pheno_data <- pData(phenoData(gset))

# Assuming 'feature_data' is your dataframe and 'gene_assignment' is the column with the data


feature_data$Gene.Symbol <- sapply(strsplit(as.character(feature_data$gene_assignment), " // "), function(x) x[2])
feature_data$Gene.ID <- sapply(strsplit(as.character(feature_data$gene_assignment), " // "), function(x) x[5])

feature_data$Gene.ID1 <- sapply(strsplit(as.character(feature_data$gene_assignment), " // "), function(x) x[1])
feature_data$Gene.Symbol1 <- sapply(strsplit(as.character(feature_data$gene_assignment), " // "), function(x) x[2])

feature_data$Gene.ID <- sapply(strsplit(as.character(feature_data$Gene.ID), "/// "), function(x) x[1])

colnames(feature_data)

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

GSE102088 <- new_gset

max(exprs(GSE102088))
min(exprs(GSE102088))

save(GSE102088, file = "/home/aiusrdata/RCode/TNBC/GSE102088.RData")
load("/home/aiusrdata/RCode/TNBC/GSE102088.RData")

summary(GSE102088)  # Gives a summary of the ExpressionSet
str(GSE102088)      # Shows the structure of the ExpressionSet







#GSE86945 ###########################################


gset <- getGEO("GSE86945", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL17586", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))



# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

ex <- exprs(gset)



gsms <- paste0(rep("1", ncol(ex)), collapse = "")

sml <- strsplit(gsms, split="")[[1]]

ex_df <- as.data.frame(ex)
feature_data <- fData(gset)
pheno_data <- pData(phenoData(gset))

# Assuming 'feature_data' is your dataframe and 'gene_assignment' is the column with the data


feature_data$Gene.Symbol <- sapply(strsplit(as.character(feature_data$gene_assignment), " // "), function(x) x[2])
feature_data$Gene.ID <- sapply(strsplit(as.character(feature_data$gene_assignment), " // "), function(x) x[5])

feature_data$Gene.ID1 <- sapply(strsplit(as.character(feature_data$gene_assignment), " // "), function(x) x[1])
feature_data$Gene.Symbol1 <- sapply(strsplit(as.character(feature_data$gene_assignment), " // "), function(x) x[2])

feature_data$Gene.ID <- sapply(strsplit(as.character(feature_data$Gene.ID), "/// "), function(x) x[1])

colnames(feature_data)

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

GSE86945 <- new_gset

max(exprs(GSE86945))
min(exprs(GSE86945))

save(GSE86945, file = "/home/aiusrdata/RCode/TNBC/GSE86945.RData")
load("/home/aiusrdata/RCode/TNBC/GSE86945.RData")

summary(GSE86945)  # Gives a summary of the ExpressionSet
str(GSE86945)      # Shows the structure of the ExpressionSet



#GSE58812 ############################################

gset <- getGEO("GSE58812", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))



# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

ex <- exprs(gset)

max(ex)
min(ex)

gsms <- paste0(rep("1", ncol(ex)), collapse = "")

sml <- strsplit(gsms, split="")[[1]]


ex_df <- as.data.frame(ex)
feature_data <- fData(gset)
pheno_data <- pData(phenoData(gset))

feature_data$Gene.Symbol <- sub("///.*", "", feature_data$Gene.symbol)
# feature_data$ENTREZ_GENE_ID <- sub(" ///.*", "", feature_data$ENTREZ_GENE_ID)
feature_data$Gene.ID <- sub("///.*", "", feature_data$Gene.ID)

colnames(feature_data)


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

GSE58812 <- new_gset

max(exprs(GSE58812))
min(exprs(GSE58812))

save(GSE58812, file = "/home/aiusrdata/RCode/TNBC/GSE58812.RData")
load("/home/aiusrdata/RCode/TNBC/GSE58812.RData")

summary(GSE58812)  # Gives a summary of the ExpressionSet
str(GSE58812)      # Shows the structure of the ExpressionSet


#GSE83937 ###########################
# treatment
gset <- getGEO("GSE83937", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))




# log2 transformation
ex <- exprs(gset)
max(ex)
min(ex)

qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

ex <- exprs(gset)
max(ex)
min(ex)

gsms <- paste0(rep("1", ncol(ex)), collapse = "")

sml <- strsplit(gsms, split="")[[1]]


ex_df <- as.data.frame(ex)
feature_data <- fData(gset)
pheno_data <- pData(phenoData(gset))

feature_data$Gene.Symbol <- sub("///.*", "", feature_data$Gene.symbol)
# feature_data$ENTREZ_GENE_ID <- sub(" ///.*", "", feature_data$ENTREZ_GENE_ID)
feature_data$Gene.ID <- sub("///.*", "", feature_data$Gene.ID)

colnames(feature_data)


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

GSE83937 <- new_gset

max(exprs(GSE83937))
min(exprs(GSE83937))

save(GSE83937, file = "/home/aiusrdata/RCode/TNBC/GSE83937.RData")
load("/home/aiusrdata/RCode/TNBC/GSE83937.RData")

summary(GSE83937)  # Gives a summary of the ExpressionSet
str(GSE83937)      # Shows the structure of the ExpressionSet




#GSE76275 ########################
# tnbc and not tnbc breast cancer



#meta ###############################################      


# Define a list of dataset filenames
dataset_files <- c("/home/aiusrdata/RCode/TNBC/GSE76250.RData",
                   "/home/aiusrdata/RCode/TNBC/GSE38959.RData",
                   "/home/aiusrdata/RCode/TNBC/GSE45827.RData",
                   "/home/aiusrdata/RCode/TNBC/GSE65194.RData")



gsms_list <- list(
  paste0("11111111111111111111111111111111111111111111111111",
         "11111111111111111111111111111111111111111111111111",
         "11111111111111111111111111111111111111111111111111",
         "111111111111111000000000000000000000000000000000"),
  paste0("11111111111111111111111111111100000000000000000"),
  paste0("1111111111111111111111111111111111111111100000000000"),
  paste0("111111111111111111111111111111111111111111111111111111100000000000")
  
)

gsms <- paste0(gsms_list, collapse = "")
sml <- strsplit(gsms, split="")[[1]]


####################
# sml <- factor(sml, levels = c("0", "1"), labels = c("Normal", "TNBC"))
# ump <- umap(t(allq), n_neighbors = 15, random_state = 123)
# ump_df <- as.data.frame(ump$layout)
# colnames(ump_df) <- c("UMAP1", "UMAP2")
# ump_df$label <- sml
# 
# # Plot using ggplot2
# ggplot(ump_df, aes(x = UMAP1, y = UMAP2, color = label)) +
#   geom_point(size = 1.5) +
#   labs(title = "UMAP plot", x = "UMAP1", y = "UMAP2") +
#   theme_minimal() +
#   theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank()) +
#   scale_color_manual(values = c("blue", "red"))
# 

#####################

# Load datasets
datasets <- lapply(dataset_files, function(file) {
  load(file)
  get(ls()[ls() != "file"])  # Assuming each .RData file loads one object
})

# Extract expression matrices, feature data, and phenotype data
expression_matrices <- lapply(datasets, exprs)
feature_data <- lapply(datasets, fData)
phenotype_data <- lapply(datasets, function(ds) pData(phenoData(ds)))

for (i in seq_along(feature_data)) {
  print(colnames(feature_data[[i]], 1))
  print(dim(expression_matrices[[i]]))
}

# for (i in seq_along(expression_matrices)) {
#   eset <- expression_matrices[[i]]
#   if (min(eset) < 0) {
#     shift_val <- abs(min(eset)) + 1
#     expression_matrices[[i]] <- eset + shift_val
#   }
# }

for (i in seq_along(expression_matrices)) {
  cat(max(expression_matrices[[i]]),min(expression_matrices[[i]]),"\n")
}


# for (i in seq_along(expression_matrices)) {
#   eset <- expression_matrices[[i]]
#   # Normalize between arrays using quantile normalization
#   eset_norm <- normalizeBetweenArrays(eset, method = "quantile")
#   expression_matrices[[i]] <- eset_norm
# }

for (i in seq_along(expression_matrices)) {
  cat(max(expression_matrices[[i]]),min(expression_matrices[[i]]),"\n")
}


# Find common genes across all datasets
common_genes <- Reduce(intersect, lapply(expression_matrices, rownames))
length(common_genes)
# Subset expression matrices to common genes
expression_matrices_common <- lapply(expression_matrices, function(em) em[common_genes, ])

# Combine expression matrices
all_expressions <- do.call(cbind, expression_matrices_common)

# Use feature data from the first dataset
all_fs <- feature_data[[1]][common_genes,]

colnames(phenotype_data[[1]])
phenotype_data[[1]] <- phenotype_data[[1]] %>%
  rename(
    age = `age (yrs):ch1`,  # Correctly reference the age column if that's its name
    disease = `tissue type:ch1`,
    gender = `gender:ch1`
  )

colnames(phenotype_data[[2]])

phenotype_data[[2]] <- phenotype_data[[2]] %>%
  rename(
    age = `age (years):ch1`,  # Ensure this is the correct column name in your dataframe
    disease = `triple-negative status:ch1`,
    gender = `gender:ch1`
  )

phenotype_data[[3]]$gender <- ''
phenotype_data[[3]] <- phenotype_data[[3]] %>%
  rename(
    disease = `diagnosis:ch1`,
    age = `age at diag:ch1`,  # Ensure this is the correct column name in your dataframe
    gender = `gender`
  )

# Ensure all phenotype data have the same columns
all_columns <- Reduce(union, lapply(phenotype_data, names))

phenotype_data_full <- lapply(phenotype_data, function(pd) {
  missing_columns <- setdiff(all_columns, names(pd))
  for (col in missing_columns) {
    pd[[col]] <- NA
  }
  pd <- pd %>% select(all_of(all_columns))
  return(pd)
})

# Combine phenotype data
all_pd <- do.call(rbind, phenotype_data_full)

# Transpose combined expression matrix for UMAP
t_all <- t(all_expressions)

# Prepare batch vector
batch <- factor(unlist(lapply(seq_along(expression_matrices_common), function(i) {
  rep(i - 1, ncol(expression_matrices_common[[i]]))
})))

# Perform UMAP dimensionality reduction
ump <- umap(t_all, n_neighbors = 15, random_state = 123)

# Plotting with adjusted point size
par(mar=c(3,3,2,6), xpd=TRUE)
plot(ump$layout, main="UMAP plot, nbrs=15", xlab="", ylab="", col=batch, pch=20, cex=0.5)  # Adjusted cex for smaller points

# Add legend with adjusted point size
legend("topright", inset=c(-0.15,0), legend=paste0("temp", seq_along(datasets)), pch=20,
       col=1:length(datasets), title="Group", pt.cex=0.5)  # Adjusted pt.cex for smaller legend points

# Adjust batch correction and normalization
allc <- ComBat(all_expressions, batch)
boxplot(allc)
allq <- normalizeQuantiles(allc)
boxplot(allq)

# Perform UMAP dimensionality reduction after normalization
ump <- umap(t(allq), n_neighbors = 5, random_state = 123)

# Plotting
par(mar=c(3,3,2,6), xpd=TRUE)
plot(ump$layout, main="UMAP plot, nbrs=5", xlab="", ylab="", col=batch, pch=20, cex=1.5)

# Add legend
legend("topright", inset=c(-0.15,0), legend=paste0("temp", seq_along(datasets)), pch=20,
       col=1:length(datasets), title="Group", pt.cex=1.5)

# Create ExpressionSet
all_gset <- ExpressionSet(
  assayData = as.matrix(allq),
  phenoData = AnnotatedDataFrame(all_pd),
  featureData = AnnotatedDataFrame(all_fs)  # Assuming all_fs is already prepared
)





save(all_gset, file = "/home/aiusrdata/RCode/TNBC/all_meta.RData")
load("/home/aiusrdata/RCode/TNBC/all_meta.RData")

summary(all_gset)  # Gives a summary of the ExpressionSet
str(all_gset)      # Shows the structure of the ExpressionSet


gset <- all_gset

ex <- exprs(gset) 
dim(ex)

# 
# temp1_gsms <- 0
# 
# 
# temp1_gsms <- "11111111111111111111111111111100000000000000000"
# 
# temp2_gsms <- paste0("11111111111111111111111111111111111111111111111111",
#                      "11111111111111111111111111111111111111111111111111",
#                      "11111111111111111111111111111111111111111111111111",
#                      "111111111111111000000000000000000000000000000000")
# 
# temp3_gsms <- paste0(rep("1", 198), collapse = "")
# 
# 
# temp4_gsms <- paste0(rep("1", 100), collapse = "")
# 
# temp5_gsms <- paste0(rep("1", 107), collapse = "")
# 
# temp6_gsms <- paste0(rep("1", 131), collapse = "") 
# 
# # temp4_gsms <-  paste0(rep("1", 100), collapse = "")
# 
# gsms <- paste0(temp1_gsms, temp2_gsms,temp3_gsms,temp4_gsms,temp5_gsms,temp6_gsms)
# 
# sml <- strsplit(gsms, split="")[[1]]
# 
sml

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

ex <- exprs(gset)
max(ex)
min(ex)

# assign samples to groups and set up design matrix
gs <- factor(sml, levels = c(0, 1), labels = c("normal", "TNBC"))

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

save(tT2, file = "/home/aiusrdata/RCode/TNBC/meta_tT2_test5.RData")



upregulated_all <- tT2[tT2$logFC >= 1 & tT2$adj.P.Val < 0.05, ]
downregulated_all <- tT2[tT2$logFC <= -1 & tT2$adj.P.Val < 0.05, ]

dim(upregulated_all);dim(downregulated_all)

p_value_threshold <- 0.05
tmp1_tT2 <- tT2[abs(tT2$logFC) >= 1 & tT2$adj.P.Val < p_value_threshold, ]
upregulated_tmp1_tT2 <- tT2[tT2$logFC >= 1 & tT2$adj.P.Val < p_value_threshold, ]
downregulated_tmp1_tT2 <- tT2[tT2$logFC <= -1 & tT2$adj.P.Val < p_value_threshold, ]
dim(tmp1_tT2 );dim(upregulated_tmp1_tT2);dim(downregulated_tmp1_tT2 )

kk_temp1 <- enrichKEGG(gene         = as.character(tmp1_tT2 $Gene.ID),
                       organism     = 'hsa',
                       pvalueCutoff = 0.05)



head(kk_temp1@result$Description,10)

hist(tT2$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
     ylab = "Number of genes", main = "P-adj value distribution")

# summarize test results as "up", "down" or "not expressed"
dT <- decideTests(fit2, adjust.method="fdr", p.value=0.05, lfc=0)

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

# MD plot (log fold change vs mean log expression)
# highlight statistically significant (p-adj < 0.05) probes
plotMD(fit2, column=ct, status=dT[,ct], legend=F, pch=20, cex=1)
abline(h=0)

# General expression data analysis
ex <- exprs(gset)

# box-and-whisker plot
dev.new(width=3+ncol(gset)/6, height=5)
ord <- order(gs)  # order samples by group
palette(c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02",
          "#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666"))
par(mar=c(7,4,2,1))
title <- paste ("temp2", "/", annotation(gset), sep ="")
boxplot(ex[,ord], boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=gs[ord])
legend("topleft", groups, fill=palette(), bty="n")
dev.off()

# expression value distribution
par(mar=c(4,4,2,1))
title <- paste ("temp2", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex, group=gs, main=title, legend ="topright")

# UMAP plot (dimensionality reduction)
ex <- na.omit(ex) # eliminate rows with NAs
ex <- ex[!duplicated(ex), ]  # remove duplicates
ump <- umap(t(ex), n_neighbors = 15, random_state = 123)
par(mar=c(3,3,2,6), xpd=TRUE)
plot(ump$layout, main="UMAP plot, nbrs=15", xlab="", ylab="", col=gs, pch=20, cex=1.5)
legend("topright", inset=c(-0.15,0), legend=levels(gs), pch=20,
       col=1:nlevels(gs), title="Group", pt.cex=1.5)
library("maptools")  # point labels without overlaps
pointLabel(ump$layout, labels = rownames(ump$layout), method="SANN", cex=0.6)

# mean-variance trend, helps to see if precision weights are needed
plotSA(fit2, main="Mean variance trend, temp2")




gene_data <- as.data.frame(t(exprs(gset)))
pheno_data <- pData(phenoData(gset))

complete <- merge(gene_data, pheno_data, by.x = "row.names", by.y = "row.names")

data_to_plot <- complete %>%
  select(TOP2A, group)  

# Create the plot
ggplot_object <- ggplot(data_to_plot, aes(x = group, y = TOP2A, color = group)) +
  geom_boxplot() +
  labs(title = "Expression of TOP2A in TNBC and Normal Samples",
       x = "Sample Group",
       y = "Expression Level") +
  theme_minimal()  # Adds a minimalistic theme to the plot

# Print the plot
print(ggplot_object)



emt <- read_excel("/home/aiusrdata/RCode/TNBC/emt.xls",col_names = TRUE)


