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

# use for Meta analysis 
#GSE65194 #########################

gset <- getGEO("GSE65194", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

gset <- getGEO("GSE65194", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

fvarLabels(gset) <- make.names(fvarLabels(gset))

gsms <- paste0("00XX0000000X0000000XXXX0000XXXXXX000XXX0000X0000XX",
               "XXXX0XXX00X00XX000000000X0XXXXXXXX000X000000XXXXXX",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXX1X1XX11XX1XX1X11XX1X11XX")
sml <- strsplit(gsms, split="")[[1]]

# filter out excluded samples (marked as "X")
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

ex <- exprs(gset)
dim(ex)
table(sml)

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

dim(GSE65194)
table(sml)


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


gsms <- paste0("000XXX00XX0X00X0X0X010100X0XXXXXXXXXXXX1X1XXXXXXXX",
               "XXXXXXXXXXXXXXXXX1X1XXXXXXXXXX11X111XXXXXXXXXXXXXX",
               "XXXXXXXXXX1XX1XXXXXXXXXXXXXXXXXXXX1XXXXXXXX11XXXXX",
               "XXXXXXXXXXX1XXXXXX1X1XXXXXXX1XXXXXX1XXX1XXXXXXXXXX",
               "XXXXXXXXXXXXXX1XX1XXXXXX11XXXXXX1XXX1X1XXXX1XXXXXX",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX1XXXXXXXX",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXXXXX")
sml <- strsplit(gsms, split="")[[1]]

# filter out excluded samples (marked as "X")
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]


ex <- exprs(gset)

table(sml)
dim(ex)

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

gsms <- paste0("00X0XX00000XXXX0X00XXX0X00000000XX000XXX00XX0X00XX",
               "X00X00XX000000000XXXXXXXXXXXXXXXXXX11111111111XXXX",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXXX")
sml <- strsplit(gsms, split="")[[1]]

# filter out excluded samples (marked as "X")
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

ex <- exprs(gset)

table(sml)
dim(ex)



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
load("/home/aiusrdata/RCode/TNBC/GSE45827.RData")

summary(GSE45827)  # Gives a summary of the ExpressionSet
str(GSE45827)      # Shows the structure of the ExpressionSet




#GSE61724 ####


gset <- getGEO("GSE61724", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL6244", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

gsms <- "XXXXXXXXX1X0XXXXXXXXXX000111XXX00XX0XX0XXXXXXX00X0XXX0X0XXXX0XX00XXX"
sml <- strsplit(gsms, split="")[[1]]

# filter out excluded samples (marked as "X")
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]



ex <- exprs(gset)
table(sml)
dim(ex)



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
feature_data$Gene.ID <- sub("///.*", "", feature_data$Gene.ID)


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


ex <- exprs(gset)
table(sml)
dim(ex)




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
feature_data$Gene.ID <- sub("///.*", "", feature_data$Gene.ID)



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
feature_data$Gene.ID <- sub("///.*", "", feature_data$Gene.ID)



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
dim(exprs(GSE37751))

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
gsms <- "0000000000000000000000000000001111111111111XXXX"
sml <- strsplit(gsms, split="")[[1]]

# filter out excluded samples (marked as "X")
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

ex <- exprs(gset)
dim(ex)
table(sml)


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

table(sml)

dim(ex)

#GSE76250 ###################################################################


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

GSE76250 <- new_gset

max(exprs(GSE76250))
min(exprs(GSE76250))

save(GSE76250, file = "/home/aiusrdata/RCode/TNBC/GSE76250.RData")
load("/home/aiusrdata/RCode/TNBC/GSE76250.RData")

summary(GSE76250)  # Gives a summary of the ExpressionSet
str(GSE76250)      # Shows the structure of the ExpressionSet

dim(ex)
table(sml)




#                     meta 




# not used in the meta





#not used in meta 

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





# meta common platforms GPL570 #####################



# Define a list of dataset filenames
dataset_files1 <- c("/home/aiusrdata/RCode/TNBC/GSE65194.RData",
                    "/home/aiusrdata/RCode/TNBC/GSE31448.RData",
                    "/home/aiusrdata/RCode/TNBC/GSE45827.RData")



gsms_list1 <- list(
  paste0("000000000000000000000000000000000000000000000000000000011111111111"),
  paste0("0000000000010100011111111111111111111111111111"),
  paste0("0000000000000000000000000000000000000000011111111111")  
)

for (gsms in gsms_list1) {
  # Split the string into individual characters
  chars <- strsplit(gsms, split = "")[[1]]
  
  # Count the occurrences of '0' and '1'
  count_0 <- sum(chars == '0')
  count_1 <- sum(chars == '1')
  
  # Print the counts
  cat("TNBCs:", count_0, "Normals:", count_1, "\n")
}


gsms <- paste0(gsms_list1, collapse = "")
sml <- strsplit(gsms, split="")[[1]]



# Load datasets
datasets <- lapply(dataset_files1, function(file) {
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


common_genes <- Reduce(intersect, lapply(expression_matrices, rownames))
length(common_genes)

expression_matrices_common <- lapply(expression_matrices, function(em) em[common_genes, ])

all_expressions <- do.call(cbind, expression_matrices_common)

all_fs <- feature_data[[1]][common_genes,]

colnames(phenotype_data[[1]])
phenotype_data[[1]] <- phenotype_data[[1]] %>%
  rename(
    age = `age (yrs):ch1`,  
    disease = `tissue type:ch1`,
    gender = `gender:ch1`
  )

colnames(phenotype_data[[2]])

phenotype_data[[2]] <- phenotype_data[[2]] %>%
  rename(
    age = `age (years):ch1`,  
    disease = `triple-negative status:ch1`,
    gender = `gender:ch1`
  )

phenotype_data[[3]]$gender <- ''
phenotype_data[[3]] <- phenotype_data[[3]] %>%
  rename(
    disease = `diagnosis:ch1`,
    age = `age at diag:ch1`, 
    gender = `gender`
  )


all_columns <- Reduce(union, lapply(phenotype_data, names))

phenotype_data_full <- lapply(phenotype_data, function(pd) {
  missing_columns <- setdiff(all_columns, names(pd))
  for (col in missing_columns) {
    pd[[col]] <- NA
  }
  pd <- pd %>% select(all_of(all_columns))
  return(pd)
})

all_pd <- do.call(rbind, phenotype_data_full)

t_all <- t(all_expressions)

batch_labels <- c("GSE65194", "GSE31448", "GSE45827")
batch <- factor(unlist(lapply(seq_along(expression_matrices_common), function(i) {
  rep(i - 1, ncol(expression_matrices_common[[i]]))
})), labels = batch_labels)


allc <- ComBat(all_expressions, batch)
allq <- normalizeQuantiles(allc)

# boxplot(allc)
# boxplot(allq)



#####################################
# bpx plot part 

library(sva)
res <- 700

batch_labels <- c("GSE65194", "GSE31448", "GSE45827")

group_info <- data.frame(
  Dataset = rep(batch_labels, times = c(66, 46, 52)),
  Group = factor(c(rep("TNBC", 55), rep("Normal", 11),
                   rep("TNBC", 15), rep("Normal", 31),
                   rep("TNBC", 41), rep("Normal", 11)))
)

ord <- order(group_info$Dataset) 

colors <- c("GSE65194" = "#1B9E77", "GSE31448" = "#7570B3", "GSE45827" = "#E7298A")

png(filename="/home/aiusrdata/RCode/TNBC/plots/gpl570_boxplot_before_normalization.png", width=12, height=6, units="in", res=res)
par(mar=c(5, 6, 2, 4))  # Similar adjustment for the second plot
par(cex.axis=0.7)
title1 <- "GPL570 - Before Quantile Normalization"
boxplot(allc[,ord], boxwex=0.6, notch=TRUE, main=title1, outline=FALSE, las=2, col=colors[group_info$Dataset[ord]])
legend("topright", inset=c(1.02,-0.05), legend=batch_labels, fill=colors, bty="c", cex=0.7, pt.cex=0.5, xpd=TRUE)
dev.off()

png(filename="/home/aiusrdata/RCode/TNBC/plots/gpl570_boxplot_after_normalization.png", width=12, height=6, units="in", res=res)
par(mar=c(5, 6, 2, 4))  # Similar adjustment for the second plot
par(cex.axis=0.7)
title2 <- "GPL570 - After Quantile Normalization"
boxplot(allq[,ord], boxwex=0.6, notch=TRUE, main=title2, outline=FALSE, las=2, col=colors[group_info$Dataset[ord]])
legend("topright", inset=c(1.02,-0.05), legend=batch_labels, fill=colors, bty="c", cex=0.7, pt.cex=0.5, xpd=TRUE)
dev.off()


#####################################
# #pca and umap and identify the outliers

# Load necessary libraries
library(ggplot2)
library(umap)

# Assuming 'allq', 'allc', 'all_expressions', 'batch', and 'group_info' are already loaded in your environment
datasets <- list(allq = allq, allc = allc, all_expressions = all_expressions)
titles <- list(
  allq = "After Combat and Quantile Normalization",
  allc = "After Combat",
  all_expressions = "Before Batch Effect Removal"
)
base_path <- "/home/aiusrdata/RCode/TNBC/plots"

if (!dir.exists(base_path)) {
  dir.create(base_path, recursive = TRUE)
}

perform_analysis <- function(data, dataset_name) {
  # Perform PCA
  pca_res <- prcomp(t(data), scale. = TRUE)
  pca_df <- as.data.frame(pca_res$x[, 1:2])
  colnames(pca_df) <- c("PC1", "PC2")
  pca_df$batch <- batch
  
  # Perform UMAP
  umap_res <- umap(t(data), n_neighbors = 15, random_state = 123)
  umap_df <- as.data.frame(umap_res$layout)
  colnames(umap_df) <- c("UMAP1", "UMAP2")
  umap_df$batch <- batch
  
  # Identify coordinates of GSM781355
  gsm_pca <- pca_df[row.names(pca_df) == "GSM781355", ]
  gsm_umap <- umap_df[row.names(umap_df) == "GSM781355", ]
  
  # Offset for label
  offset <- 1
  
  # Theme for white background and publication quality
  publication_theme <- theme_minimal() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.border = element_rect(color = "black", fill = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(hjust = 0.5, color = "black")
    )
  
  # Plot PCA
  pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, shape = batch, color = group_info$Group)) +
    geom_point(size = 2) +
    geom_text(aes(x = gsm_pca$PC1 + offset, y = gsm_pca$PC2 + offset, label = "GSM781355"), color = "red", hjust = 0) +
    geom_segment(aes(x = gsm_pca$PC1 + offset, y = gsm_pca$PC2 + offset, xend = gsm_pca$PC1, yend = gsm_pca$PC2), 
                 arrow = arrow(length = unit(0.3, "cm")), color = "red") +
    publication_theme +
    labs(title = paste("PCA -", titles[[dataset_name]]), x = "PC1", y = "PC2") +
    scale_color_discrete(name = "Batch")
  ggsave(sprintf("%s/GPL_570_%s_PCA.png", base_path, dataset_name), plot = pca_plot, width = 10, height = 8, dpi = 300)
  
  # Plot UMAP
  umap_plot <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, shape = batch, color = group_info$Group)) +
    geom_point(size = 2) +
    geom_text(aes(x = gsm_umap$UMAP1 + offset, y = gsm_umap$UMAP2 + offset, label = "GSM781355"), color = "red", hjust = 0) +
    geom_segment(aes(x = gsm_umap$UMAP1 + offset, y = gsm_umap$UMAP2 + offset, xend = gsm_umap$UMAP1, yend = gsm_umap$UMAP2), 
                 arrow = arrow(length = unit(0.3, "cm")), color = "red") +
    publication_theme +
    labs(title = paste("UMAP -", titles[[dataset_name]]), x = "UMAP1", y = "UMAP2") +
    scale_color_discrete(name = "Batch")
  ggsave(sprintf("%s/GPL_570_%s_UMAP.png", base_path, dataset_name), plot = umap_plot, width = 10, height = 8, dpi = 300)
}

for (dataset_name in names(datasets)) {
  perform_analysis(datasets[[dataset_name]], dataset_name)
}


# 
# datasets <- list(allq = allq, allc = allc, all_expressions = all_expressions)
# titles <- list(
#   allq = "After Combat and Quantile Normalization",
#   allc = "After Combat",
#   all_expressions = "Before Batch Effect Removal"
# )
# base_path <- "/home/aiusrdata/RCode/TNBC/plots"
# 
# if (!dir.exists(base_path)) {
#   dir.create(base_path, recursive = TRUE)
# }
# 
# perform_analysis <- function(data, dataset_name) {
#   # Perform PCA
#   pca_res <- prcomp(t(data), scale. = TRUE)
#   pca_df <- as.data.frame(pca_res$x[, 1:2])
#   colnames(pca_df) <- c("PC1", "PC2")
#   pca_df$batch <- batch
#   
#   # Perform UMAP
#   umap_res <- umap(t(data), n_neighbors = 15, random_state = 123)
#   umap_df <- as.data.frame(umap_res$layout)
#   colnames(umap_df) <- c("UMAP1", "UMAP2")
#   umap_df$batch <- batch
#   
#   # Theme for white background and publication quality
#   publication_theme <- theme_minimal() +
#     theme(
#       plot.background = element_rect(fill = "white", color = NA),
#       panel.background = element_rect(fill = "white", color = NA),
#       panel.border = element_rect(color = "black", fill = NA),
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       legend.background = element_rect(fill = "white", color = NA),
#       plot.title = element_text(hjust = 0.5, color = "black")
#     )
#   
#   # Plot PCA
#   pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, shape = batch, color = group_info$Group)) +
#     geom_point(size = 2) +
#     publication_theme +
#     labs(title = paste("PCA -", titles[[dataset_name]]), x = "PC1", y = "PC2") +
#     scale_color_discrete(name = "Batch")
#   ggsave(sprintf("%s/GPL_570_%s_PCA.png", base_path, dataset_name), plot = pca_plot, width = 10, height = 8, dpi = 300)
#   
#   # Plot UMAP
#   umap_plot <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, shape = batch, color = group_info$Group)) +
#     geom_point(size = 2) +
#     publication_theme +
#     labs(title = paste("UMAP -", titles[[dataset_name]]), x = "UMAP1", y = "UMAP2") +
#     scale_color_discrete(name = "Batch")
#   ggsave(sprintf("%s/GPL_570_%s_UMAP.png", base_path, dataset_name), plot = umap_plot, width = 10, height = 8, dpi = 300)
# }
# 
# for (dataset_name in names(datasets)) {
#   perform_analysis(datasets[[dataset_name]], dataset_name)
# }


#####################################




# Create ExpressionSet
all_gset <- ExpressionSet(
  assayData = as.matrix(allq),
  phenoData = AnnotatedDataFrame(all_pd),
  featureData = AnnotatedDataFrame(all_fs)  # Assuming all_fs is already prepared
)





save(all_gset, file = "/home/aiusrdata/RCode/TNBC/meta_gpl570.RData")
load("/home/aiusrdata/RCode/TNBC/meta_gpl570.RData")

dim(all_gset)
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

###############################################
#Find and Remove Outlier 

datExpr <- t(ex)
dim(datExpr)

sampleTree = hclust(dist(datExpr), method = "average");

# Load required packages
library(WGCNA)

# Define a function to adjust cex based on the number of samples
adjust_cex <- function(n_samples) {
  if (n_samples > 500) {
    return(0.3)
  } else if (n_samples > 100) {
    return(0.5)
  } else {
    return(0.7)
  }
}

# Assuming datExpr, sampleTree, and gset are already defined

# Calculate appropriate cex for the given number of samples
cex_value <- adjust_cex(nrow(datExpr))

# Plot the initial dendrogram and save as PDF
pdf("/home/aiusrdata/git_projects/bioinformatics/GPL570_initial_sample_outliers_tree.pdf", width = 8, height = 6)
par(cex = cex_value)  # Adjust text size based on sample size
par(mar = c(5, 4, 4, 2) + 0.1)  # Adjust margins
plot(sampleTree, main = "Sample Clustering to Detect Outliers",
     sub = "", xlab = "", ylab = "Height",
     cex.lab = cex_value, cex.axis = cex_value, cex.main = cex_value,
     cex = cex_value)  # Adjust text size
abline(h = 150, col = "red", lwd = 2)  # Add a horizontal line
dev.off()

# Cut the tree to identify clusters with a cut height of 120
clust = cutreeStatic(sampleTree, cutHeight = 150, minSize = 10)
table(clust)

# Identify outliers (samples assigned to cluster 0)
outliers = which(clust == 0)
outlier_samples = rownames(datExpr)[outliers]
outlier_samples
# Replace outliers in sml with "X"
sml[outliers] <- "X"

# Print the modified sml to verify
print(sml)

# Remove outliers from the dataset
datExpr_clean = datExpr[-outliers, ]

# Re-cluster the samples without outliers
sampleTree_clean = hclust(dist(datExpr_clean), method = "average")

# Calculate appropriate cex for the cleaned dataset
cex_value_clean <- adjust_cex(nrow(datExpr_clean))

# Plot the cleaned dendrogram and save as PDF
pdf("/home/aiusrdata/git_projects/bioinformatics/GPL570_cleaned_sample_outliers_tree.pdf", width = 8, height = 6)
par(cex = cex_value_clean)  # Adjust text size based on sample size
par(mar = c(5, 4, 4, 2) + 0.1)  # Adjust margins
plot(sampleTree_clean, main = "Sample Clustering After Removing Outliers",
     sub = "", xlab = "", ylab = "Height",
     cex.lab = cex_value_clean, cex.axis = cex_value_clean, cex.main = cex_value_clean,
     cex = cex_value_clean)  # Adjust text size
abline(h = 150, col = "red", lwd = 2)  # Add a horizontal line
dev.off()

# Cut the tree to identify clusters with a cut height of 120 for the cleaned data
clust_clean = cutreeStatic(sampleTree_clean, cutHeight = 150, minSize = 10)
table(clust_clean)

# Identify outliers (samples assigned to cluster 0) in the cleaned data
outliers_clean = which(clust_clean == 0)
outlier_samples_clean = rownames(datExpr_clean)[outliers_clean]
outlier_samples_clean

# Remove columns corresponding to outlier samples from the `ex` dataframe
ex_clean <- ex[, !colnames(ex) %in% outlier_samples]

# Print the first few rows of the cleaned dataframe to verify
head(ex_clean)

# Remove "X" values from sml
sml <- sml[sml != "X"]

# Print the modified sml to verify
print(sml)

# Remove columns corresponding to outlier samples from the `gset` dataset
gset_clean <- gset[, !colnames(gset) %in% outlier_samples]

# Print the first few rows of the cleaned gset to verify
head(gset_clean)

ex <- ex_clean
dim(ex)
gset<-gset_clean
###############################################


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

save(tT2, file = "/home/aiusrdata/RCode/TNBC/meta_tT2_gpl570.RData")

load("/home/aiusrdata/RCode/TNBC/meta_tT2_gpl570.RData")


upregulated_all <- tT2[tT2$logFC >= 1.5 & tT2$adj.P.Val <= 0.01, ]
downregulated_all <- tT2[tT2$logFC <= -1.5 & tT2$adj.P.Val <= 0.01, ]

dim(upregulated_all);dim(downregulated_all)

p_value_threshold <- 0.01
tmp1_tT2 <- tT2[abs(tT2$logFC) >= 1.5 & tT2$adj.P.Val < p_value_threshold, ]
upregulated_tmp1_tT2 <- tT2[tT2$logFC >= 1.5 & tT2$adj.P.Val < p_value_threshold, ]
downregulated_tmp1_tT2 <- tT2[tT2$logFC <= -1.5 & tT2$adj.P.Val < p_value_threshold, ]
dim(tmp1_tT2 );dim(upregulated_tmp1_tT2);dim(downregulated_tmp1_tT2 )

kk_temp1 <- enrichKEGG(gene         = as.character(tmp1_tT2 $Gene.ID),
                       organism     = 'hsa',
                       pvalueCutoff = 0.05)



head(kk_temp1@result$Description,10)



#############################################

hist_data <- hist(tT2$adj.P.Val, plot = FALSE)

colors <- ifelse(hist_data$breaks[-length(hist_data$breaks)] <= 0.01, "green", "grey")

plot(hist_data, col = colors, border = "white", xlab = "P-adj",
     ylab = "Number of genes", main = "P-adjusted Value Distribution")

legend("topright", legend = c("p-adjust <= 0.01", "p-adjust > 0.01"), 
       fill = c("green", "grey"), title = "Legend", cex = 0.8, bty = "n")


range_logFC <- range(tT2$logFC, na.rm = TRUE)  # Remove NA values if any
bin_width <- 0.1  # Define the bin width

breaks <- seq(from = floor(range_logFC[1] / bin_width) * bin_width,
              to = ceiling(range_logFC[2] / bin_width) * bin_width,
              by = bin_width)

max_count <- max(hist(tT2$logFC, plot = FALSE, breaks = breaks)$counts)

tT2$color <- ifelse(tT2$logFC >= 1.5, "red", ifelse(tT2$logFC <= -1.5, "blue", "grey"))
ggplot(tT2, aes(x = logFC, fill = color)) +
  geom_histogram(binwidth = 0.1, color = "white") +  # Adjust binwidth as needed
  scale_fill_identity() +
  geom_vline(xintercept = c(1.5, -1.5), color = c("red", "blue"), linetype = "dashed", linewidth = 1) +
  scale_y_continuous(limits = c(0, max_count * 1.2)) +  # Adjust y limits to provide space for text
  geom_text(aes(label = paste("n =", sum(logFC >= 1.5)), y = max_count * 1.1, x = 1.75), color = "red", vjust = 0, hjust = 0) +
  geom_text(aes(label = paste("n =", sum(logFC <= -1.5)), y = max_count * 1.1, x = -1.75), color = "blue", vjust = 0, hjust = 1) +
  labs(x = "Log Fold Change", y = "Number of Genes", title = "Log Fold Change Distribution") +
  theme_minimal()

#############################################


dT <- decideTests(fit2, adjust.method="fdr", p.value=0.05, lfc=1.5)

# Venn diagram of results
vennDiagram(dT, circle.col=palette())

# create Q-Q plot for t-statistic
t.good <- which(!is.na(fit2$F)) # filter out bad probes
qqt(fit2$t[t.good], fit2$df.total[t.good], main="Moderated t statistic")

# volcano plot (log P-value vs log fold change)
colnames(fit2) # list contrast names
ct <- 1        # choose contrast of interest


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








# meta common platforms GPL6244 #####################



# Define a list of dataset filenames
dataset_files2 <- c("/home/aiusrdata/RCode/TNBC/GSE61724.RData",
                    "/home/aiusrdata/RCode/TNBC/GSE36295.RData",
                    "/home/aiusrdata/RCode/TNBC/GSE37751.RData")



gsms_list2 <- list(
  paste0("10000111000000000000"),
  paste0("1111100000000000"),
  paste0("1111111111111111111111111111111111111111111111100000000000000")  
)

for (gsms in gsms_list2) {
  # Split the string into individual characters
  chars <- strsplit(gsms, split = "")[[1]]
  
  # Count the occurrences of '0' and '1'
  count_0 <- sum(chars == '0')
  count_1 <- sum(chars == '1')
  
  # Print the counts
  cat("TNBCs:", count_0, "Normals:", count_1, "\n")
}


gsms <- paste0(gsms_list2, collapse = "")
sml <- strsplit(gsms, split="")[[1]]



# Load datasets
datasets <- lapply(dataset_files2, function(file) {
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


allc <- ComBat(all_expressions, batch)
allq <- normalizeQuantiles(allc)
# boxplot(allq)
# boxplot(allc)

#####################################
# box plot part 

library(sva)

res <- 700

batch_labels <- c("GSE61724", "GSE36295", "GSE37751")

group_info <- data.frame(
  Dataset = rep(batch_labels, times = c(20, 16, 61)),
  Group = factor(c(rep("TNBC", 16), rep("Normal", 4),
                   rep("TNBC", 11), rep("Normal", 5),
                   rep("TNBC", 14), rep("Normal", 47)))
)


ord <- order(group_info$Dataset) 

colors <- c("GSE61724" = "#1B9E77", "GSE36295" = "#7570B3", "GSE37751" = "#E7298A")

png(filename="/home/aiusrdata/RCode/TNBC/plots/gpl6244_boxplot_before_combat.png", width=12, height=6, units="in", res=res)
par(mar=c(5, 6, 2, 4))  # Similar adjustment for the second plot
par(cex.axis=0.7)
title1 <- "GPL6244 - Before Combat"
boxplot(all_expressions[,ord], boxwex=0.6, notch=TRUE, main=title1, outline=FALSE, las=2, col=colors[group_info$Dataset[ord]])
legend("topright", inset=c(1.02,-0.05), legend=batch_labels, fill=colors, bty="c", cex=0.7, pt.cex=0.5, xpd=TRUE)
dev.off()



png(filename="/home/aiusrdata/RCode/TNBC/plots/gpl6244_boxplot_before_normalization.png", width=12, height=6, units="in", res=res)
par(mar=c(5, 6, 2, 4))  # Similar adjustment for the second plot
par(cex.axis=0.7)
title1 <- "GPL6244 - Before Quantile Normalization"
boxplot(allc[,ord], boxwex=0.6, notch=TRUE, main=title1, outline=FALSE, las=2, col=colors[group_info$Dataset[ord]])
legend("topright", inset=c(1.02,-0.05), legend=batch_labels, fill=colors, bty="c", cex=0.7, pt.cex=0.5, xpd=TRUE)
dev.off()

png(filename="/home/aiusrdata/RCode/TNBC/plots/gpl6244_boxplot_after_normalization.png", width=12, height=6, units="in", res=res)
par(mar=c(5, 6, 2, 4))  # Similar adjustment for the second plot
par(cex.axis=0.7)
title2 <- "GPL6244 - After Quantile Normalization"
boxplot(allq[,ord], boxwex=0.6, notch=TRUE, main=title2, outline=FALSE, las=2, col=colors[group_info$Dataset[ord]])
legend("topright", inset=c(1.02,-0.05), legend=batch_labels, fill=colors, bty="c", cex=0.7, pt.cex=0.5, xpd=TRUE)
dev.off()



#####################################
#pca and umap

library(ggplot2)
library(umap)
library(dplyr)

# Assuming 'batch' and 'group_info' are pre-defined and correctly aligned with your samples
# Make sure that 'batch' is a factor and 'group_info' is a dataframe with a 'Group' column

datasets <- list(allq = allq, allc = allc)
titles <- list(
  allq = "After Combat and Quantile Normalization",
  allc = "After Combat",
  all_expressions = "Before Batch Effect Removal"
)
base_path <- "/home/aiusrdata/RCode/TNBC/plots"

if (!dir.exists(base_path)) {
  dir.create(base_path, recursive = TRUE)
}

perform_analysis <- function(data, dataset_name) {
  # Perform PCA
  pca_res <- prcomp(t(data), scale. = TRUE)
  pca_df <- as.data.frame(pca_res$x[, 1:2])
  colnames(pca_df) <- c("PC1", "PC2")
  pca_df$batch <- batch
  pca_df$sample <- colnames(data)
  
  # Perform UMAP
  umap_res <- umap(t(data), n_neighbors = 15, random_state = 123)
  umap_df <- as.data.frame(umap_res$layout)
  colnames(umap_df) <- c("UMAP1", "UMAP2")
  umap_df$batch <- batch
  umap_df$sample <- colnames(data)
  
  # Theme for white background and publication quality
  publication_theme <- theme_minimal() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.border = element_rect(color = "black", fill = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(hjust = 0.5, color = "black")
    )
  
  # Highlight specific samples
  highlight_samples <- c("GSM1511827", "GSM927002", "GSM927015", "GSM927019", "GSM927048", "GSM927059")
  
  # Plot PCA
  pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, shape = batch, color = group_info$Group, label = sample)) +
    geom_point(size = 2) +
    geom_text(data = filter(pca_df, sample %in% highlight_samples), vjust = -1, color = "black", size = 3) +
    publication_theme +
    labs(title = paste("PCA -", titles[[dataset_name]]), x = "PC1", y = "PC2") +
    scale_color_discrete(name = "Group")
  ggsave(sprintf("%s/GPL_6244_%s_PCA.png", base_path, dataset_name), plot = pca_plot, width = 10, height = 8, dpi = 300)
  
  # Plot UMAP
  umap_plot <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, shape = batch, color = group_info$Group, label = sample)) +
    geom_point(size = 2) +
    geom_text(data = filter(umap_df, sample %in% highlight_samples), vjust = -1, color = "black", size = 3) +
    publication_theme +
    labs(title = paste("UMAP -", titles[[dataset_name]]), x = "UMAP1", y = "UMAP2") +
    scale_color_discrete(name = "Group")
  ggsave(sprintf("%s/GPL_6244_%s_UMAP.png", base_path, dataset_name), plot = umap_plot, width = 10, height = 8, dpi = 300)
}

for (dataset_name in names(datasets)) {
  perform_analysis(datasets[[dataset_name]], dataset_name)
}


datasets <- list(all_expressions = all_expressions)
titles <- list(
  allq = "After Combat and Quantile Normalization",
  allc = "After Combat",
  all_expressions = "Before Batch Effect Removal"
)
base_path <- "/home/aiusrdata/RCode/TNBC/plots"

if (!dir.exists(base_path)) {
  dir.create(base_path, recursive = TRUE)
}

perform_analysis <- function(data, dataset_name) {
  # Perform PCA
  pca_res <- prcomp(t(data), scale. = TRUE)
  pca_df <- as.data.frame(pca_res$x[, 1:2])
  colnames(pca_df) <- c("PC1", "PC2")
  pca_df$batch <- batch
  
  # Perform UMAP
  umap_res <- umap(t(data), n_neighbors = 15, random_state = 123)
  umap_df <- as.data.frame(umap_res$layout)
  colnames(umap_df) <- c("UMAP1", "UMAP2")
  umap_df$batch <- batch
  
  # Theme for white background and publication quality
  publication_theme <- theme_minimal() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.border = element_rect(color = "black", fill = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(hjust = 0.5, color = "black")
    )
  
  # Plot PCA
  pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, shape = batch, color = group_info$Group)) +
    geom_point(size = 2) +
    publication_theme +
    labs(title = paste("PCA -", titles[[dataset_name]]), x = "PC1", y = "PC2") +
    scale_color_discrete(name = "Batch")
  ggsave(sprintf("%s/GPL_6244_%s_PCA.png", base_path, dataset_name), plot = pca_plot, width = 10, height = 8, dpi = 300)
  
  # Plot UMAP
  umap_plot <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, shape = batch, color = group_info$Group)) +
    geom_point(size = 2) +
    publication_theme +
    labs(title = paste("UMAP -", titles[[dataset_name]]), x = "UMAP1", y = "UMAP2") +
    scale_color_discrete(name = "Batch")
  ggsave(sprintf("%s/GPL_6244_%s_UMAP.png", base_path, dataset_name), plot = umap_plot, width = 10, height = 8, dpi = 300)
}

for (dataset_name in names(datasets)) {
  perform_analysis(datasets[[dataset_name]], dataset_name)
}
#####################################




# Create ExpressionSet
all_gset <- ExpressionSet(
  assayData = as.matrix(allq),
  phenoData = AnnotatedDataFrame(all_pd),
  featureData = AnnotatedDataFrame(all_fs)  # Assuming all_fs is already prepared
)





save(all_gset, file = "/home/aiusrdata/RCode/TNBC/meta_gpl6244.RData")
load("/home/aiusrdata/RCode/TNBC/meta_gpl6244.RData")

dim(all_gset)
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

###############################################
#Find and Remove Outlier 

datExpr <- t(ex)
dim(datExpr)

sampleTree = hclust(dist(datExpr), method = "average");

# Load required packages
library(WGCNA)

# Define a function to adjust cex based on the number of samples
adjust_cex <- function(n_samples) {
  if (n_samples > 500) {
    return(0.3)
  } else if (n_samples > 100) {
    return(0.5)
  } else {
    return(0.7)
  }
}

# Assuming datExpr, sampleTree, and gset are already defined

# Calculate appropriate cex for the given number of samples
cex_value <- adjust_cex(nrow(datExpr))

# Plot the initial dendrogram and save as PDF
pdf("/home/aiusrdata/git_projects/bioinformatics/GPL6244_initial_sample_outliers_tree.pdf", width = 8, height = 6)
par(cex = cex_value)  # Adjust text size based on sample size
par(mar = c(5, 4, 4, 2) + 0.1)  # Adjust margins
plot(sampleTree, main = "Sample Clustering to Detect Outliers",
     sub = "", xlab = "", ylab = "Height",
     cex.lab = cex_value, cex.axis = cex_value, cex.main = cex_value,
     cex = cex_value)  # Adjust text size
abline(h = 100, col = "red", lwd = 2)  # Add a horizontal line
dev.off()

# Cut the tree to identify clusters with a cut height of 120
clust = cutreeStatic(sampleTree, cutHeight = 100, minSize = 10)
table(clust)

# Identify outliers (samples assigned to cluster 0)
outliers = which(clust == 0)
outlier_samples = rownames(datExpr)[outliers]
outlier_samples

# Replace outliers in sml with "X"
sml[outliers] <- "X"

# Print the modified sml to verify
print(sml)

# Remove outliers from the dataset
datExpr_clean = datExpr[-outliers, ]

# Re-cluster the samples without outliers
sampleTree_clean = hclust(dist(datExpr_clean), method = "average")

# Calculate appropriate cex for the cleaned dataset
cex_value_clean <- adjust_cex(nrow(datExpr_clean))

# Plot the cleaned dendrogram and save as PDF
pdf("/home/aiusrdata/git_projects/bioinformatics/GPL6244_cleaned_sample_outliers_tree.pdf", width = 8, height = 6)
par(cex = cex_value_clean)  # Adjust text size based on sample size
par(mar = c(5, 4, 4, 2) + 0.1)  # Adjust margins
plot(sampleTree_clean, main = "Sample Clustering After Removing Outliers",
     sub = "", xlab = "", ylab = "Height",
     cex.lab = cex_value_clean, cex.axis = cex_value_clean, cex.main = cex_value_clean,
     cex = cex_value_clean)  # Adjust text size
abline(h = 100, col = "red", lwd = 2)  # Add a horizontal line
dev.off()

# Cut the tree to identify clusters with a cut height of 120 for the cleaned data
clust_clean = cutreeStatic(sampleTree_clean, cutHeight = 100, minSize = 10)
table(clust_clean)

# Identify outliers (samples assigned to cluster 0) in the cleaned data
outliers_clean = which(clust_clean == 0)
outlier_samples_clean = rownames(datExpr_clean)[outliers_clean]
outlier_samples_clean

# Remove columns corresponding to outlier samples from the `ex` dataframe
ex_clean <- ex[, !colnames(ex) %in% outlier_samples]

# Print the first few rows of the cleaned dataframe to verify
head(ex_clean)

# Remove "X" values from sml
sml <- sml[sml != "X"]

# Print the modified sml to verify
print(sml)

# Remove columns corresponding to outlier samples from the `gset` dataset
gset_clean <- gset[, !colnames(gset) %in% outlier_samples]

# Print the first few rows of the cleaned gset to verify
head(gset_clean)

ex <- ex_clean
dim(ex)
gset<-gset_clean
###############################################


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

save(tT2, file = "/home/aiusrdata/RCode/TNBC/meta_tT2_meta_gpl6244.RData")



upregulated_all <- tT2[tT2$logFC >= 1.5 & tT2$adj.P.Val <+ 0.01, ]
downregulated_all <- tT2[tT2$logFC <= -1.5 & tT2$adj.P.Val <+ 0.01, ]

dim(upregulated_all);dim(downregulated_all)

p_value_threshold <- 0.01
tmp1_tT2 <- tT2[abs(tT2$logFC) >= 1.5 & tT2$adj.P.Val < p_value_threshold, ]
upregulated_tmp1_tT2 <- tT2[tT2$logFC >= 1.5 & tT2$adj.P.Val < p_value_threshold, ]
downregulated_tmp1_tT2 <- tT2[tT2$logFC <= -1.5 & tT2$adj.P.Val < p_value_threshold, ]
dim(tmp1_tT2 );dim(upregulated_tmp1_tT2);dim(downregulated_tmp1_tT2 )

kk_temp1 <- enrichKEGG(gene         = as.character(tmp1_tT2 $Gene.ID),
                       organism     = 'hsa',
                       pvalueCutoff = 0.05)

head(kk_temp1@result$Description,10)


gene_vector1 <- setNames(tmp1_tT2$logFC, tmp1_tT2$ID)

sorted_gene_vector1 <- sort(gene_vector1, decreasing = TRUE)
kk3 <- gseKEGG(geneList     = sorted_gene_vector1,
               organism     = 'hsa',
               minGSSize    = 120,
               pvalueCutoff = 0.05,
               verbose      = FALSE)


###

gseKEGG(geneList     = sorted_gene_vector1,
        organism     = 'hsa',
        minGSSize    = 120,
        pvalueCutoff = 0.05,
        verbose      = FALSE)

#############################################

hist_data <- hist(tT2$adj.P.Val, plot = FALSE)

colors <- ifelse(hist_data$breaks[-length(hist_data$breaks)] <= 0.01, "green", "grey")

plot(hist_data, col = colors, border = "white", xlab = "P-adj",
     ylab = "Number of genes", main = "P-adjusted Value Distribution")

legend("topright", legend = c("p-adjust <= 0.01", "p-adjust > 0.01"), 
       fill = c("green", "grey"), title = "Legend", cex = 0.8, bty = "n")


range_logFC <- range(tT2$logFC, na.rm = TRUE)  # Remove NA values if any
bin_width <- 0.1  # Define the bin width

breaks <- seq(from = floor(range_logFC[1] / bin_width) * bin_width,
              to = ceiling(range_logFC[2] / bin_width) * bin_width,
              by = bin_width)

max_count <- max(hist(tT2$logFC, plot = FALSE, breaks = breaks)$counts)

tT2$color <- ifelse(tT2$logFC >= 1.5, "red", ifelse(tT2$logFC <= -1.5, "blue", "grey"))
ggplot(tT2, aes(x = logFC, fill = color)) +
  geom_histogram(binwidth = 0.1, color = "white") +  # Adjust binwidth as needed
  scale_fill_identity() +
  geom_vline(xintercept = c(1.5, -1.5), color = c("red", "blue"), linetype = "dashed", linewidth = 1) +
  scale_y_continuous(limits = c(0, max_count * 1.2)) +  # Adjust y limits to provide space for text
  geom_text(aes(label = paste("n =", sum(logFC >= 1.5)), y = max_count * 1.1, x = 1.75), color = "red", vjust = 0, hjust = 0) +
  geom_text(aes(label = paste("n =", sum(logFC <= -1.5)), y = max_count * 1.1, x = -1.75), color = "blue", vjust = 0, hjust = 1) +
  labs(x = "Log Fold Change", y = "Number of Genes", title = "Log Fold Change Distribution") +
  theme_minimal()

#############################################







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








#meta whole ###############################################      


# Define a list of dataset filenames
dataset_files <- c("/home/aiusrdata/RCode/TNBC/GSE65194.RData",
                   "/home/aiusrdata/RCode/TNBC/GSE31448.RData",
                   "/home/aiusrdata/RCode/TNBC/GSE45827.RData",
                   "/home/aiusrdata/RCode/TNBC/GSE61724.RData",
                   "/home/aiusrdata/RCode/TNBC/GSE36295.RData",
                   "/home/aiusrdata/RCode/TNBC/GSE37751.RData",
                   "/home/aiusrdata/RCode/TNBC/GSE38959.RData",
                   "/home/aiusrdata/RCode/TNBC/GSE76250.RData")



gsms_list <- list(
  paste0("000000000000000000000000000000000000000000000000000000011111111111"),
  paste0("0000000000010100011111111111111111111111111111"),
  paste0("0000000000000000000000000000000000000000011111111111"),
  paste0("10000111000000000000"),
  paste0("1111100000000000"),
  paste0("1111111111111111111111111111111111111111111111100000000000000"),
  paste0("0000000000000000000000000000001111111111111"),
  paste0("000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000111111111111111111111111111111111")
  
)

for (gsms in gsms_list) {
  # Split the string into individual characters
  chars <- strsplit(gsms, split = "")[[1]]
  
  # Count the occurrences of '0' and '1'
  count_0 <- sum(chars == '0')
  count_1 <- sum(chars == '1')
  
  # Print the counts
  cat("TNBCs:", count_0, "Normals:", count_1, "\n")
}


gsms <- paste0(gsms_list, collapse = "")
sml <- strsplit(gsms, split="")[[1]]


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


# Adjust batch correction and normalization
allc <- ComBat(all_expressions, batch)
allq <- normalizeQuantiles(allc)
# boxplot(allc)
# boxplot(allq)


#####################################
# box plot part 

library(sva)

res <- 700

batch_labels <- c("GSE65194", "GSE31448", "GSE45827","GSE61724","GSE36295","GSE37751","GSE38959","GSE76250")

group_info <- data.frame(
  Dataset = rep(batch_labels, times = c(66,46,52,20,16,61,43,198)),
  Group = factor(c(rep("TNBC", 55), rep("Normal", 11),
                   rep("TNBC", 15), rep("Normal", 31),
                   rep("TNBC", 41), rep("Normal", 11),
                   rep("TNBC", 16), rep("Normal", 4),
                   rep("TNBC", 11), rep("Normal", 5),
                   rep("TNBC", 14), rep("Normal", 47),
                   rep("TNBC", 30), rep("Normal", 13),
                   rep("TNBC", 165), rep("Normal", 33)))
)


ord <- order(group_info$Dataset) 

colors <- c(
  "GSE65194" = "#1B9E77",  # Existing color
  "GSE31448" = "#D95F02",  # New color
  "GSE45827" = "#7570B3",  # New color
  "GSE61724" = "#E7298A",  # New color
  "GSE36295" = "#66A61E",  # Existing color
  "GSE37751" = "#E6AB02",  # Existing color
  "GSE38959" = "#A6761D",  # New color
  "GSE76250" = "#666666"   # New color
)

png(filename="/home/aiusrdata/RCode/TNBC/plots/all_datasets_boxplot_before_normalization.png", width=12, height=6, units="in", res=res)
par(mar=c(5, 6, 2, 4))  # Similar adjustment for the second plot
par(cex.axis=0.7)
title1 <- "All Datasets - Before Quantile Normalization"
boxplot(allc[,ord], boxwex=0.6, notch=TRUE, main=title1, outline=FALSE, las=2, col=colors[group_info$Dataset[ord]])
legend("topright", inset=c(1.02,-0.05), legend=batch_labels, fill=colors, bty="c", cex=0.7, pt.cex=0.5, xpd=TRUE)
dev.off()

png(filename="/home/aiusrdata/RCode/TNBC/plots/all_datasets_boxplot_after_normalization.png", width=12, height=6, units="in", res=res)
par(mar=c(5, 6, 2, 4))  # Similar adjustment for the second plot
par(cex.axis=0.7)
title2 <- "All Datasets - After Quantile Normalization"
boxplot(allq[,ord], boxwex=0.6, notch=TRUE, main=title2, outline=FALSE, las=2, col=colors[group_info$Dataset[ord]])
legend("topright", inset=c(1.02,-0.05), legend=batch_labels, fill=colors, bty="c", cex=0.7, pt.cex=0.5, xpd=TRUE)
dev.off()



#####################################
#pca and umap

datasets <- list(allq = allq, allc = allc, all_expressions = all_expressions)
titles <- list(
  allq = "After Combat and Quantile Normalization",
  allc = "After Combat",
  all_expressions = "Before Batch Effect Removal"
)
base_path <- "/home/aiusrdata/RCode/TNBC/plots"

if (!dir.exists(base_path)) {
  dir.create(base_path, recursive = TRUE)
}

perform_analysis <- function(data, dataset_name) {
  # Perform PCA
  pca_res <- prcomp(t(data), scale. = TRUE)
  pca_df <- as.data.frame(pca_res$x[, 1:2])
  colnames(pca_df) <- c("PC1", "PC2")
  pca_df$batch <- batch
  
  # Perform UMAP
  umap_res <- umap(t(data), n_neighbors = 15, random_state = 123)
  umap_df <- as.data.frame(umap_res$layout)
  colnames(umap_df) <- c("UMAP1", "UMAP2")
  umap_df$batch <- batch
  
  # Theme for white background and publication quality
  publication_theme <- theme_minimal() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.border = element_rect(color = "black", fill = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(hjust = 0.5, color = "black")
    )
  
  # Plot PCA
  pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, shape = batch, color = group_info$Group)) +
    geom_point(size = 2) +
    publication_theme +
    labs(title = paste("PCA -", titles[[dataset_name]]), x = "PC1", y = "PC2") +
    scale_color_discrete(name = "Batch")
  ggsave(sprintf("%s/All_Datasets_%s_PCA.png", base_path, dataset_name), plot = pca_plot, width = 10, height = 8, dpi = 300)
  
  # Plot UMAP
  umap_plot <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, shape = batch, color = group_info$Group)) +
    geom_point(size = 2) +
    publication_theme +
    labs(title = paste("UMAP -", titles[[dataset_name]]), x = "UMAP1", y = "UMAP2") +
    scale_color_discrete(name = "Batch")
  ggsave(sprintf("%s/All_Datasets_%s_UMAP.png", base_path, dataset_name), plot = umap_plot, width = 10, height = 8, dpi = 300)
}

for (dataset_name in names(datasets)) {
  perform_analysis(datasets[[dataset_name]], dataset_name)
}


#####################################



# Create ExpressionSet
all_gset <- ExpressionSet(
  assayData = as.matrix(allq),
  phenoData = AnnotatedDataFrame(all_pd),
  featureData = AnnotatedDataFrame(all_fs)  # Assuming all_fs is already prepared
)





save(all_gset, file = "/home/aiusrdata/RCode/TNBC/all_meta.RData")
load("/home/aiusrdata/RCode/TNBC/all_meta.RData")

dim(all_gset)
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

###############################################
#Find and Remove Outlier 

datExpr <- t(ex)
dim(datExpr)

sampleTree = hclust(dist(datExpr), method = "average");

# Load required packages
library(WGCNA)

# Define a function to adjust cex based on the number of samples
adjust_cex <- function(n_samples) {
  if (n_samples > 500) {
    return(0.3)
  } else if (n_samples > 100) {
    return(0.5)
  } else {
    return(0.7)
  }
}

# Assuming datExpr, sampleTree, and gset are already defined

# Calculate appropriate cex for the given number of samples
cex_value <- adjust_cex(nrow(datExpr))

# Plot the initial dendrogram and save as PDF
pdf("/home/aiusrdata/git_projects/bioinformatics/initial_sample_outliers_tree.pdf", width = 8, height = 6)
par(cex = cex_value)  # Adjust text size based on sample size
par(mar = c(5, 4, 4, 2) + 0.1)  # Adjust margins
plot(sampleTree, main = "Sample Clustering to Detect Outliers",
     sub = "", xlab = "", ylab = "Height",
     cex.lab = cex_value, cex.axis = cex_value, cex.main = cex_value,
     cex = cex_value)  # Adjust text size
abline(h = 90, col = "red", lwd = 2)  # Add a horizontal line
dev.off()

# Cut the tree to identify clusters with a cut height of 120
clust = cutreeStatic(sampleTree, cutHeight = 120, minSize = 10)
table(clust)

# Identify outliers (samples assigned to cluster 0)
outliers = which(clust == 0)
outlier_samples = rownames(datExpr)[outliers]

# Replace outliers in sml with "X"
sml[outliers] <- "X"

# Print the modified sml to verify
print(sml)

# Remove outliers from the dataset
datExpr_clean = datExpr[-outliers, ]

# Re-cluster the samples without outliers
sampleTree_clean = hclust(dist(datExpr_clean), method = "average")

# Calculate appropriate cex for the cleaned dataset
cex_value_clean <- adjust_cex(nrow(datExpr_clean))

# Plot the cleaned dendrogram and save as PDF
pdf("/home/aiusrdata/git_projects/bioinformatics/cleaned_sample_outliers_tree.pdf", width = 8, height = 6)
par(cex = cex_value_clean)  # Adjust text size based on sample size
par(mar = c(5, 4, 4, 2) + 0.1)  # Adjust margins
plot(sampleTree_clean, main = "Sample Clustering After Removing Outliers",
     sub = "", xlab = "", ylab = "Height",
     cex.lab = cex_value_clean, cex.axis = cex_value_clean, cex.main = cex_value_clean,
     cex = cex_value_clean)  # Adjust text size
abline(h = 90, col = "red", lwd = 2)  # Add a horizontal line
dev.off()

# Cut the tree to identify clusters with a cut height of 120 for the cleaned data
clust_clean = cutreeStatic(sampleTree_clean, cutHeight = 120, minSize = 80)
table(clust_clean)

# Identify outliers (samples assigned to cluster 0) in the cleaned data
outliers_clean = which(clust_clean == 0)
outlier_samples_clean = rownames(datExpr_clean)[outliers_clean]
outlier_samples_clean

# Remove columns corresponding to outlier samples from the `ex` dataframe
ex_clean <- ex[, !colnames(ex) %in% outlier_samples]

# Print the first few rows of the cleaned dataframe to verify
head(ex_clean)

# Remove "X" values from sml
sml <- sml[sml != "X"]

# Print the modified sml to verify
print(sml)

# Remove columns corresponding to outlier samples from the `gset` dataset
gset_clean <- gset[, !colnames(gset) %in% outlier_samples]

# Print the first few rows of the cleaned gset to verify
head(gset_clean)

ex <- ex_clean
dim(ex)
gset<-gset_clean


save(gset, file = "/home/aiusrdata/RCode/TNBC/ro_all_meta.RData")
load("/home/aiusrdata/RCode/TNBC/ro_all_meta.RData")


###############################################

# assign samples to groups and set up design matrix

gs <- factor(sml)
groups <- make.names(c("TNBC","normal"))
levels(gs) <- groups
gset$group <- gs

################################################################################

#export dataset as csv file for ML 

# Transpose the ex dataframe
csv_file <- t(ex)

# Ensure the length of sml matches the number of rows in csv_file
if (length(sml) == nrow(csv_file)) {
  # Add sml as a new column to csv_file
  csv_file <- as.data.frame(csv_file)
  csv_file$Label <- sml
  
  # Save the csv_file dataframe as a CSV file
  write.csv(csv_file, "/home/aiusrdata/git_projects/bioinformatics/cleaned_data_with_labels.csv", row.names = TRUE)
} else {
  stop("The length of sml does not match the number of rows in csv_file.")
}

# Print the first few rows of the modified csv_file to verify
head(csv_file)
################################################################################

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



upregulated_all <- tT2[tT2$logFC >= 1.5 & tT2$adj.P.Val < 0.05, ]
downregulated_all <- tT2[tT2$logFC <= -1.5 & tT2$adj.P.Val < 0.05, ]

dim(upregulated_all);dim(downregulated_all)

p_value_threshold <- 0.01
tmp1_tT2 <- tT2[abs(tT2$logFC) >= 1.5 & tT2$adj.P.Val < p_value_threshold, ]
upregulated_tmp1_tT2 <- tT2[tT2$logFC >= 1.5 & tT2$adj.P.Val < p_value_threshold, ]
downregulated_tmp1_tT2 <- tT2[tT2$logFC <= -1.5 & tT2$adj.P.Val < p_value_threshold, ]
dim(tmp1_tT2 );dim(upregulated_tmp1_tT2);dim(downregulated_tmp1_tT2 )

kk_temp1 <- enrichKEGG(gene         = as.character(tmp1_tT2 $Gene.ID),
                       organism     = 'hsa',
                       pvalueCutoff = 0.05)



head(kk_temp1@result$Description,10)

#############################################

hist_data <- hist(tT2$adj.P.Val, plot = FALSE)

colors <- ifelse(hist_data$breaks[-length(hist_data$breaks)] <= 0.01, "green", "grey")

plot(hist_data, col = colors, border = "white", xlab = "P-adj",
     ylab = "Number of genes", main = "P-adjusted Value Distribution")

legend("topright", legend = c("p-adjust <= 0.01", "p-adjust > 0.01"), 
       fill = c("green", "grey"), title = "Legend", cex = 0.8, bty = "n")


range_logFC <- range(tT2$logFC, na.rm = TRUE)  # Remove NA values if any
bin_width <- 0.1  # Define the bin width

breaks <- seq(from = floor(range_logFC[1] / bin_width) * bin_width,
              to = ceiling(range_logFC[2] / bin_width) * bin_width,
              by = bin_width)

max_count <- max(hist(tT2$logFC, plot = FALSE, breaks = breaks)$counts)

tT2$color <- ifelse(tT2$logFC >= 1.5, "red", ifelse(tT2$logFC <= -1.5, "blue", "grey"))
ggplot(tT2, aes(x = logFC, fill = color)) +
  geom_histogram(binwidth = 0.1, color = "white") +  # Adjust binwidth as needed
  scale_fill_identity() +
  geom_vline(xintercept = c(1.5, -1.5), color = c("red", "blue"), linetype = "dashed", linewidth = 1) +
  scale_y_continuous(limits = c(0, max_count * 1.2)) +  # Adjust y limits to provide space for text
  geom_text(aes(label = paste("n =", sum(logFC >= 1.5)), y = max_count * 1.1, x = 1.75), color = "red", vjust = 0, hjust = 0) +
  geom_text(aes(label = paste("n =", sum(logFC <= -1.5)), y = max_count * 1.1, x = -1.75), color = "blue", vjust = 0, hjust = 1) +
  labs(x = "Log Fold Change", y = "Number of Genes", title = "Log Fold Change Distribution") +
  theme_minimal()

#############################################



# summarize test results as "up", "down" or "not expressed"
dT <- decideTests(fit2, adjust.method="fdr", p.value=0.01, lfc=1.5)

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
















#analysis the result #####################

env1 <- new.env()
env2 <- new.env()
env3 <- new.env()
env4 <- new.env()

load("/home/aiusrdata/RCode/TNBC/results/GSE76250_tT2.RData", envir = env1)
load("/home/aiusrdata/RCode/TNBC/results/GSE38959_tT2.RData", envir = env2)
load("/home/aiusrdata/RCode/TNBC/meta_tT2_gpl570.RData", envir = env3)
load("/home/aiusrdata/RCode/TNBC/meta_tT2_meta_gpl6244.RData", envir = env4)
                                 
process_dataset <- function(env, dataset_name) {
  tT2 <- env[[dataset_name]]
  p_value_threshold <- 0.01
  
  tmp1_tT2 <- tT2[abs(tT2$logFC) >= 1.5 & tT2$adj.P.Val < p_value_threshold, ]
  upregulated_tmp1_tT2 <- tT2[tT2$logFC >= 1.5 & tT2$adj.P.Val < p_value_threshold, ]
  downregulated_tmp1_tT2 <- tT2[tT2$logFC <= -1.5 & tT2$adj.P.Val < p_value_threshold, ]
  cat(dim(tmp1_tT2))
  return(list(
    tmp1 = tmp1_tT2,
    upregulated = upregulated_tmp1_tT2,
    downregulated = downregulated_tmp1_tT2
  ))
}

results_env1 <- process_dataset(env1, "tT2")
results_env2 <- process_dataset(env2, "tT2")
results_env3 <- process_dataset(env3, "tT2")
results_env4 <- process_dataset(env4, "tT2")

dim(results_env1$tmp1)
dim(results_env2$tmp1)
dim(results_env3$tmp1)
dim(results_env4$tmp1)


intersect_genes <- inner_join(results_env1$tmp1, results_env2$tmp1, by = "Gene.ID")
length(intersect_genes)


intersect_genes_all <- inner_join(results_env1$tmp1, results_env2$tmp1, by = "Gene.ID") %>%
  inner_join(results_env3$tmp1, by = "Gene.ID") %>%
  inner_join(results_env4$tmp1, by = "Gene.ID")

length(intersect_genes_all)

kk_temp1 <- enrichKEGG(gene         = as.character(intersect_genes_all$Gene.ID),
                       organism     = 'hsa',
                       pvalueCutoff = 0.05)

head(kk_temp1@result$Description,10)


# Assuming the gene IDs are in a column named 'Gene.ID'
gene_ids_env1 <- results_env1$tmp1$Gene.ID
gene_ids_env2 <- results_env2$tmp1$Gene.ID
gene_ids_env3 <- results_env3$tmp1$Gene.ID
gene_ids_env4 <- results_env4$tmp1$Gene.ID

# Corrected code to draw a Venn diagram for four groups
venn.plot <- draw.quad.venn(
  area1 = length(gene_ids_env1),
  area2 = length(gene_ids_env2),
  area3 = length(gene_ids_env3),
  area4 = length(gene_ids_env4),
  n12 = length(intersect(gene_ids_env1, gene_ids_env2)),
  n13 = length(intersect(gene_ids_env1, gene_ids_env3)),
  n14 = length(intersect(gene_ids_env1, gene_ids_env4)),
  n23 = length(intersect(gene_ids_env2, gene_ids_env3)),
  n24 = length(intersect(gene_ids_env2, gene_ids_env4)),
  n34 = length(intersect(gene_ids_env3, gene_ids_env4)),
  n123 = length(Reduce(intersect, list(gene_ids_env1, gene_ids_env2, gene_ids_env3))),
  n124 = length(Reduce(intersect, list(gene_ids_env1, gene_ids_env2, gene_ids_env4))),
  n134 = length(Reduce(intersect, list(gene_ids_env1, gene_ids_env3, gene_ids_env4))),
  n234 = length(Reduce(intersect, list(gene_ids_env2, gene_ids_env3, gene_ids_env4))),
  n1234 = length(intersect(intersect(intersect(gene_ids_env1, gene_ids_env2), gene_ids_env3), gene_ids_env4)),
  category = c("GSE76250", "GSE38959", "GPL570", "GPL622"),
  fill = c("skyblue", "pink1", "mediumorchid", "orange"),
  label.col = "black",
  cex = 2,
  cat.cex = 2,
  cat.col = c("darkblue", "darkred", "darkgreen", "darkorange")
)

# Draw the Venn diagram
grid.draw(venn.plot)



# End Code ###############
