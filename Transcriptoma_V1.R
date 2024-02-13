if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.18")
BiocManager::install("TCGAbiolinks")
install.packages("DT")
install.packages("glmnet")
if (!require("caret", quietly = TRUE)) install.packages("caret")
if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx")
if (!requireNamespace("biomaRt", quietly = TRUE))  install.packages("biomaRt")

# Required libraries:
library(TCGAbiolinks)
library(dplyr)
library(tibble)
library(DT)
library(SummarizedExperiment)
library(survival)
library(survminer)
library(readr)
library(glmnet)
library(caret)
library(openxlsx)
library(biomaRt)




#----------------------------recolher os dados

#recolher os dados de miRNA (not working)
# query_miRNA <- GDCquery(
#   project = c("TCGA-GBM", "TCGA-LGG"),
#   data.category = "Transcriptome Profiling",
#   data.type = "Gene Expression Quantification",
#   workflow.type = "STAR - Counts",
#   experimental.strategy = ("miRNA-Seq")
# )
# GDCdownload(query = query_miRNA)
# TCGA_miRNA <- GDCprepare(query = query_miRNA)

#recolher os dados de RNA 
query_RNA <- GDCquery(
  project = c("TCGA-GBM", "TCGA-LGG"),
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  experimental.strategy = ("RNA-Seq")
)
GDCdownload(query = query_RNA)
TCGA_RNA <- GDCprepare(query = query_RNA)

#----------------------------informação sobre os dados

# nrow(TCGA_RNA) #numero de genes
# assayNames(TCGA_RNA) #tipos de dados sobre os genes
# TCGA_RNA$patient #id dos pacientes
# table(TCGA_RNA@rowRanges@elementMetadata@listData[["gene_type"]]) #contagem dos tipos de genes
# TCGA_RNA@colData@listData$ #dados clinicos disponiveis

#----------------------------Curvas de sobrevivência

# Extract the relevant columns
patient <- TCGA_RNA@colData@listData[["patient"]]
project_id <- TCGA_RNA@colData@listData[["project_id"]]
vital_status <- TCGA_RNA@colData@listData[["vital_status"]]
days_to_death <- TCGA_RNA@colData@listData[["days_to_death"]]
days_to_last_follow_up <- TCGA_RNA@colData@listData[["days_to_last_follow_up"]]

# Create a data frame
data_RNA <- data.frame(
  patient = patient,
  project_id = project_id,
  vital_status = vital_status,
  days_to_death = days_to_death,
  days_to_last_follow_up = days_to_last_follow_up
)

# Remove rows where "days to death" is negative and "days to last follow up" is negative
#data_RNA <- data_RNA[data_RNA$days_to_death >= 0 & data_RNA$days_to_last_follow_up >= 0, ]

# apagar as linhas que têm NA no vital status (eram 7)
data_RNA <- data_RNA[complete.cases(data_RNA$vital_status), ]

# apagar as linhas que têm Not reported no vital status (eram 3)
data_RNA <- data_RNA[data_RNA$vital_status != "Not Reported", ]

# alterar o vital_status para 1 em caso de Dead e 0 se Alive
data_RNA <- data_RNA %>%
  mutate(vital_status = ifelse(vital_status == "Alive",0,1))

# criar a coluna days, assumindo os dias até à morte em caso de Dead e os dias até ao último follow up em caso de Alive
data_RNA <- data_RNA %>% mutate(
  days = case_when(
    data_RNA$vital_status == 1 ~ data_RNA$days_to_death,
    data_RNA$vital_status == 0 ~ data_RNA$days_to_last_follow_up
  )
)

# apagar as linhas que têm NA no days (era 1)
data_RNA <- data_RNA[complete.cases(data_RNA$days), ]

# deixar apenas as linhas que têm valores positivos no days
data_RNA <- data_RNA[data_RNA$days > 0, ]

# Keep only distinct rows based on the "patient" column
data_RNA <- distinct(data_RNA, patient, .keep_all = TRUE)

# Kaplan-Meier com a separação por doença
fit_disease_RNA <- survfit(Surv(data_RNA$days, data_RNA$vital_status) ~ data_RNA$project_id, data = data_RNA)
survdiff(Surv(data_RNA$days, data_RNA$vital_status) ~ data_RNA$project_id, data = data_RNA)
ggsurvplot(fit_disease_RNA, data = data_RNA)

#----------------------------pacientes em comum com os da Survival curve feita com os dados clinicos

#dados da curva de sobrevivência
barcode_survival <- read_csv("barcode_survival.csv", col_types = cols(...1 = col_skip(), ))
barcode_survival <- barcode_survival[["x"]]

#dados do transcriptoma
barcode_RNA <- data_RNA$patient

# do código do Gabriel
library(readr)
Data_Gabriel <- read_delim("Glioma_merged_classification_mutations_survival.csv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
barcode_Gabriel<- Data_Gabriel[['Patient_ID']]


#diagrama de Venn
if (!require(devtools)) install.packages("devtools")
devtools::install_github("gaospecial/ggVennDiagram")

library("ggVennDiagram")
ggVennDiagram(x <- list(barcode_survival, barcode_RNA),
              category.names = list("Survival", "RNA"))

library("ggVennDiagram")
ggVennDiagram(x <- list(barcode_Gabriel, barcode_RNA),
              category.names = list("Gabriel", "RNA"))

#---------------------------- Selecionar os genes de interesse

# contagem dos tipos de genes
# table(TCGA_RNA@rowRanges@elementMetadata@listData[["gene_type"]])

# transformar a informação sobre os genes numa dataframe
genes_table <- as.data.frame(TCGA_RNA@rowRanges@elementMetadata@listData)
genes_list <- genes_table$gene_id

# selecionar os tipos de genes de interesse
protein_coding_genes_table <- genes_table[genes_table$gene_type == "protein_coding", ]
long_n_coding_genes_table <- genes_table[genes_table$gene_type == "lncRNA", ]

# selecionar apenas as colunas gene_id e gene_type dos tipos de genes de interesse
protein_coding_genes <- data.frame(gene_id = protein_coding_genes_table$gene_id, gene_type = protein_coding_genes_table$gene_type)
long_n_coding_genes <- data.frame(gene_id = long_n_coding_genes_table$gene_id, gene_type = long_n_coding_genes_table$gene_type)

# dataframe com a expressão de todos os genes
expression_data <- as.data.frame(TCGA_RNA@assays@data@listData$tpm_unstrand)

# adicionar os barcodes ao dataframe com a expressão dos genes de interesse
colnames(expression_data) <- patient

#remover individuos em duplicado
expression_data <- expression_data[, !duplicated(names(expression_data))]

# adicionar o nome dos genes ao dataframe com a expressão de todos os genes
rownames(expression_data) <- genes_list

# dataframe com a expressão dos genes apenas com os individuos dos quais temos informação sobre sobrevivência
expression_data <- expression_data[,colnames(expression_data) %in% data_RNA$patient]

#transpor a dataframe
expression_data <- as.data.frame(t(expression_data))

# dataframe com a expressão dos genes de interesse
protein_coding_expression_data <- expression_data[, colnames(expression_data) %in% protein_coding_genes$gene_id]
long_n_coding_expression_data <- expression_data[, colnames(expression_data) %in% long_n_coding_genes$gene_id]

# remove genes where all expression values are 0 
#no protein coding deixamos de ter 19962 -> 19515, no long non coding deixamos de ter 19601 -> 16278)
protein_coding_expression_data <- protein_coding_expression_data[, colSums(protein_coding_expression_data) != 0]
long_n_coding_expression_data <- long_n_coding_expression_data[, colSums(long_n_coding_expression_data) != 0]

# transform row names into a column
protein_coding_expression_data <- protein_coding_expression_data %>% rownames_to_column(var = "patient")
long_n_coding_expression_data <- long_n_coding_expression_data %>% rownames_to_column(var = "patient")

#---------------------------- Preparar os dados para os testes

# create new dataframe with survival selected data
suvival_data_GBM <- data_RNA[data_RNA$project_id == "TCGA-GBM", ]
suvival_data_GBM <- suvival_data_GBM[, c('patient', 'vital_status', 'days')]

# select only the patients from the gene expression dataframe with survival data 
merged_data <- merge(suvival_data_GBM, protein_coding_expression_data, by = "patient")
protein_coding_expression_data <- merged_data[, c(-2,-3)]

# remove genes where the sum of the gene expression is 0 (necessary to repeat this step, since there is a selection of patients) 
column_sums <- colSums(protein_coding_expression_data[, -1])
columns_to_remove <- names(column_sums[column_sums == 0])
protein_coding_expression_data <- protein_coding_expression_data[, !(names(protein_coding_expression_data) %in% columns_to_remove)]

#create a list of genes to analyse
gene_list_protein_coding <- colnames(protein_coding_expression_data[,-1])

# Create a separate data frame for categorized gene expression
protein_coding_expression_data_categories <- data.frame(patient = protein_coding_expression_data$patient)

# Loop through each gene and categorize based on the median
for (gene in gene_list_protein_coding) {
  # Calculate the median expression for the gene
  median_expression <- median(protein_coding_expression_data[[gene]], na.rm = TRUE)
  
  # Create a grouping variable for the gene based on the median
  protein_coding_expression_data_categories[[gene]] <- ifelse(protein_coding_expression_data[[gene]] > median_expression, "high", "low")
}

# Merge survival data with gene expression data
merged_data_GBM_prot_cod <- merge(suvival_data_GBM, protein_coding_expression_data_categories, by = "patient")

#---------------------------- Log-rank test para protein coding genes na GBM

# Define the survival object
survival_object <- Surv(time = merged_data_GBM_prot_cod$days, event = merged_data_GBM_prot_cod$vital_status)

# Create a dataframe to store log-rank test results for each gene
log_rank_results <- data.frame(Gene = character(0), P_Value = numeric(0))

# Create a list to store genes where an error occurred
genes_with_errors <- list()

# Loop through each gene and perform log-rank test
for (gene in gene_list_protein_coding) {
  result <- tryCatch({
    log_rank_test <- survdiff(survival_object ~ merged_data_GBM_prot_cod[[gene]], data = merged_data_GBM_prot_cod)
    p_value <- log_rank_test$chisq[1]

    # Store results in the data frame
    data.frame(Gene = gene, P_Value = p_value)
  }, error = function(e) {
    cat("Error occurred for gene:", gene, "\n")
    cat("Error message:", conditionMessage(e), "\n")

    # Store the gene where an error occurred in the list
    genes_with_errors <<- c(genes_with_errors, gene)

    # Return a placeholder result to keep the loop going
    data.frame(Gene = character(0), P_Value = numeric(0))
  })

  # If the result is not empty, bind it to the log_rank_results data frame
  if (!is.null(result)) {
    log_rank_results <- rbind(log_rank_results, result)
  }
}

# Adjust p-values for multiple testing (e.g., Bonferroni correction)
log_rank_results$Adjusted_P_Value <- p.adjust(log_rank_results$P_Value, method = "bonferroni")

# Print the results
print(log_rank_results)

#---------------------------- Extract patient ID of GBM > RNA > protein coding with survival data

# # Specify the Excel file path
# excel_file <- "output.xlsx"
# 
# # Create a workbook
# wb <- createWorkbook()
# 
# # Add a sheet to the workbook
# addWorksheet(wb, "Sheet1")
# 
# # Write the dataframe to the sheet
# writeData(wb, sheet = "Sheet1", x = merged_data_GBM_prot_cod$patient, startRow = 1, startCol = 1, rowNames = FALSE)
# 
# # Save the workbook to an Excel file
# saveWorkbook(wb, excel_file)
# 
# # Print a message indicating the file has been created
# cat("Excel file created:", excel_file, "\n")

#---------------------------- Regularized Cox Regression para protein coding genes na GBM

#prepare RNA seq data
RNA_seq_data_cox_regression <- merged_data_GBM_prot_cod
RNA_seq_data_cox_regression <- RNA_seq_data_cox_regression[order(RNA_seq_data_cox_regression$patient), ] #novo
rownames(RNA_seq_data_cox_regression) <- RNA_seq_data_cox_regression$patient
RNA_seq_data_cox_regression <- RNA_seq_data_cox_regression[, -(1:3)]
RNA_seq_data_cox_regression <- data.matrix(RNA_seq_data_cox_regression)

# Define the survival object
merged_data_GBM_prot_cod <- merged_data_GBM_prot_cod[order(merged_data_GBM_prot_cod$patient), ] #novo
survival_object_2 <- Surv(time = merged_data_GBM_prot_cod$days, event = merged_data_GBM_prot_cod$vital_status)

# Use cv.glmnet for cross-validated model selection~
set.seed(1010)
fit <- cv.glmnet(RNA_seq_data_cox_regression, survival_object_2, family = "cox", alpha = 1)  # alpha = 1 for LASSO, 0 for ridge

# Extract coefficients from the fitted model
coefficients <- data.frame(as.matrix(coef(fit, s = "lambda.min")))  # Choose the lambda that minimizes cross-validated error
colnames(coefficients)[colnames(coefficients) == 'X1'] <- "Coef_value"
coefficients_not_zero <- subset(coefficients, Coef_value != 0)
rownames(coefficients_not_zero) <- sub("\\..*", "", rownames(coefficients_not_zero)) # Remove the version of the gene (everything after the dot including the dot)
coefficients_not_zero$ensembl_gene_id <- rownames(coefficients_not_zero)
rownames(coefficients_not_zero) <- NULL

# Print min lambda
cat("lambda min:", "\n", fit$lambda.min)

plot(fit)
plot(fit, xlim = c(log(0.14), log(0.35)), ylim = c(9.1,10))

#---------------------------- Get Information about the genes whose coefficient is not zero


# Function to get gene information
get_gene_info <- function(gene_ids, dataset = "hsapiens_gene_ensembl") {
  
  # Connect to the Ensembl database
  ensembl <- useMart("ensembl", dataset = dataset)
  
  # Get information about the genes
  gene_info <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description"),
                     filters = "ensembl_gene_id",
                     values = gene_ids$ensembl_gene_id,
                     mart = ensembl)
  
  return(gene_info)
}

# Get gene information
gene_info_result <- get_gene_info(gene_ids = coefficients_not_zero)

# Merge gene coefficents and gene information
genes_result <- merge(coefficients_not_zero, gene_info_result, by = "ensembl_gene_id")
genes_result <- genes_result[order(-abs(genes_result$Coef_value)), ]


print(genes_result)

