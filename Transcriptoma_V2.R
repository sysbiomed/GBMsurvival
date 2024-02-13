if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.18")
BiocManager::install("TCGAbiolinks")
BiocManager::install("DESeq2")                      #novo
BiocManager::install("edgeR")                       #novo
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
library(writexl)
library(edgeR)

#-----------------------------------------------FUNÇÕES-------------------------

#---------------------------- Tratar os dados de sobrevivência

# Extract survival data of interest
createSurvivalDataFrame <- function(column_names, source_data) {
  
  data_frame <- data.frame(matrix(ncol = length(column_names), nrow = nrow(source_data@colData)))
  colnames(data_frame) <- column_names

  for (col_name in column_names) {
    # Extract the column from the source and add it to the dataframe
    data_frame[[col_name]] <- source_data@colData@listData[[col_name]]
  }

  return(data_frame)
}

# Clean and prepare survival data
cleanSurvivalData <- function(data) {
  
  # apagar as linhas que têm NA no vital status
  data <- data[complete.cases(data$vital_status), ]
  
  # apagar as linhas que têm Not reported no vital status
  data <- data[data$vital_status != "Not Reported", ]
  
  # alterar o vital_status para 1 em caso de Dead e 0 se Alive
  data <- data %>%
    mutate(vital_status = ifelse(vital_status == "Alive",0,1))
  
  # criar a coluna days, assumindo os dias até à morte em caso de Dead e os dias até ao último follow up em caso de Alive
  data <- data %>% mutate(
    days = case_when(
      data$vital_status == 1 ~ data$days_to_death,
      data$vital_status == 0 ~ data$days_to_last_follow_up
    )
  )
  
  # apagar as linhas que têm NA no days
  data <- data[complete.cases(data$days), ]
  
  # deixar apenas as linhas que têm valores positivos no days
  data <- data[data$days > 0, ]
   
  # Keep only distinct rows based on the "patient" column
  data <- distinct(data, patient, .keep_all = TRUE)
  
  return(data)
}

#---------------------------- Tratar os dados de expressão génica

# Extrair dataframe com a expressão dos genes de interesse
createGeneExpressionDataframe <- function(genes_of_interest, expression_measure, source_data, survivalData) {
  
  
  genes_table <- as.data.frame(source_data@rowRanges@elementMetadata@listData) # Tranformar a informação sobre os genes numa dataframe
  genes_table <- genes_table[genes_table$gene_type == genes_of_interest, ] #selecionar os genes de interesse
  
  expression_data <- as.data.frame(source_data@assays@data@listData[[expression_measure]]) #dataframe com a expressão de todos os genes
  colnames(expression_data) <- source_data@colData@listData[["patient"]] # adicionar os barcodes ao dataframe
  expression_data <- expression_data[, !duplicated(names(expression_data))] #remover individuos em duplicado
  rownames(expression_data) <- source_data@rowRanges@elementMetadata@listData$gene_id # adicionar o nome dos genes ao dataframe
  expression_data <- expression_data[,colnames(expression_data) %in% survivalData$patient] # selecionar apenas os individuos dos quais temos informação sobre sobrevivência
  expression_data <- as.data.frame(t(expression_data)) #transpor a dataframe
  
  data_frame <- expression_data[, colnames(expression_data) %in% genes_table$gene_id] #seleciinar a expressão dos genes de interesse
  
  return(data_frame)
}

# Clean the gene expression dataframe
cleanGeneExpressionData <- function(expression_data) {
  
  # Check if "patient" column exists in the dataframe
  if ("patient" %in% colnames(expression_data)) {
    # If "patient" column exists, use function #1
    rownames(expression_data) <- expression_data$patient # transform patients column in rowsnames
    expression_data <- expression_data[, -1] # delete the column patients 
    expression_data <- expression_data[, colSums(expression_data) != 0] #remove genes where all expression values are 0
    expression_data <- expression_data[rowSums(expression_data) != 0, ] #remove samples where the sum of the gene expression is 0
    expression_data <- expression_data %>% rownames_to_column(var = "patient") # transform row names into a column again
  } else {
    # If "patient" column doesn't exist, use function #2
    expression_data <- expression_data[, colSums(expression_data) != 0] #remove genes where all expression values are 0
    expression_data <- expression_data[rowSums(expression_data) != 0, ] #remove samples where the sum of the gene expression is 0
    expression_data <- expression_data %>% rownames_to_column(var = "patient") # transform row names into a column
    
  }
  
  return(expression_data)
}

#---------------------------- Preparar os dados para os testes

# selecionar apenas os dados da doença de interesse
selectDataPerDisease <- function(disease, expression_data, survivalData) {
  
  #survival data
  survivalData <- survivalData[survivalData$project_id == disease, ] # select the survival data from the patients with the disease
  survivalData <- survivalData[, c('patient', 'vital_status', 'days')] # keep only the columns of interest
  #expression data
  merged_data <- merge(survivalData, expression_data, by = "patient") # select only the patients from the gene expression dataframe with the desired disease and select only the patients whose the sum of the gene expression is not zero
  expression_data <- merged_data[, c(-2,-3)] # delete the columns added above
  
  result <- list(survival_data = survivalData, genes_expression = expression_data)
  
  return(result)
}


# ordenar os dados por paciente
orderDataByPatient <- function(data) {
  
  data <- data[order(data$patient), ]
  
  return(data)
}

# create a dataframe with categorized gene expression (high and low)
categorizeGeneExpressionData <- function(expression_data) {
  
  categorized_expression_data <- data.frame(patient = expression_data$patient) # Create a data frame for categorized gene expression
  genes_list <- colnames(expression_data[,-1]) #create a list of genes
  
  # Loop through each gene and categorize based on the median
  for (gene in genes_list) {
    # Calculate the median expression for the gene
    median_expression <- median(expression_data[[gene]])
    
    # Create a grouping variable for the gene based on the median
    categorized_expression_data[[gene]] <- ifelse(expression_data[[gene]] > median_expression, "high", "low")
  }
  
  return(categorized_expression_data)
}

#---------------------------- Regularized Cox Regression 

#preparar os dados para a Regularized Cox Regression
prepareDataForCoxRegression <- function(data) { 
  
  rownames(data) <- data$patient
  data <- data[, -1]
  data <- data.matrix(data)
  
  return(data)
}

# Extract coefficients from the fitted model
extractCoxRegressionCoefficients <- function(cox_fit) { 

  coefficients <- data.frame(as.matrix(coef(cox_fit, s = "lambda.min")))  # Choose the lambda that minimizes cross-validated error
  colnames(coefficients)[colnames(coefficients) == 'X1'] <- "Coef_value" # Change the column name 
  coefficients <- subset(coefficients, Coef_value != 0) # select only the coefficients that are not equal to zero
  rownames(coefficients) <- sub("\\..*", "", rownames(coefficients)) # Remove the version of the gene from its name(everything after the dot including the dot)
  coefficients$ensembl_gene_id <- rownames(coefficients)
  rownames(coefficients) <- NULL
  
  # Print min lambda
  cat("Lambda min:", "\n", cox_fit$lambda.min)

  return(coefficients)
}

#---------------------------- Get Information about the genes whose coefficient is not zero

# Function to get gene information
getGenesInfo <- function(coefficients) {
  
  # Connect to the Ensembl database
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  # Get information about the genes
  gene_info <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description"),
                     filters = "ensembl_gene_id",
                     values = coefficients$ensembl_gene_id,
                     mart = ensembl)
  
  # Merge gene coefficents and gene information
  gene_info <- merge(coefficients, gene_info, by = "ensembl_gene_id")
  gene_info <- gene_info[order(-abs(gene_info$Coef_value)), ]
  
  print(gene_info)
  
  return(gene_info)
}

#---------------------------- Extract patient ID 

# Function to extract a column from a data frame and write it to an Excel file
extractSampleIDColumnToExcel <- function(data, column_name, output_file) {
  column <- as.data.frame(data[[column_name]])
  
  # Write to Excel
  write_xlsx(column, output_file)
}

#------------------------------------APLICAR-AS-FUNÇÕES-------------------------

#---------------------------- recolher todos dados

query_RNA <- GDCquery(
  project = c("TCGA-GBM", "TCGA-LGG"),
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  experimental.strategy = ("RNA-Seq")
)
GDCdownload(query = query_RNA)
TCGA_RNA <- GDCprepare(query = query_RNA)

#---------------------------- Tratar os dados de sobrevivência

# Extract survival data of interest
columns_to_extract <- c("patient", "project_id", "vital_status", "days_to_death", "days_to_last_follow_up")
survival_data <- createSurvivalDataFrame(columns_to_extract, TCGA_RNA)


# Clean and prepare survival data
survival_data <- cleanSurvivalData(survival_data)
print(survival_data)

#---------------------------- Plot da curva de sobrevivência

# Kaplan-Meier com a separação por doença
fit_disease_RNA <- survfit(Surv(survival_data$days, survival_data$vital_status) ~ survival_data$project_id)
survdiff(Surv(survival_data$days, survival_data$vital_status) ~ survival_data$project_id)
ggsurvplot(fit_disease_RNA, data = survival_data)

#---------------------------- Tratar os dados de expressão génica

# Extrair dataframe com a expressão dos genes de interesse
genes_expression <- createGeneExpressionDataframe("protein_coding", "tpm_unstrand", TCGA_RNA, survival_data)
#genes_expression <- createGeneExpressionDataframe("lncRNA", "tpm_unstrand", TCGA_RNA, survival_data)       #<------------ indicar os dados

# Clean the gene expression dataframe
genes_expression <- cleanGeneExpressionData(genes_expression)

#---------------------------- Preparar os dados para os testes

# selecionar apenas os dados da doença de interesse
data <- selectDataPerDisease("TCGA-GBM", genes_expression, survival_data)
#data <- selectDataPerDisease("TCGA-LGG", genes_expression, survival_data)                                   #<------------ indicar a doença
genes_expression <- data$genes_expression
survival_data <- data$survival_data

# remove genes where the sum of the gene expression is 0 (necessary to repeat this step, since there is a selection of patients) 
genes_expression <- cleanGeneExpressionData(genes_expression)

# ordenar os dados por paciente
genes_expression <- orderDataByPatient(genes_expression)
survival_data <- orderDataByPatient(survival_data)

# create a dataframe with categorized gene expression (high and low)
categorized_genes_expression <- categorizeGeneExpressionData(genes_expression)

#---------------------------- Regularized Cox Regression 

# preparar os dados para a Regularized Cox Regression
#dataCox <- prepareDataForCoxRegression(genes_expression)
dataCox <- prepareDataForCoxRegression(categorized_genes_expression)                                        #<------------ indicar se queremos os dados categorizados

# Define survival object
survival_object <- Surv(time = survival_data$days, event = survival_data$vital_status)

# set seed
set.seed(1010)

# fit the Regularized Cox Regression
fit <- cv.glmnet(dataCox, survival_object, family = "cox", alpha = 1)  # alpha = 1 for LASSO, 0 for ridge

# plot the fit
plot(fit)
plot(fit, xlim = c(log(0.14), log(0.35)), ylim = c(9.1,10))

# Extract coefficients from the fitted model
CoxCoefficients <- extractCoxRegressionCoefficients(fit)

#---------------------------- Get Information about the genes whose coefficient is not zero

# Function to get gene information
gene_info_result <- getGenesInfo(CoxCoefficients)

#---------------------------- Extract patient ID 

# Function to extract patient ID to an Excel
extractSampleIDColumnToExcel(survival_data, "patient", "output_file.xlsx")

