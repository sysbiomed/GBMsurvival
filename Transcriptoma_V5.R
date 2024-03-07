# Required libraries:
library(TCGAbiolinks)
library(dplyr)
library(tibble)
library(tidyverse)
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
library(pROC)
library(caTools)
library(survivalROC)
library(risksetROC)
library(lattice)
library(survMisc)

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

# EdgeR normalization
normalizationEdgeR <- function(data) { 
  
  patients_IDs <- data[, 1]
  genes_expression <- t(data[, -1])
  dge <- DGEList(counts = genes_expression)
  dge <- normLibSizes(dge, method = "TMM")
  
  # usando cpm
  #normcounts <- cpm(dge)
  
  # usando o voom
  y <- voom(dge, plot=T)
  lognormcounts <- as.data.frame(y$E)
  normcounts <- as.data.frame(2^lognormcounts)
  
  genes_expression <- t(normcounts)
  genes_expression <- as.data.frame(genes_expression)
  rownames(genes_expression) <- NULL
  genes_expression <- cbind(patient = patients_IDs, genes_expression)

  return(genes_expression)
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

#---------------------------- Split the data in train and test

# Split the data in train and test
splitTestAndTrain <- function(expression, survival, percentage) { 
  
  merged_data <- merge(survival, expression, by = "patient")
  
  # creating a sample diving into the ratio defined 
  sample <- sample.split(merged_data$patient, SplitRatio = percentage)
  
  # creating training dataset 
  train_data  <- subset(merged_data, sample == TRUE) 
  
  # creating testing dataset 
  test_data <- subset(merged_data, sample == FALSE) 
  
  # Separate survival from expression
  survival_train <- train_data[, c(1:3)]
  survival_test <- test_data[, c(1:3)]
  expression_train <- train_data[, -c(2:3)]
  expression_test <- test_data[, -c(2:3)]
  
  result <- list(survival_train = survival_train, survival_test = survival_test,
                 expression_train = expression_train, expression_test = expression_test)
  
  return(result)
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
  #rownames(coefficients) <- sub("\\..*", "", rownames(coefficients)) # Remove the version of the gene from its name(everything after the dot including the dot)
  coefficients$ensembl_gene_id <- rownames(coefficients)
  rownames(coefficients) <- NULL
  
  # Print min lambda
  #cat("Lambda min:", "\n", cox_fit$lambda.min)

  return(coefficients)
}


#---------------------------- Get Information about the genes whose coefficient is not zero

# Function to get gene information
getGenesInfo <- function(coefficients) {
  
  # Remove the version of the gene from its name(everything after the dot including the dot)
  coefficients$ensembl_gene_id <- sub("\\..*", "", coefficients$ensembl_gene_id) 
  
  # Connect to the Ensembl database
  #ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl") #mirror = "www"
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

#---------------------------- Cumulative case/dynamic control ROC, fonte: https://datascienceplus.com/time-dependent-roc-for-survival-prediction-models-in-r/
# The cumulative sensitivity considers those who have died by time t
# The dynamic specificity regards those who are still alive at time t 
# This is one gives the prediction performance for the risk (cumulative incidence) of events over the t-year period.

CumulativeCaseDynamicControlROC <- function(survival_data, survival_probabilities) {
  
  ## Define a helper function to evaluate at various t
  survivalROC_helper <- function(t) {
    survivalROC(Stime        = survival_data$days,
                status       = survival_data$vital_status,
                marker       = survival_probabilities,
                predict.time = t,
                method       = "NNE",
                span = 0.25 * nrow(survival_data)^(-0.20))
  }
  ## Evaluate every X days
  survivalROC_data <- data_frame(t = ceiling(max(survival_data$days)/6) * c(1,2,3,4,5,6)) %>%
    mutate(survivalROC = map(t, survivalROC_helper),
           ## Extract scalar AUC
           auc = map_dbl(survivalROC, magrittr::extract2, "AUC"),
           ## Put cut off dependent values in a data_frame
           df_survivalROC = map(survivalROC, function(obj) {
             as_data_frame(obj[c("cut.values","TP","FP")])
           })) %>%
    dplyr::select(-survivalROC) %>%
    unnest() %>%
    arrange(t, FP, TP)
  ## Plot
  survivalROC_data %>%
    ggplot(mapping = aes(x = FP, y = TP)) +
    geom_point() +
    geom_line() +
    geom_label(data = survivalROC_data %>% dplyr::select(t,auc) %>% unique,
               mapping = aes(label = sprintf("%.3f", auc)), x = 0.5, y = 0.5) +
    facet_wrap( ~ t) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          legend.key = element_blank(),
          plot.title = element_text(hjust = 0.5),
          strip.background = element_blank())
}

#---------------------------- Incident case/dynamic control ROC, fonte: https://datascienceplus.com/time-dependent-roc-for-survival-prediction-models-in-r/
# The incident sensitivity considers those who die at time t
# The dynamic specificity regards those who are still alive at time t 
# This one gives the prediction performance for the hazard (incidence in the risk set) of events at t-year among those who are in the risk set at t.

IncidentCaseDynamicControlROC <- function(survival_data, survival_probabilities) {
  
  ## Define a helper function to evaluate at various t
  risksetROC_helper <- function(t) {
    risksetROC(Stime        = survival_data$days,
               status       = survival_data$vital_status,
               marker       = survival_probabilities,
               predict.time = t,
               method       = "Cox",
               plot         = FALSE)
  }
  ## Evaluate every 180 days
  risksetROC_data <- data_frame(t = ceiling(max(survival_data$days)/6) * c(1,2,3,4,5,6)) %>%
    mutate(risksetROC = map(t, risksetROC_helper),
           ## Extract scalar AUC
           auc = map_dbl(risksetROC, magrittr::extract2, "AUC"),
           ## Put cut off dependent values in a data_frame
           df_risksetROC = map(risksetROC, function(obj) {
             ## marker column is too short!
             marker <- c(-Inf, obj[["marker"]], Inf)
             bind_cols(data_frame(marker = marker),
                       as_data_frame(obj[c("TP","FP")]))
           })) %>%
    dplyr::select(-risksetROC) %>%
    unnest() %>%
    arrange(t, FP, TP)
  ## Plot
  risksetROC_data %>%
    ggplot(mapping = aes(x = FP, y = TP)) +
    geom_point() +
    geom_line() +
    geom_label(data = risksetROC_data %>% dplyr::select(t,auc) %>% unique,
               mapping = aes(label = sprintf("%.3f", auc)), x = 0.5, y = 0.5) +
    facet_wrap( ~ t) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          legend.key = element_blank(),
          plot.title = element_text(hjust = 0.5),
          strip.background = element_blank())
}



#---------------------------- Find the best threshold

# Find the best threshold
# findBestThreshold <- function(survival_train, survival_probabilities) { 
#   
#   # Create ROC curve
#   roc_curve <- roc(survival_train$vital_status, survival_probabilities)
#   
#   # Plot ROC curve and add threshold to the ROC curve plot
#   par(pty="s")
#   plot(roc_curve, main = "ROC Curve", col = "blue", lwd = 2, print.auc=TRUE,
#        print.thres=TRUE, legacy.axes= TRUE, add= FALSE, xlim = c(1,0))
#   
#   # Find threshold based on Youden's J statistic 
#   best_threshold <- coords(roc_curve, "best")$threshold
#   
#   return(best_threshold)
# }

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

#---------------------------- Plot da curva de sobrevivência

# Kaplan-Meier com a separação por doença
# fit_disease_RNA <- survfit(Surv(survival_data$days, survival_data$vital_status) ~ survival_data$project_id)
# survdiff(Surv(survival_data$days, survival_data$vital_status) ~ survival_data$project_id)
# ggsurvplot(fit_disease_RNA, data = survival_data)

#---------------------------- Tratar os dados de expressão génica

# Extrair dataframe com a expressão dos genes de interesse
genes_expression <- createGeneExpressionDataframe("protein_coding", "unstranded", TCGA_RNA, survival_data)
#genes_expression <- createGeneExpressionDataframe("lncRNA", "unstranded", TCGA_RNA, survival_data)       #<------------ indicar os dados

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

# EdgeR+voom normalization
genes_expression <- normalizationEdgeR(genes_expression)

# Clean the gene expression dataframe
genes_expression <- cleanGeneExpressionData(genes_expression)

# ordenar os dados por paciente
genes_expression <- orderDataByPatient(genes_expression)
survival_data <- orderDataByPatient(survival_data)

# create a dataframe with categorized gene expression (high and low)
#genes_expression <- categorizeGeneExpressionData(genes_expression)                                        #<------------ indicar se queremos os dados categorizados

#---------------------------- Find the best alpha

# List of alpha values
alpha_values <- seq(0.1, 1, by = 0.1)

# List of seed values
seed_values <- seq(1000, 1100, by = 10)


# Initialize an empty dataframe to store results
results_df <- data.frame(matrix(nrow = length(alpha_values), ncol = 2))
colnames(results_df) <- c("Alpha", "Average_C_Index")

# Create an empty list to store c-index values for each alpha
c_index_list <- vector("list", length(alpha_values))

# Iterate over each alpha value
for (i in seq_along(alpha_values)) {
  
  cat("alpha value: ", "\n", alpha_values[i], "\n")
  
  # Store c-index values for each seed
  c_index_values <- numeric(length(seed_values))

  # Iterate over each seed value
  for (j in seq_along(seed_values)) {

    cat("seed value: ", "\n", seed_values[j], "\n")
    
    # Split the data in train and test
    set.seed(seed_values[j])
    splited <- splitTestAndTrain(genes_expression, survival_data, 0.7)
    
    genes_expression_train <- splited$expression_train
    survival_train <- splited$survival_train
    genes_expression_test <- splited$expression_test
    survival_test <- splited$survival_test

    # Define survival object
    survival_object <- Surv(time = survival_train$days, event = survival_train$vital_status)
    
    # fit the Regularized Cox Regression
    set.seed(1012)
    fit <- cv.glmnet(prepareDataForCoxRegression(genes_expression_train),
                     survival_object, 
                     family = "cox",
                     type.measure = "C",
                     alpha = alpha_values[i])  # alpha = 1 for LASSO, 0 for ridge
    
    # Extract coefficients from the fitted model
    CoxCoefficients <- extractCoxRegressionCoefficients(fit)
    
    if (length(CoxCoefficients$ensembl_gene_id) == 0) {
      # If no coefficients found, set c-index to 0.5 and skip to the next iteration,
      cat("No coefficients found or only one coefficient found", "\n")
      c_index_values[j] <- 0.5
      # Store c-index values for the current alpha in the list
      c_index_list[[i]] <- c_index_values
      next
    }
    
    print(CoxCoefficients)
    
    # Fit a Cox regression model using the covariates
    fit <- coxph(survival_object ~ .,
                 data = subset(genes_expression_train, select = CoxCoefficients$ensembl_gene_id), # specify coefficients
                 init = as.numeric(CoxCoefficients$Coef_value), # specify coefficient values
                 iter.max = 0) # force the software to keep those values
    
    #---------------------------- Test the coefficients (predictor) found
    # Construct a risk score based on the linear predictor on the test data
    survival_probabilities_test <- predict(fit, newdata = subset(genes_expression_test, select = CoxCoefficients$ensembl_gene_id), type = "lp")
   
    # Plot AUC based on incident/dynamic definition
    riskAUC = risksetAUC(Stime=survival_test$days,
                      status = survival_test$vital_status,
                      marker = survival_probabilities_test,
                      method = "Cox",
                      tmax = ceiling(max(survival_data$days)),
                      plot = FALSE)
    
    # Store c-index value for the current seed, explicação: https://www.youtube.com/watch?v=rRYfWAsG4RI
    c_index_values[j] <- max(0.5, riskAUC$Cindex) # Ensure c-index is not less than 0.5
  
  }
  
  # Store c-index values for the current alpha in the list
  c_index_list[[i]] <- c_index_values
  
  # Calculate average c-index for the current alpha
  avg_c_index <- mean(c_index_values)

  # Store results in dataframe
  results_df[i, "Alpha"] <- alpha_values[i]
  results_df[i, "Average_C_Index"] <- avg_c_index
}

# Add c-index list to results_df
results_df$c_index_values <- c_index_list

# Print the results dataframe
print(results_df)

#---------------------------- Test the best alpha

# Extract the best alpha value
best_alpha <- results_df$Alpha[which.max(results_df$Average_C_Index)]

# Split the data in train and test
set.seed(1012)
splited <- splitTestAndTrain(genes_expression, survival_data, 0.7)

genes_expression_train <- splited$expression_train
survival_train <- splited$survival_train
genes_expression_test <- splited$expression_test
survival_test <- splited$survival_test

# Define survival object
survival_object <- Surv(time = survival_train$days, event = survival_train$vital_status)

# fit the Regularized Cox Regression
set.seed(1012)
fit <- cv.glmnet(prepareDataForCoxRegression(genes_expression_train),
                 survival_object, 
                 family = "cox",
                 type.measure = "C",
                 alpha = best_alpha)  # alpha = 1 for LASSO, 0 for ridge

# plot the fit
plot(fit)

# Extract coefficients from the fitted model
CoxCoefficients <- extractCoxRegressionCoefficients(fit)

# Function to get gene information
gene_info_result <- getGenesInfo(CoxCoefficients)

# Fit a Cox regression model using the covariates
fit <- coxph(survival_object ~ ., 
             data = subset(genes_expression_train, select = CoxCoefficients$ensembl_gene_id), # specify coefficients 
             init = as.numeric(CoxCoefficients$Coef_value), # specify coefficient values
             iter.max = 0) # force the software to keep those values

# Test the proportional hazards assumption for a Cox regression model fit (coxph)
Propo_hazards <- cox.zph(fit, transform="km", terms=TRUE, singledf=FALSE, global=TRUE)
Propo_hazards_global_p_value <- Propo_hazards$table["GLOBAL", "p"]

#---------------------------- Test the coefficients (predictor) found

# Construct a risk score based on the linear predictor on the test data
survival_probabilities_test <- predict(fit, newdata = subset(genes_expression_test, select = CoxCoefficients$ensembl_gene_id), type = "lp")

#---------------------------- Cumulative case/dynamic control ROC, fonte: https://datascienceplus.com/time-dependent-roc-for-survival-prediction-models-in-r/
# The cumulative sensitivity considers those who have died by time t
# The dynamic specificity regards those who are still alive at time t 
# This is one gives the prediction performance for the risk (cumulative incidence) of events over the t-year period.
CumulativeCaseDynamicControlROC(survival_test, survival_probabilities_test)

#---------------------------- Incident case/dynamic control ROC, fonte: https://datascienceplus.com/time-dependent-roc-for-survival-prediction-models-in-r/
# The incident sensitivity considers those who die at time t
# The dynamic specificity regards those who are still alive at time t 
# This one gives the prediction performance for the hazard (incidence in the risk set) of events at t-year among those who are in the risk set at t.
IncidentCaseDynamicControlROC(survival_test, survival_probabilities_test)

#---------------------------- Plot AUC based on incident/dynamic definition

# Plot AUC
riskAUC = risksetAUC(Stime=survival_test$days,
                     status = survival_test$vital_status,
                     marker = survival_probabilities_test,
                     method = "Cox",
                     tmax = ceiling(max(survival_data$days)),
                     plot = TRUE)

# Store c-index value for the current seed, explicação: https://www.youtube.com/watch?v=rRYfWAsG4RI
c_index_test <- riskAUC$Cindex

#---------------------------- Find the best cutoff value on the train data

# Construct a risk score based on the linear predictor on the train data
survival_probabilities_train <- predict(fit, newdata = subset(genes_expression_train, select = CoxCoefficients$ensembl_gene_id), type = "lp")

# Find the best threshold
find_best_cutoff <- cutp(survfit(Surv(survival_train$days, survival_train$vital_status) ~ survival_probabilities_train))
best_threshold <- find_best_cutoff$survival_probabilities_train$survival_probabilities_train[which.max(find_best_cutoff$survival_probabilities_train$U)]

#---------------------------- Apply the best cutoff value on the test data

# Categorize individuals of the test data based on the threshold
risk_groups <- ifelse(survival_probabilities_test > best_threshold, "High", "Low")

# Kaplan-Meier com a separação por High/ Low
fit_surv <- survfit(Surv(survival_test$days, survival_test$vital_status) ~ risk_groups)
survdiff(Surv(survival_test$days, survival_test$vital_status) ~ risk_groups)
ggsurvplot(fit_surv, data = survival_test)

#---------------------------- Extract patient ID 

# Function to extract patient ID to an Excel
#extractSampleIDColumnToExcel(survival_data, "patient", "output_file.xlsx")

