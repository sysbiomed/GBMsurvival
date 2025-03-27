## Prepare survival data

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

## Prepare gene expression data

# Extrair dataframe com a expressão dos genes de interesse
createGeneExpressionDataframe <- function(genes_of_interest, expression_measure, source_data, survivalData) {


  genes_table <- as.data.frame(source_data@rowRanges@elementMetadata@listData) # Tranformar a informação sobre os genes num dataframe
  genes_table <- genes_table[genes_table$gene_type == genes_of_interest, ] #selecionar os genes de interesse

  expression_data <- as.data.frame(source_data@assays@data@listData[[expression_measure]]) #dataframe com a expressão de todos os genes
  colnames(expression_data) <- source_data@colData@listData[["patient"]] # adicionar os barcodes ao dataframe
  expression_data <- expression_data[, !duplicated(names(expression_data))] #remover individuos em duplicado
  rownames(expression_data) <- source_data@rowRanges@elementMetadata@listData$gene_id # adicionar o nome dos genes ao dataframe
  expression_data <- expression_data[,colnames(expression_data) %in% survivalData$patient] # selecionar apenas os individuos dos quais temos informação sobre sobrevivência
  expression_data <- as.data.frame(t(expression_data)) #transpor a dataframe

  data_frame <- expression_data[, colnames(expression_data) %in% genes_table$gene_id] #selecionar a expressão dos genes de interesse

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

## Select the disease and normalize the data

# selecionar apenas os dados da doença de interesse
selectDataPerDisease <- function(disease, expression_data, survivalData) {

  #survival data
  survivalData <- survivalData[survivalData$project_id %in% disease, ] # select the survival data from the patients with the disease
  survivalData <- survivalData[, c('patient', 'vital_status', 'days')] # keep only the columns of interest
  #expression data
  merged_data <- merge(survivalData, expression_data, by = "patient") # select only the patients from the gene expression dataframe with the desired disease and select only the patients whose the sum of the gene expression is not zero
  expression_data <- merged_data %>% dplyr::select(-vital_status, -days) # delete the columns added above

  result <- list(survival_data = survivalData, genes_expression = expression_data)

  return(result)
}

# EdgeR normalization
normalizationEdgeR <- function(data) {

  patients_IDs <- data[, 1]
  genes_expression <- t(data[, -1])
  dge <- DGEList(counts = genes_expression)
  dge <- normLibSizes(dge, method = "TMM")

  # usando o voom
  y <- voom(dge, plot=F)
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
  rownames(data) <- NULL

  return(data)
}


## Split the data in train and test

splitTestAndTrain <- function(expression, survival, percentage, seed_split) {

  merged_data <- merge(survival, expression, by = "patient")

  # # creating a sample diving into the ratio defined
  # sample <- sample.split(merged_data$patient, SplitRatio = percentage)
  #
  # # creating training dataset
  # train_data  <- subset(merged_data, sample == TRUE)
  #
  # # creating testing dataset
  # test_data <- subset(merged_data, sample == FALSE)

  # survival_train <- train_data[, c('patient', 'vital_status', 'days')]
  # survival_test <- test_data[, c('patient', 'vital_status', 'days')]
  # expression_train <- train_data %>% dplyr::select(-vital_status, -days)
  # expression_test <- test_data %>% dplyr::select(-vital_status, -days)

  # Use the `vital_status` as the strata variable for stratified sampling
  set.seed(seed_split)
  split <- initial_split(merged_data, prop = percentage, strata = "vital_status")

  # Extract train and test sets from the split
  survival_train <- training(split)[, c('patient', 'vital_status', 'days')]
  expression_train <- training(split) %>% dplyr::select(-vital_status, -days)
  survival_test <- testing(split)[, c('patient', 'vital_status', 'days')]
  expression_test <- testing(split) %>% dplyr::select(-vital_status, -days)

  result <- list(survival_train = survival_train, survival_test = survival_test,
                 expression_train = expression_train, expression_test = expression_test)

  return(result)
}
