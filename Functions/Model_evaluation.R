# Function to calculate the C-index

calculate_c_index_p_value <- function(survival_object_train, genes_expression_train, models_coefficients, genes_expression_test, survival_test, survival_data) {

  # Fit a Cox regression model using the covariates
  fit <- coxph(survival_object_train ~ .,
               data = subset(genes_expression_train, select = models_coefficients$ensembl_gene_id), # specify coefficients
               init = as.numeric(models_coefficients$Coef_value), # specify coefficient values
               iter.max = 0) # force the software to keep those values

  # Construct a risk score based on the linear predictor on the test data
  survival_probabilities_test <- predict(fit, newdata = subset(genes_expression_test, select = models_coefficients$ensembl_gene_id), type ="lp")

  c_index <- concordance(Surv(survival_test$days, survival_test$vital_status) ~ survival_probabilities_test, reverse=TRUE)$concordance

  # Categorize individuals of the test data based on the median of the train data
  survival_probabilities_train <- predict(fit, newdata = subset(genes_expression_train, select = models_coefficients$ensembl_gene_id), type ="lp")
  risk <- ifelse(survival_probabilities_test > median(survival_probabilities_train), "High", "Low")

  # P-value da Kaplan-Meier com a separação por High/ Low
  p_value <- survdiff(Surv(survival_test$days/365.25, survival_test$vital_status) ~ risk)$pvalue

  result <- list(c_index = c_index, p_value = p_value)

  # Return the list of C-index values
  return(result)
}


# Function to calculate the Integrated Brier Score (IBS)

calculate_integrated_brier_score <- function(survival_object_train, genes_expression_train, models_coefficients, genes_expression_test, survival_test) {

  # Fit a Cox regression model using the covariates
  fit <- coxph(survival_object_train ~ .,
               data = subset(genes_expression_train, select = models_coefficients$ensembl_gene_id), # specify coefficients
               init = as.numeric(models_coefficients$Coef_value), # specify coefficient values
               iter.max = 0, # force the software to keep those values
               x = TRUE)

  # Prepare test data by adding survival time and status as columns
  genes_expression_test$days <- survival_test$days
  genes_expression_test$vital_status <- survival_test$vital_status

  # Calculate the Brier Score over time and the Integrated Brier Score (IBS)
  ibs_result <- pec(object = list("Cox" = fit),
                    formula = Surv(days, vital_status) ~ 1,
                    data = genes_expression_test,
                    times = sort(unique(survival_test$days)), # time points for evaluation
                    start = min(survival_test$days),
                    maxtime = max(survival_test$days))

  # Extract the Integrated Brier Score
  ibs <- crps(ibs_result)[2]

  return(ibs)
}


# Function to calculate Mean Reciprocal Rank (MRR) for each model
calculate_MRR <- function(metric_df, higher_is_better = TRUE) { #TRUE for C-index, FALSE for IBS

  #notas: para a definição do rank os NA é como se não existissem e para o calculo do MRR não entram na média.
  #Exemplo (menor é melhor): values: 3, NA, 5, 2, NA, 4; rank: 2, NA, 4, 1, NA, 3

  # Determine the ranks based on the `higher_is_better` parameter
  ranked_df <- if (higher_is_better) {
    # Higher values better: rank in descending order
    apply(metric_df, 2, function(x) rank(-x, ties.method = "average", na.last = "keep")) # 2 means it operates over columns
  } else {
    # Lower values better: rank in ascending order
    apply(metric_df, 2, function(x) rank(x, ties.method = "average", na.last = "keep"))
  }

  # Find the reciprocal rank of each model in each seed
  reciprocal_ranks <- 1 / ranked_df

  # Calculate the Mean Reciprocal Rank (MRR) for each model, , excluding NAs
  mrr_values <- rowMeans(reciprocal_ranks, na.rm = TRUE)

  # Create a data frame with MRR results
  mrr_results <- data.frame(MRR = mrr_values)

  return(mrr_results)
}

# Get Information about the genes whose coefficient is not zero

getGenesInfo <- function(coefficients) {

  # Remove the version of the gene from its name(everything after the dot including the dot)
  coefficients$ensembl_gene_id <- sub("\\..*", "", coefficients$ensembl_gene_id)

  # Connect to the Ensembl database
  ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "www")
  # Get information about the genes
  gene_info <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description"),
                     filters = "ensembl_gene_id",
                     values = coefficients$ensembl_gene_id,
                     mart = ensembl)

  # Merge gene coefficents and gene information
  gene_info <- merge(coefficients, gene_info, by = "ensembl_gene_id")
  gene_info <- gene_info[order(-abs(gene_info$Coef_value)), ]

  return(gene_info)
}

# Cumulative case/dynamic control ROC

# fonte: https://datascienceplus.com/time-dependent-roc-for-survival-prediction-models-in-r/
# The cumulative sensitivity considers those who have died by time t
# The dynamic specificity regards those who are still alive at time t
# This is one gives the prediction performance for the risk (cumulative incidence) of events over the t-year period.

CumulativeCaseDynamicControlROC <- function(survival_data, survival_probabilities, title) {

  ## Define a helper function to evaluate at various t
  survivalROC_helper <- function(t) {
    survivalROC(Stime        = survival_data$days/365.25,
                status       = survival_data$vital_status,
                marker       = survival_probabilities,
                predict.time = t,
                method       = "NNE",
                span = 0.25 * sum(survival_data$vital_status)^(-0.20))
  }
  ## Evaluate every X days
  survivalROC_data <- data_frame(t = c(1,2,3,4,5,6)) %>%
    mutate(survivalROC = purrr::map(t, survivalROC_helper),
           ## Extract scalar AUC
           auc = purrr::map_dbl(survivalROC, magrittr::extract2, "AUC"),
           ## Put cut off dependent values in a data_frame
           df_survivalROC = purrr::map(survivalROC, function(obj) {
             as_tibble(obj[c("cut.values","TP","FP")])
           })) %>%
    dplyr::select(-survivalROC) %>%
    tidyr::unnest() %>%
    arrange(t, FP, TP)

  ## Define custom labels for each time point
  custom_labels <- c("1" = "Year 1",
                     "2" = "Year 2",
                     "3" = "Year 3",
                     "4" = "Year 4",
                     "5" = "Year 5",
                     "6" = "Year 6")

  ## Plot
  survivalROC_data %>%
    ggplot(mapping = aes(x = FP, y = TP)) +
    geom_point() +
    geom_line() +
    geom_label(data = survivalROC_data %>% dplyr::select(t,auc) %>% unique,
               mapping = aes(label = sprintf("%.3f", auc)), x = 0.5, y = 0.5) +
    facet_wrap(~ t, labeller = labeller(t = custom_labels)) +  # Apply custom labels
    theme_bw() +
    ggtitle(title) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          legend.key = element_blank(),
          plot.title = element_text(hjust = 0.5),
          strip.background = element_blank())
}

# Incident case/dynamic control ROC

# fonte: https://datascienceplus.com/time-dependent-roc-for-survival-prediction-models-in-r/
# The incident sensitivity considers those who die at time t
# The dynamic specificity regards those who are still alive at time t
# This one gives the prediction performance for the hazard (incidence in the risk set) of events at t-year among those who are in the risk set at t.

IncidentCaseDynamicControlROC <- function(survival_data, survival_probabilities, title) {

  ## Define a helper function to evaluate at various t
  risksetROC_helper <- function(t) {
    risksetROC(Stime        = survival_data$days/365.25,
               status       = survival_data$vital_status,
               marker       = survival_probabilities,
               predict.time = t,
               method       = "Cox",
               plot         = FALSE)
  }
  ## Evaluate every year
  risksetROC_data <- data_frame(t = c(1,2,3,4,5,6)) %>%
     mutate(risksetROC = purrr::map(t, risksetROC_helper),
           ## Extract scalar AUC
           auc = purrr::map_dbl(risksetROC, magrittr::extract2, "AUC"),
           ## Put cut off dependent values in a data_frame
           df_risksetROC = purrr::map(risksetROC, function(obj) {
             ## marker column is too short!
             marker <- c(-Inf, obj[["marker"]], Inf)
             bind_cols(data_frame(marker = marker),
                       as_tibble(obj[c("TP","FP")]))
           })) %>%
    dplyr::select(-risksetROC) %>%
    tidyr::unnest() %>%
    arrange(t, FP, TP)

  ## Define custom labels for each time point
  custom_labels <- c("1" = "Year 1",
                     "2" = "Year 2",
                     "3" = "Year 3",
                     "4" = "Year 4",
                     "5" = "Year 5",
                     "6" = "Year 6")

  ## Plot
  risksetROC_data %>%
    ggplot(mapping = aes(x = FP, y = TP)) +
    geom_point() +
    geom_line() +
    geom_label(data = risksetROC_data %>% dplyr::select(t,auc) %>% unique,
               mapping = aes(label = sprintf("%.3f", auc)), x = 0.5, y = 0.5) +
    facet_wrap(~ t, labeller = labeller(t = custom_labels)) +  # Apply custom labels
    theme_bw() +
    ggtitle(title) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          legend.key = element_blank(),
          plot.title = element_text(hjust = 0.5),
          strip.background = element_blank())
}



# Calculate martingale residuals

calculate_martingale_residuals <- function(models_coefficients_list, genes_expression, survival_data) {

  survival_object_all <- Surv(time = survival_data$days, event = survival_data$vital_status)
  residuals_list <- list()

  models_coefficients <- unlist(models_coefficients_list, recursive = FALSE)

  for (i in seq_along(models_coefficients)) {
    if (length(models_coefficients[[i]]$ensembl_gene_id) == 0) {
      residuals_list[[i]] <- NULL
      next
    }

    fit <- coxph(survival_object_all ~ .,
                 data = subset(genes_expression, select = models_coefficients[[i]]$ensembl_gene_id),
                 init = as.numeric(models_coefficients[[i]]$Coef_value),
                 iter.max = 0)

    martingale_residuals <- residuals(fit, type = "martingale")
    residuals_list[[i]] <- abs(martingale_residuals)
  }

  # Remove NULL elements
  residuals_list <- residuals_list[!sapply(residuals_list, is.null)] #These NULL elements are cases where the model had no coefficients

  # Transform in a matrix
  residuals_matrix <- do.call(cbind, residuals_list)

  return(residuals_matrix)
}

# Calculate rank product

calculate_rank_product <- function(data_matrix) {

  ranks <- apply(data_matrix, 2, rank, ties.method = "first")
  rank_product <- apply(ranks, 1, prod)

  return(rank_product)
}

# Make a df with the genes that are selected in more than one seed

get_gene_counts_with_names <- function(genes_list) {
  # Flatten the list of all genes and count their occurrences
  all_genes <- unlist(genes_list)
  gene_counts <- table(all_genes)

  # Filter genes that appear in more than one seed
  genes_in_multiple_seeds <- names(gene_counts[gene_counts > 1])
  counts_in_multiple_seeds <- as.numeric(gene_counts[gene_counts > 1])

  # Create a data frame to display the genes and their counts
  df_genes_in_multiple_seeds <- data.frame(Gene = genes_in_multiple_seeds, Count = counts_in_multiple_seeds)

  # Initialize the connection to the Ensembl database
  ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

  # List of Ensembl gene IDs (remove the version numbers)
  gene_ids <- gsub("\\..*", "", df_genes_in_multiple_seeds$Gene)

  # Get information about the genes
  gene_info <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                     filters = "ensembl_gene_id",
                     values = gene_ids,
                     mart = ensembl)

  # Merge the original data frame with gene symbols
  df_genes_in_multiple_seeds$Gene_Name <- gene_info$external_gene_name[match(gene_ids, gene_info$ensembl_gene_id)]

  # Order the data frame by the 'Count' column in descending order
  df_genes_in_multiple_seeds <- df_genes_in_multiple_seeds[order(-df_genes_in_multiple_seeds$Count), ]

  # Select only the columns to "Gene_Name" and "Count"
  df_genes_in_multiple_seeds <- df_genes_in_multiple_seeds[, c("Gene_Name", "Count")]

  # Return the resulting data frame
  return(df_genes_in_multiple_seeds)
}


