# Function to calculate the C-index

calculate_c_index <- function(survival_object_train, genes_expression_train, models_coefficients, genes_expression_test, survival_test, survival_data) {
  
  # Fit a Cox regression model using the covariates
  fit <- coxph(survival_object_train ~ ., 
               data = subset(genes_expression_train, select = models_coefficients$ensembl_gene_id), # specify coefficients 
               init = as.numeric(models_coefficients$Coef_value), # specify coefficient values
               iter.max = 0) # force the software to keep those values
  
  # Construct a risk score based on the linear predictor on the test data
  survival_probabilities_test <- predict(fit, newdata = subset(genes_expression_test, select = models_coefficients$ensembl_gene_id), type ="lp")
  
  riskAUC = risksetAUC(Stime=survival_test$days,
                       status = survival_test$vital_status,
                       marker = survival_probabilities_test,
                       method = "Cox",
                       tmax = ceiling(max(survival_data$days)),
                       plot = FALSE)
  
  # Store c-index value for the current model
  c_index <-riskAUC$Cindex
  
  # Return the list of C-index values
  return(c_index)
}


# Get Information about the genes whose coefficient is not zero

getGenesInfo <- function(coefficients) {
  
  # Remove the version of the gene from its name(everything after the dot including the dot)
  coefficients$ensembl_gene_id <- sub("\\..*", "", coefficients$ensembl_gene_id) 
  
  # Connect to the Ensembl database
  #ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl") #mirror = "www")
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
  residuals_list <- residuals_list[!sapply(residuals_list, is.null)]
  
  # Transform in a matrix
  residuals_matrix <- do.call(cbind, residuals_list)
  
  return(residuals_matrix)
}

# Calculate rank product

calculate_rank_product <- function(data_matrix) {
  n_samples <- nrow(data_matrix)
  n_models <- ncol(data_matrix)
  
  ranks <- apply(data_matrix, 2, rank, ties.method = "first")
  rank_product <- apply(ranks, 1, prod)
  
  return(rank_product)
}
