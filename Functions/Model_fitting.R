## Regularized Cox Regression functions


#univariate Cox regression for each variable
univariate_cox <- function(survivaldata, expressiondata, significant_threshold) {

  results <- data.frame()

  # Define survival object
  survival_object <- Surv(time = survivaldata$days, event = survivaldata$vital_status)

  # Define the covariates to study
  covariates <- colnames(expressiondata)[colnames(expressiondata) != "patient"]

  for (cov in covariates) {

    # Fit a Cox regression model using the covariates
    set.seed(1)
    uni_cox <- coxph(survival_object ~ ., data = subset(expressiondata, select = cov))

    #summary_uni_cox <- summary(uni_cox)

    #p_value <- summary_uni_cox$coefficients[5]  # p-value

    # Construct a risk score based on the linear predictor on the test data
    survival_probabilities <- predict(uni_cox, newdata = subset(expressiondata, select = cov), type = "lp")

    # Categorize individuals of the test data based on the median
    risk <- ifelse(survival_probabilities > median(survival_probabilities), "High", "Low")

    # P-value da Kaplan-Meier com a separação por High/ Low
    p_value <- survdiff(Surv(survivaldata$days/365.25, survivaldata$vital_status) ~ risk)$pvalue

    # Adicionar o p-value desta variável ao dataframe
    results <- rbind(results, data.frame(variable = cov, p_value = p_value))

  }

  # Extract the names of the significant variables
  significant_variable_names <- results$variable[results$p_value < significant_threshold]

  return(significant_variable_names)
}



#preparar os dados para a Regularized Cox Regression
prepareDataForCoxRegression <- function(data) {

  rownames(data) <- data$patient
  data <- data %>% dplyr::select(-patient)
  data <- data.matrix(data)

  return(data)
}

# Extract coefficients from the fitted model
extractCoxRegressionCoefficients <- function(cox_fit) {

  coefficients <- data.frame(as.matrix(coef(cox_fit, s = "lambda.min")))  # Choose the lambda that minimizes cross-validated error
  colnames(coefficients)[colnames(coefficients) == 'X1'] <- "Coef_value" # Change the column name
  coefficients <- subset(coefficients, Coef_value != 0) # select only the coefficients that are not equal to zero
  coefficients$ensembl_gene_id <- rownames(coefficients)
  rownames(coefficients) <- NULL

  return(coefficients)
}


## Find the best alpha for Regularized Cox Regression model glmnet

find_best_alpha_for_glmnet <- function(alpha_vector, seed_values, genes_expression, survival_data) {

  # Initialize an empty dataframe to store results
  results_df <- data.frame(matrix(nrow = length(alpha_vector), ncol = 2))
  colnames(results_df) <- c("Alpha", "Average_C_Index")

  # Create an empty list to store c-index values for each alpha
  c_index_list <- vector("list", length(alpha_vector))

  # Iterate over each alpha value
  for (i in seq_along(alpha_vector)) {
    cat("alpha: ",alpha_vector[i], "\n")
    # Store c-index values for each seed
    c_index_values <- numeric(length(seed_values))

    # Iterate over each seed value
    for (j in seq_along(seed_values)) {
      cat("seed: ",seed_values[j], "\n")
      # Split the data in train and test
      set.seed(seed_values[j])
      splited <- splitTestAndTrain(genes_expression, survival_data, 0.7)

      genes_expression_train <- splited$expression_train
      survival_train <- splited$survival_train
      genes_expression_test <- splited$expression_test
      survival_test <- splited$survival_test

      # Run univariate Cox regression for each variable
      significant_variable_names <- univariate_cox(survival_train, genes_expression_train, 0.05)

      # reduce the genes expression data frame to only the colunmns of interest
      genes_expression_train <- genes_expression_train[, c("patient", significant_variable_names)]
      print(ncol(genes_expression_train))

      # Define survival object
      survival_object <- Surv(time = survival_train$days, event = survival_train$vital_status)

      # fit the Regularized Cox Regression GLMNET
      set.seed(1012)
      fit <- cv.glmnet(prepareDataForCoxRegression(genes_expression_train),
                       survival_object,
                       family = "cox",
                       #type.measure = "C",
                       alpha = alpha_vector[i])  # alpha = 1 for LASSO, 0 for ridge

      # Extract coefficients from the fitted model
      CoxCoefficients <- extractCoxRegressionCoefficients(fit)

      if (length(CoxCoefficients$ensembl_gene_id) == 0) {
        # If no coefficients found, set c-index to NA and skip to the next iteration,
        c_index_values[j] <- NA

        next
      }


      # Fit a Cox regression model using the covariates
      fit <- coxph(survival_object ~ .,
                   data = subset(genes_expression_train, select = CoxCoefficients$ensembl_gene_id), # specify coefficients
                   init = as.numeric(CoxCoefficients$Coef_value), # specify coefficient values
                   iter.max = 0) # force the software to keep those values

      # Construct a risk score based on the linear predictor on the test data
      survival_probabilities_test <- predict(fit, newdata = subset(genes_expression_test, select = CoxCoefficients$ensembl_gene_id), type = "lp")



      # Store c-index value for the current seed, explicação: https://www.youtube.com/watch?v=rRYfWAsG4RI
      c_index_values[j] <- concordance(Surv(survival_test$days, survival_test$vital_status) ~ survival_probabilities_test, reverse=TRUE)$concordance

      # c_index_values[j] <- concordance.index(x = survival_probabilities_test,
      #                                        surv.time = survival_test$days,
      #                                        surv.event = survival_test$vital_status,
      #                                        method = "noether")$c.index

      # # Plot AUC based on incident/dynamic definition
      # riskAUC = risksetAUC(Stime=survival_test$days,
      #                      status = survival_test$vital_status,
      #                      marker = survival_probabilities_test,
      #                      method = "Cox",
      #                      tmax = ceiling(max(survival_data$days)),
      #                      plot = FALSE)
      #
      # c_index_values[j] <- riskAUC$Cindex

      # c_index_values[j] <- summary(fit)$concordance[1]

    }

    # Store c-index values for the current alpha in the list
    c_index_list[[i]] <- c_index_values

    # Calculate average c-index for the current alpha
    avg_c_index <- mean(c_index_values, na.rm = TRUE)

    # Store results in dataframe
    results_df[i, "Alpha"] <- alpha_vector[i]
    results_df[i, "Average_C_Index"] <- avg_c_index
  }

  # Add c-index list to results_df
  results_df$c_index_values <- c_index_list

  return(results_df)
}


## Random Survival Forest Function to analyse different thresholds

fit_RSF <- function(genes_expression_train, survival_object, best_alpha, survival_train, thresholds_list, importance_rsf) {
  models_coefficients <- list()

  for (k in thresholds_list) {
    # Select the top features based on importance
    top_features_indices_rsf <- which(importance_rsf >= k)

    z <- as.matrix(genes_expression_train[, -which(names(genes_expression_train) %in% c("patient"))])
    top_features_rsf <- colnames(z)[top_features_indices_rsf]
    top_features_rsf <- top_features_rsf[!is.na(top_features_rsf)]

    # Fit Cox with selected features
    set.seed(1012)
    fit_rsf_model <- cv.glmnet(as.matrix(subset(genes_expression_train, select = top_features_rsf)),
                               survival_object,
                               family = "cox",
                               alpha = best_alpha)

    RandomForest_Coefficients <- extractCoxRegressionCoefficients(fit_rsf_model)
    models_coefficients[[as.character(k)]] <- RandomForest_Coefficients
  }

  return(list(coefficients = models_coefficients))
}


## Boruta Function to analyse different thresholds

fit_Boruta <- function(genes_expression_train, survival_object, best_alpha, survival_train, thresholds_list) {

  ############# Fit the Boruta model

  models_coefficients <- list()

  for (k in thresholds_list) {

  set.seed(1012) # For reproducibility
  boruta_output <- Boruta(survival_object ~ .,
                          data = as.data.frame(prepareDataForCoxRegression(genes_expression_train)),
                          doTrace = 2)
                          #maxRuns = 600,
                          #pValue = k) #0.0001)

  # Get the significant features
  significant_features <- getSelectedAttributes(boruta_output,
                                                withTentative = k)


  # Fit a Cox regression model to avoid having predictors that are highly correlated with others
  set.seed(1012)
  fit_boruta_model <- cv.glmnet(as.matrix(subset(genes_expression_train, select = significant_features)),
                                survival_object,
                                family = "cox",
                                alpha = best_alpha)

  # Extract coefficients from the fitted model
  Boruta_Coefficients <- extractCoxRegressionCoefficients(fit_boruta_model)


  models_coefficients[[as.character(k)]] <- Boruta_Coefficients
  }


# Return the dataframe and the fitted models
return(models_coefficients)

}



## Fit models

fit_models <- function(genes_expression_train, survival_object, best_alpha, survival_train) {

  ############# Fit the Regularized Cox Regression GLMNET with the best alpha
  set.seed(1012)
  fit_COX <- cv.glmnet(prepareDataForCoxRegression(genes_expression_train),
                       survival_object,
                       family = "cox",
                       #type.measure = "C",
                       alpha = best_alpha)  # alpha = 1 for LASSO, 0 for ridge

  # Extract coefficients from the fitted model
  Cox_Coefficients <- extractCoxRegressionCoefficients(fit_COX)

  # ############# fit the Regularized Cox Regression SIS
  # set.seed(1012)
  # fit_SIS <- SIS(prepareDataForCoxRegression(genes_expression_train),
  #                survival_object,
  #                family = "cox",
  #                penalty='lasso', # Cox model currently not implemented with the 'SCAD' (the default), 'MCP', or 'lasso' are provided
  #                tune='cv',
  #                #type.measure='deviance', #For penalty='SCAD' and penalty='MCP', only type.measure='deviance' is available
  #                iter = FALSE, #Specifies whether to perform iterative SIS.
  #                nsis = 100,
  #                seed = 1000)
  #
  #
  # # Extract coefficients from the fitted model
  # SIS_Coefficients <- data.frame(Coef_value = as.data.frame(fit_SIS$coef.est)[,1],
  #                                ensembl_gene_id = colnames(prepareDataForCoxRegression(genes_expression_train)[, fit_SIS$ix, drop = FALSE]))
  #
  # ############# fit the Regularized Cox Regression (I)SIS
  # set.seed(1012)
  # invisible(capture.output({fit_ISIS <- SIS(prepareDataForCoxRegression(genes_expression_train),
  #                                           survival_object,
  #                                           family = "cox",
  #                                           penalty='lasso', # Cox model currently not implemented with the 'SCAD' (the default), 'MCP', or 'lasso' are provided
  #                                           tune='cv',
  #                                           #type.measure='deviance', #For penalty='SCAD' and penalty='MCP', only type.measure='deviance' is available
  #                                           iter = TRUE, #Specifies whether to perform iterative SIS.
  #                                           nsis = 25,
  #                                           seed = 100,
  #                                           varISIS='cons') #Specifies whether to perform any of the two ISIS variants based on randomly splitting the sample into two groups
  #
  # }))
  #
  # # Extract coefficients from the fitted model
  # ISIS_Coefficients <- data.frame(Coef_value = as.data.frame(fit_ISIS$coef.est)[,1],
  #                                 ensembl_gene_id = colnames(prepareDataForCoxRegression(genes_expression_train)[, fit_ISIS$ix, drop = FALSE]))
  #

  ############# Fit the Fast Unified Random Forest

  # Prepare the data for the random forest
  data_randomforest <- data.frame(
    days = survival_train$days,
    vital_status = survival_train$vital_status,
    genes_expression_train)

  # Prepare the data for the random forest
  data_randomforest <- data_randomforest[, -which(names(data_randomforest) %in% c("patient"))]

  # Fit the Random Survival Forest model
  set.seed(1012)
  rsf_model <- rfsrc(Surv(days, vital_status) ~ .,
                     data = data_randomforest, seed = -123)

  # Extract feature importance from the causal forest
  importance_rsf <- vimp(rsf_model, seed = -123)$importance

  # Select the top features based on importance
  top_features_indices_rsf <- which(importance_rsf >= 0.001)

  # Convert the 'genes_expression_train' data frame to a matrix, excluding the column "patient" by finding its index and removing it
  z <- as.matrix(genes_expression_train[, -which(names(genes_expression_train) %in% c("patient"))])

  # Retrieve the column names (i.e., gene names) from the matrix 'z' and subset them using 'top_features_indices_rsf', which contains the indices of the top features (genes)
  top_features_rsf <- colnames(z)[top_features_indices_rsf]

  # remove NA features
  top_features_rsf <- top_features_rsf[!is.na(top_features_rsf)]

  # Fit a Cox regression model to avoid having predictors that are highly correlated with others
  set.seed(1012)
  fit_rsf_model <- cv.glmnet(as.matrix(subset(genes_expression_train, select = top_features_rsf)),
                             survival_object,
                             family = "cox",
                             alpha = best_alpha)

  # Extract coefficients from the fitted model
  RandomForest_Coefficients <- extractCoxRegressionCoefficients(fit_rsf_model)


  # ############# Fit the Causal forest model
  #
  # # Convert the 'genes_expression_train' data frame to a matrix, excluding the column "patient" by finding its index and removing it
  # X <- as.matrix(genes_expression_train[, -which(names(genes_expression_train) %in% c("patient"))])
  #
  # # Train the causal forest
  # set.seed(1012)
  # causal_forest_model <- causal_survival_forest(X,
  #                                               survival_train$days,
  #                                               rep(1, length(survival_train$days)),
  #                                               survival_train$vital_status,
  #                                               horizon = 1000, #max(survival_train$days, na.rm = TRUE))
  #                                               target = "survival.probability")
  #
  # # Extract feature importance from the causal forest
  # importance_causal_forest <- variable_importance(causal_forest_model)
  # cat("aqui: ", mean(importance_causal_forest))
  # # Select the top features based on importance
  # top_features_indices_causal_forest <- order(importance_causal_forest, decreasing = TRUE)[1:100]
  #
  # # Retrieve the column names (i.e., gene names) from the matrix 'x' and subset them using 'top_features_indices_causal_forest', which contains the indices of the top features (genes)
  # top_features_causal_forest <- colnames(X)[top_features_indices_causal_forest]
  #
  # # remove NA features
  # top_features_causal_forest <- top_features_causal_forest[!is.na(top_features_causal_forest)]
  #
  # # Fit a Cox regression model to avoid having predictors that are highly correlated with others
  # set.seed(1012)
  # fit_causal_forest_model <- cv.glmnet(as.matrix(subset(genes_expression_train, select = top_features_causal_forest)),
  #                                      survival_object,
  #                                      family = "cox",
  #                                      alpha = best_alpha)
  #
  #
  # # Extract coefficients from the fitted model
  # CausalForest_Coefficients <- extractCoxRegressionCoefficients(fit_causal_forest_model)
  #
  #
  ############# Fit the Boruta model
  set.seed(1012) # For reproducibility
  boruta_output <- Boruta(survival_object ~ .,
                          data = as.data.frame(prepareDataForCoxRegression(genes_expression_train)),
                          doTrace = 2)

# Get the significant features
significant_features <- getSelectedAttributes(boruta_output,
                                              withTentative = TRUE)


# Fit a Cox regression model to avoid having predictors that are highly correlated with others
set.seed(1012)
fit_boruta_model <- cv.glmnet(as.matrix(subset(genes_expression_train, select = significant_features)),
                              survival_object,
                              family = "cox",
                              alpha = best_alpha)

# Extract coefficients from the fitted model
Boruta_Coefficients <- extractCoxRegressionCoefficients(fit_boruta_model)



models_coefficients <- list(Cox_Coefficients, RandomForest_Coefficients, Boruta_Coefficients) #SIS_Coefficients, ISIS_Coefficients, CausalForest_Coefficients,

# Return the dataframe and the fitted models
return(models_coefficients)

}
