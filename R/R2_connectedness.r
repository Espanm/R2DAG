#library(pcalg)
#library(relaimpo)
#source("lingam.r")

R2_network <- function(data, method, directed=TRUE, amat=FALSE) {

  if (amat==FALSE || directed==TRUE){
    lingam_amat <- data2amat(data)
    amat <- lingam_amat$amat
  }

  # Assign column names to adjacency matrix
  colnames(amat) <- rownames(amat) <- colnames(data)

  # Get the number of variables
  p <- ncol(data)

  # Initialize result matrix
  result_matrix <- matrix(0, nrow = p, ncol = p)
  colnames(result_matrix) <- rownames(result_matrix) <- colnames(data)

  # Initialize tci vector
  tci <- c()

  if (directed==TRUE){
    # Compute reachability paths
    paths <- compute_reachability(amat)
  }

  # Iterate through nodes
  for (i in seq_len(p)) {

    if (directed==TRUE){
    # Identify parent nodes (nodes with paths leading to i)
      selected_cols <- paths[[colnames(data)[i]]]
    } else {
      selected_cols <-  colnames(data[, -i, drop=FALSE])
    }


    # Ensure at least one predictor exists
    if (length(selected_cols) > 0) {
      selected_cols <- c(colnames(data)[i], selected_cols)  # Add target variable

      # Subset data for regression
      filtered_data <- data[, selected_cols, drop = FALSE]

      # Fit linear model: target variable as dependent, others as independent
      model <- lm(filtered_data[, 1] ~ ., data = filtered_data[, -1, drop = FALSE])

      # Compute R-squared value
      r2 <- summary(model)$r.squared

      tci <- c(tci, r2)

      # If only one predictor, assign R² directly
      if (length(selected_cols) == 2) {
        result_matrix[i, selected_cols[2]] <- r2
      } else {
        # Compute relative importance only if multiple predictors exist
        rel_importance <- calc.relimp(model, type = method, rela = TRUE)

        # Extract slot values dynamically
        slot_values <- slot(rel_importance, tolower(method)) * r2

        # Store values in result matrix
        valid_cols <- intersect(names(slot_values), colnames(result_matrix))
        result_matrix[i, valid_cols] <- slot_values[valid_cols]
      }
    }
  }

  result <- list()
  result$table <- result_matrix
  result$from <- colSums(result_matrix)
  result$to <- rowSums(result_matrix)
  result$tci <- mean(tci)
  result$lingam <- lingam_amat$lingam

  return(result)
}





