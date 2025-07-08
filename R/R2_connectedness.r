#library(pcalg)
#library(relaimpo)
#source("lingam.r")


R2_network <- function(data, method="genizi", directed=TRUE, amat=FALSE) {

  if ((is.logical(amat) && !amat) && directed) {
    lingam_amat <- data2amat(data)
    amat <- lingam_amat$amat
  }

  if (!is.logical(amat)){
    # Assign column names to adjacency matrix
    colnames(amat) <- rownames(amat) <- colnames(data)
    s <- count_isolated_nodes(amat)
  }

  # Get the number of variables
  p <- ncol(data)

  # Initialize result matrix
  result_matrix <- matrix(0, nrow = p, ncol = p)
  colnames(result_matrix) <- rownames(result_matrix) <- colnames(data)

  # Initialize tci vector
  tci <- c(0)

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
  result$to <- colSums(result_matrix)
  result$from <- rowSums(result_matrix)
  result$net <- result$to - result$from

  valid_tci <- tci[is.finite(tci)]
  nonzero_from <- sum(result$from != 0)

  if (length(valid_tci) == 0 || nonzero_from == 0) {
    result$tci <- 0
  } else {
    if (directed){
      result$tci <- sum(valid_tci) / (p-1)
    }
    else{
      result$tci <- sum(valid_tci) / p
    }
  }

  result$amat <- amat
  result$npdc <- result$table - t(result$table)


  return(result)
}



rolling_network <- function(data, block_size, method="genizi", directed=TRUE, dynamic_structure=TRUE){

  if (dynamic_structure==FALSE){

  lingam_model <- lingam(data, verbose=FALSE)
  # Extract adjacency matrix
  coefs <- as.matrix(lingam_model$Bpruned)
  # Convert to binary adjacency matrix
  amat <- ifelse(coefs != 0, 1, 0)

  n <- nrow(data)
  p <- ncol(data)

  tci_vector <- c()
  to_df <- data.frame(matrix(ncol = p, nrow = 0))
  colnames(to_df) <- colnames(data)
  from_df <- data.frame(matrix(ncol = p, nrow = 0))
  colnames(from_df) <- colnames(data)

  for (i in 1:(n-block_size+1)){
    rolling_data <- data[i:(i+block_size-1),]
    network <- R2_network(rolling_data, method=method, directed=directed, amat=amat)
    tci_vector <- c(network$tci, tci_vector)
    to_df <- rbind(to_df, as.data.frame(as.list(network$to)))
    from_df <- rbind(from_df, as.data.frame(as.list(network$from)))

  }
  }
  else{
    n <- nrow(data)
    p <- ncol(data)

    tci_vector <- c()
    to_df <- data.frame(matrix(ncol = p, nrow = 0))
    colnames(to_df) <- colnames(data)
    from_df <- data.frame(matrix(ncol = p, nrow = 0))
    amat_list <- list()
    colnames(from_df) <- colnames(data)

    for (i in 1:(n-block_size)){
      rolling_data <- data[i:(i+block_size),]
      network <- R2_network(rolling_data, method=method, directed=directed, amat=FALSE)
      tci_vector <- c(network$tci, tci_vector)
      to_df <- rbind(to_df, as.data.frame(as.list(network$to)))
      from_df <- rbind(from_df, as.data.frame(as.list(network$from)))
      amat_list[[i]] <- network$amat

    }
  }

  result <- c()
  result$tci_vector <- tci_vector
  result$to_df <- to_df
  result$from_df <- from_df
  result$amat_list <- amat_list

  return(result)
}

theoretical_calc <- function(A0, colnames = paste0("X", 1:nrow(A0)), directed = T) {
  p <- nrow(A0)
  A0_inv <- solve(A0)
  cov_Y <- A0_inv %*% t(A0_inv)
  corr_Y <- cov2cor(cov_Y)

  result_matrix <- matrix(0, nrow = p, ncol = p)
  colnames(result_matrix) <- rownames(result_matrix) <- colnames

  amat <- matrix(0, nrow = p, ncol = p)
  for (i in 1:p) {
    for (j in 1:p) {
      if (i != j && A0[i, j] != 0) {
        amat[i, j] <- 1
      }
    }
  }
  colnames(amat) <- rownames(amat) <- colnames

  if (directed) {
    paths <- compute_reachability(amat)
  }

  tci_values <- numeric(p)

  for (i in 1:p) {
    varname <- colnames[i]

    if (directed) {
      parent_names <- paths[[varname]]
    } else {
      parent_names <- colnames[-i]  # all other variables
    }

    parents <- which(colnames %in% parent_names)

    if (length(parents) == 0) {
      tci_values[i] <- 0
      next
    }

    r_yx <- corr_Y[i, parents]
    R_xx <- corr_Y[parents, parents, drop = FALSE]

    if (length(parents) == 1) {
      contrib <- r_yx^2
      names(contrib) <- colnames[parents]
    } else {
      eig <- eigen(R_xx)
      if (any(eig$values <= 0)) {
        warning(paste("R_xx is not positive definite at node", i))
        next
      }

      R_inv_sqrt <- eig$vectors %*% diag(1 / sqrt(eig$values)) %*% t(eig$vectors)
      c_vec <- R_inv_sqrt %*% r_yx
      contrib <- c_vec^2
      names(contrib) <- colnames[parents]
    }

    result_matrix[i, names(contrib)] <- contrib
    tci_values[i] <- sum(contrib)
  }

  result <- list()
  result$table <- result_matrix
  result$to <- colSums(result_matrix)
  result$from <- rowSums(result_matrix)
  result$net <- result$to - result$from
  if (directed){
    result$tci <- sum(tci_values) / (p-1)
  }
  else{
    result$tci <- sum(tci_values) / p
  }
  result$amat <- amat
  result$npdc <- result_matrix - t(result_matrix)

  return(result)
}






