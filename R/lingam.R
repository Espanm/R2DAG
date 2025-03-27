#library(pcalg)

data2amat <- function(data){

  # Learn the Bayesian network structure using LiNGAM
  lingam_model <- lingam(data, verbose=FALSE)
  # Extract adjacency matrix
  coefs <- as.matrix(lingam_model$Bpruned)
  # Convert to binary adjacency matrix
  amat <- ifelse(coefs != 0, 1, 0)

  res <- c()
  res$lingam <- lingam_model
  res$amat <- amat

  return(res)
}

compute_reachability <- function(amat) {
  # Ensure adjacency matrix is logical (TRUE/FALSE)
  reachability <- amat > 0

  n <- nrow(amat)
  node_names <- colnames(amat)

  # Compute transitive closure
  for (k in seq_len(n)) {
    reachability <- reachability | (reachability %*% reachability > 0)
  }

  # Initialize inverse reachability dictionary
  reach_dict <- setNames(vector("list", n), node_names)

  for (j in seq_len(n)) {
    # For node j, find all i such that i can reach j
    reachable_from <- which(reachability[j, ])
    reach_dict[[node_names[j]]] <- node_names[reachable_from]
  }

  return(reach_dict)
}

count_isolated_nodes <- function(adj_mat) {
  # Check if it's a square matrix
  if (!is.matrix(adj_mat) || nrow(adj_mat) != ncol(adj_mat)) {
    stop("Input must be a square adjacency matrix.")
  }

  # An isolated node has no incoming and no outgoing edges
  row_sums <- rowSums(adj_mat)
  col_sums <- colSums(adj_mat)

  isolated <- (row_sums + col_sums) == 0

  return(sum(isolated))
}
