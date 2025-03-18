#library(pcalg)

data2amat <- function(data){

  # Learn the Bayesian network structure using LiNGAM
  lingam_model <- lingam(data, verbose=FALSE)
  # Extract adjacency matrix
  coefs <- t(as.matrix(lingam_model$Bpruned))
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

    # Compute transitive closure using matrix multiplication for efficiency
    for (k in seq_len(n)) {
      reachability <- reachability | (reachability %*% reachability > 0)
    }

    # Convert reachability matrix into a dictionary with column names
    reach_dict <- setNames(vector("list", n), node_names)

    for (i in seq_len(n)) {
      reachable_nodes <- which(reachability[, i])  # Get indices of reachable nodes
      reach_dict[[node_names[i]]] <- node_names[reachable_nodes]  # Map indices to column names
    }

    return(reach_dict)
  }
