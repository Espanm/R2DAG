# =========================
# Speed helper functions
# =========================

# Fast Genizi "direct" contributions for one row i using only the correlation matrix.
.direct_row_genizi_fast <- function(i, corr_Y, amat) {
  p <- nrow(corr_Y)
  out <- stats::setNames(numeric(p), colnames(corr_Y))

  parents <- which(amat[i, ] == 1L)
  if (length(parents) == 0) return(out)

  r_yx <- corr_Y[i, parents, drop = TRUE]

  # Single parent: contribution is simply corr^2
  if (length(parents) == 1L) {
    out[parents] <- r_yx^2
    return(out)
  }

  # Multiple parents: Genizi via whitening of R_xx
  R_xx <- corr_Y[parents, parents, drop = FALSE]
  eig  <- eigen(R_xx, symmetric = TRUE)

  # Numerical floor on eigenvalues
  eig$values[eig$values < .Machine$double.eps] <- .Machine$double.eps

  R_inv_sqrt <- eig$vectors %*% diag(1 / sqrt(eig$values)) %*% t(eig$vectors)
  c_vec <- as.numeric(R_inv_sqrt %*% r_yx)

  out[parents] <- c_vec^2
  out
}

# Indirect contributions using absolute path products: |B^2| + |B^3| + ...
.indirect_from_powers_abs <- function(Bstd) {
  p <- nrow(Bstd)
  Ind <- matrix(0, p, p)
  if (p <= 1) return(Ind)

  Pk <- Bstd %*% Bstd
  Ind <- Ind + abs(Pk)

  for (k in 3:p) {
    Pk <- Pk %*% Bstd
    if (all(Pk == 0)) break
    Ind <- Ind + abs(Pk)
  }

  diag(Ind) <- 0

  Ind
}


# Compute standardized beta matrix (stable and fast).
compute_standardized_betas <- function(corr_Y, amat) {
  p <- nrow(amat)
  Bstd <- matrix(0, p, p)

  for (i in seq_len(p)) {
    parents <- which(amat[i, ] == 1L)
    if (length(parents) == 0) next

    r_yx <- corr_Y[i, parents, drop = TRUE]

    # Single predictor: standardized beta equals correlation
    if (length(parents) == 1L) {
      Bstd[i, parents] <- as.numeric(r_yx)
      next
    }

    R_xx <- corr_Y[parents, parents, drop = FALSE]
    beta_std <- tryCatch(solve(R_xx, r_yx), error = function(e) NA)

    if (all(is.finite(beta_std))) {
      Bstd[i, parents] <- beta_std
    } else {
      # Fallback: eigen-based inverse if solve() is unstable
      eig <- eigen(R_xx, symmetric = TRUE)
      eig$values[eig$values < .Machine$double.eps] <- .Machine$double.eps
      R_inv <- eig$vectors %*% diag(1 / eig$values) %*% t(eig$vectors)
      Bstd[i, parents] <- as.numeric(R_inv %*% r_yx)
    }
  }

  colnames(Bstd) <- rownames(Bstd) <- colnames(corr_Y)
  Bstd
}

# =========================
# LiNGAM -> adjacency matrix
# =========================
data2amat <- function(data) {
  # pcalg is in Imports, so it is available when the package loads
  lingam_model <- pcalg::lingam(as.matrix(data), verbose = FALSE)

  coefs <- as.matrix(lingam_model$Bpruned)
  amat <- ifelse(coefs != 0, 1L, 0L)
  diag(amat) <- 0L

  list(
    lingam = lingam_model,
    amat = amat
  )
}

# =========================
# Reachability (transitive closure)
# =========================
compute_reachability <- function(amat) {
  # amat: adjacency matrix where amat[i, j] = 1 means i -> j (row -> col)
  reachability <- amat > 0
  n <- nrow(amat)
  node_names <- colnames(amat)

  # Transitive closure (boolean matrix multiplication)
  for (k in seq_len(n)) {
    reachability <- reachability | (reachability %*% reachability > 0)
  }

  # For each target node j: list nodes i such that i can reach j
  reach_dict <- stats::setNames(vector("list", n), node_names)

  for (j in seq_len(n)) {
    reachable_from <- which(reachability[, j])
    reach_dict[[node_names[j]]] <- node_names[reachable_from]
  }

  reach_dict
}

# =========================
# Simple graph diagnostic
# =========================
count_isolated_nodes <- function(adj_mat) {
  if (!is.matrix(adj_mat) || nrow(adj_mat) != ncol(adj_mat)) {
    stop("Input must be a square adjacency matrix.")
  }

  row_sums <- rowSums(adj_mat)
  col_sums <- colSums(adj_mat)

  isolated <- (row_sums + col_sums) == 0
  sum(isolated)
}

# =========================
# Bootstrap LiNGAM structures
# =========================
bootstrap_lingam <- function(data, B = 100) {
  n <- nrow(data)
  p <- ncol(data)

  amat_keys <- vector("character", B)

  for (b in seq_len(B)) {
    idx <- sample.int(n, size = n, replace = TRUE)
    resampled_data <- data[idx, , drop = FALSE]

    lingam_fit <- pcalg::lingam(as.matrix(resampled_data), verbose = FALSE)
    amat <- ifelse(as.matrix(lingam_fit$Bpruned) != 0, 1L, 0L)
    diag(amat) <- 0L

    # Store as a string key for frequency counting
    amat_keys[b] <- paste(amat, collapse = "")
  }

  tab <- sort(table(amat_keys), decreasing = TRUE)
  keys <- names(tab)

  unique_amats <- lapply(keys, function(s) {
    matrix(as.integer(strsplit(s, "", fixed = TRUE)[[1]]),
           nrow = p, byrow = FALSE)
  })

  list(
    matrices = unique_amats,
    frequencies = as.numeric(tab) / B
  )
}
