toposort_adj <- function(adj) {
  k <- nrow(adj)
  adj <- adj != 0
  indeg <- colSums(adj)

  order <- integer(0)
  queue <- which(indeg == 0)

  while (length(queue) > 0) {
    v <- queue[1]
    queue <- queue[-1]
    order <- c(order, v)

    children <- which(adj[v, ])
    for (c in children) {
      indeg[c] <- indeg[c] - 1
      if (indeg[c] == 0) queue <- c(queue, c)
    }
  }

  if (length(order) != k) {
    stop("Adjacency matrix contains a cycle; DAG required.")
  }
  order
}

construct_svar_template <- function(adj) {
  k <- nrow(adj)
  causal_order <- toposort_adj(adj)

  A_causal <- matrix(NA_real_, k, k)
  diag(A_causal) <- 1
  A_causal[upper.tri(A_causal)] <- 0

  inv_order <- order(causal_order)
  A_original <- A_causal[inv_order, inv_order]

  A_original
}

FEVD_TCI <- function(Phi,
                     Sigma,
                     nfore = 10) {

  if (!is.numeric(nfore) || length(nfore) != 1L || nfore <= 0) {
    stop("`nfore` must be a single positive integer.")
  }
  if (!is.matrix(Sigma) || nrow(Sigma) != ncol(Sigma)) {
    stop("`Sigma` must be a square covariance matrix.")
  }

  k <- nrow(Sigma)

  # ---- Parse Phi into lag matrices A1..Ap
  if (is.list(Phi)) {
    if (length(Phi) == 0L) stop("`Phi` cannot be an empty list.")
    ok <- vapply(Phi, function(M) is.matrix(M) && all(dim(M) == c(k, k)), logical(1))
    if (!all(ok)) stop("All elements of list `Phi` must be (k x k) matrices.")
    A_list <- Phi
  } else if (is.matrix(Phi)) {
    if (nrow(Phi) != k || ncol(Phi) %% k != 0) {
      stop("If `Phi` is a matrix, it must have dimension k x (k*p).")
    }
    p <- ncol(Phi) / k
    A_list <- lapply(seq_len(p), function(j) {
      Phi[, ((j - 1L) * k + 1L):(j * k), drop = FALSE]
    })
  } else {
    stop("`Phi` must be a list of matrices or a matrix.")
  }

  p <- length(A_list)

  # ---- MA coefficients Psi_h, h = 0..(nfore-1)
  Psi <- vector("list", nfore)
  Psi[[1L]] <- diag(k)

  if (nfore > 1L) {
    for (h in 2L:nfore) {
      Psi_h <- matrix(0, k, k)
      for (j in 1L:min(p, h - 1L)) {
        Psi_h <- Psi_h + A_list[[j]] %*% Psi[[h - j]]
      }
      Psi[[h]] <- Psi_h
    }
  }

  # ---- Generalized FEVD components
  denom_mat <- Reduce(`+`, lapply(Psi, function(P) P %*% Sigma %*% t(P)))
  denom <- diag(denom_mat)

  sig_diag <- diag(Sigma)
  if (any(!is.finite(sig_diag)) || any(sig_diag <= 0)) {
    stop("Diagonal of Sigma must be finite and strictly positive.")
  }

  enum <- matrix(0, k, k)
  for (h in seq_len(nfore)) {
    M <- Psi[[h]] %*% Sigma
    enum <- enum + (M ^ 2)
  }
  enum <- sweep(enum, 2, sig_diag, `/`)

  FEVD <- matrix(NA_real_, k, k)
  for (i in seq_len(k)) {
    if (is.finite(denom[i]) && denom[i] > 0) {
      FEVD[i, ] <- enum[i, ] / denom[i]
    }
  }

  # ---- Safe row normalization
  rs <- rowSums(FEVD, na.rm = TRUE)
  for (i in seq_len(k)) {
    if (!is.finite(rs[i]) || rs[i] <= 0) {
      FEVD[i, ] <- NA_real_
    } else {
      FEVD[i, ] <- FEVD[i, ] / rs[i]
    }
  }

  # ---- Connectedness Table → TCI
  CT <- FEVD * 100
  diag(CT) <- 0
  TO <- colSums(CT, na.rm = TRUE)
  mean(TO)
}

off_diag_percentage <- function(mat) {
  if (!is.matrix(mat)) stop("Input must be a matrix.")

  total_sum <- sum(mat)
  diag_sum <- sum(diag(mat))

  if (total_sum == 0) return(NA_real_)
  1 - (diag_sum / total_sum)
}

calculate_non_diag_percentage <- function(mat, direction = "row") {
  if (!is.matrix(mat)) stop("Input must be a matrix.")
  if (!(direction %in% c("row", "col"))) stop("Direction must be either 'row' or 'col'.")

  nrow_mat <- nrow(mat)
  ncol_mat <- ncol(mat)
  result <- numeric(if (direction == "row") nrow_mat else ncol_mat)

  if (direction == "row") {
    for (i in seq_len(nrow_mat)) {
      row_sum <- sum(mat[i, ])
      if (row_sum == 0) {
        result[i] <- NA_real_
        next
      }
      non_diag_sum <- sum(mat[i, -i])
      result[i] <- non_diag_sum / row_sum
    }
  } else {
    for (j in seq_len(ncol_mat)) {
      col_sum <- sum(mat[, j])
      if (col_sum == 0) {
        result[j] <- NA_real_
        next
      }
      non_diag_sum <- sum(mat[-j, j])
      result[j] <- non_diag_sum / col_sum
    }
  }

  result
}

extract_lag_matrices <- function(var_model) {
  K <- length(var_model$varresult)
  p <- var_model$p

  # Qualify coef to avoid R CMD check NOTE
  coef_vecs <- lapply(var_model$varresult, stats::coef)

  varnames <- colnames(var_model$y)
  A_list <- vector("list", length = p)

  for (j in seq_len(p)) {
    A_j <- matrix(NA_real_, nrow = K, ncol = K)
    coef_names <- paste0(varnames, ".", "l", j)

    for (eq in seq_len(K)) {
      A_j[, eq] <- coef_vecs[[eq]][coef_names]
    }
    A_list[[j]] <- t(A_j)
  }

  names(A_list) <- paste0("A", seq_len(p))
  A_list
}

calculate_irf <- function(var_model, n.ahead, ortho = FALSE, shock = "None") {
  K <- length(var_model$varresult)
  p <- var_model$p
  varnames <- colnames(var_model$y)

  A_list <- extract_lag_matrices(var_model)

  Sigma_u <- summary(var_model)$covres

  if (is.character(shock) && shock == "None") {
    shock_mat <- if (ortho) chol(Sigma_u) else diag(K)
  } else {
    shock_mat <- shock
  }

  IRF <- array(0, dim = c(K, K, n.ahead + 1))
  IRF[, , 1] <- shock_mat

  if (n.ahead > 0) {
    for (h in seq_len(n.ahead)) {
      IRF_h <- matrix(0, nrow = K, ncol = K)
      for (j in seq_len(min(h, p))) {
        IRF_h <- IRF_h + A_list[[j]] %*% IRF[, , h - j + 1]
      }
      IRF[, , h + 1] <- IRF_h
    }
  }

  dimnames(IRF) <- list(
    response = varnames,
    shock = varnames,
    horizon = 0:n.ahead
  )

  abs(IRF)
}

cumsum_irf <- function(IRF, h) {
  if (h > dim(IRF)[3]) {
    stop("h is larger than the number of horizons in IRF.")
  }

  n_row <- dim(IRF)[1]
  n_col <- dim(IRF)[2]
  A <- matrix(0, nrow = n_row, ncol = n_col)

  for (t in seq_len(h + 1)) {
    A <- A + IRF[, , t]
  }

  A
}
