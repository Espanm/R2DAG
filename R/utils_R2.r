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

# ============================================================
# Matrix helpers
# ============================================================

sym_sqrt <- function(M, tol = 1e-10) {
  ee <- eigen(M, symmetric = TRUE)
  vals <- pmax(ee$values, tol)
  ee$vectors %*% diag(sqrt(vals), length(vals)) %*% t(ee$vectors)
}

sym_inv_sqrt <- function(M, tol = 1e-10) {
  ee <- eigen(M, symmetric = TRUE)
  vals <- pmax(ee$values, tol)
  ee$vectors %*% diag(1 / sqrt(vals), length(vals)) %*% t(ee$vectors)
}

# ============================================================
# Genizi decomposition from a predictor set
# ============================================================

genizi_from_predictors <- function(y, X, tol = 1e-10) {
  k <- ncol(X)
  out <- rep(0, k)
  names(out) <- colnames(X)

  if (k == 0) return(out)

  sdy <- stats::sd(y)
  sdx <- apply(X, 2, stats::sd)

  keep <- is.finite(sdx) & (sdx > tol)
  if (!is.finite(sdy) || sdy <= tol || !any(keep)) return(out)

  Xk <- X[, keep, drop = FALSE]
  keep_names <- colnames(Xk)

  Rxx <- stats::cor(Xk, use = "pairwise.complete.obs")
  rxy <- stats::cor(Xk, y, use = "pairwise.complete.obs")

  if (is.null(dim(Rxx))) {
    out[keep_names] <- as.numeric(rxy)^2
    return(out)
  }

  Rxx <- (Rxx + t(Rxx)) / 2
  diag(Rxx) <- 1

  P_half <- sym_sqrt(Rxx, tol = tol)
  P_inv_half <- sym_inv_sqrt(Rxx, tol = tol)

  z <- as.vector(P_inv_half %*% matrix(rxy, ncol = 1))
  G <- P_half %*% diag(z, nrow = length(z))
  gij <- rowSums(G^2)
  gij[abs(gij) < tol] <- 0

  out[keep_names] <- gij
  out
}

# ============================================================
# Graph helpers
# ============================================================

get_ancestors <- function(A, node) {
  # A is child x parent
  anc <- integer(0)
  frontier <- which(A[node, ] != 0)

  while (length(frontier) > 0) {
    new_nodes <- setdiff(frontier, anc)
    if (length(new_nodes) == 0) break

    anc <- union(anc, new_nodes)

    parents_of_new <- integer(0)
    for (v in new_nodes) {
      parents_of_new <- union(parents_of_new, which(A[v, ] != 0))
    }

    frontier <- setdiff(parents_of_new, anc)
  }

  sort(unique(anc))
}

all_directed_paths <- function(A, from, to) {
  # A is child x parent, so edges are parent -> child
  children_list <- lapply(seq_len(nrow(A)), function(j) which(A[, j] != 0))
  paths <- list()

  dfs <- function(current, target, visited, path) {
    if (current == target) {
      paths[[length(paths) + 1]] <<- path
      return(NULL)
    }

    nxt <- children_list[[current]]
    nxt <- setdiff(nxt, visited)

    if (length(nxt) == 0) return(NULL)

    for (v in nxt) {
      dfs(v, target, c(visited, v), c(path, v))
    }

    NULL
  }

  dfs(from, to, visited = from, path = from)
  paths
}

path_product_B <- function(path, B) {
  if (length(path) < 2) return(0)

  out <- 1
  for (h in seq_len(length(path) - 1)) {
    parent <- path[h]
    child  <- path[h + 1]
    out <- out * B[child, parent]
  }

  out
}

# ============================================================
# DAG estimation
# ============================================================

estimate_amat <- function(data,
                          dag_method = c("lingam", "notears"),
                          standardize = FALSE,
                          threshold = 1e-15) {

  dag_method <- match.arg(dag_method)

  stopifnot(is.data.frame(data) || is.matrix(data))
  data <- as.data.frame(data)
  stopifnot(!is.null(colnames(data)))

  X_raw <- as.matrix(data)
  X_dag <- if (standardize) scale(X_raw) else X_raw

  p <- ncol(X_raw)
  var_names <- colnames(X_raw)

  if (dag_method == "lingam") {
    set.seed(1119)
    fit_lingam <- pcalg::lingam(X_dag, verbose = FALSE)

    Bm <- if (!is.null(fit_lingam$Bpruned)) {
      fit_lingam$Bpruned
    } else if (!is.null(fit_lingam$B)) {
      fit_lingam$B
    } else {
      stop("No 'B' or 'Bpruned' found in pcalg::lingam output.")
    }

    A <- matrix(0L, p, p)
    A[abs(Bm) > threshold] <- 1L
    diag(A) <- 0L

  } else if (dag_method == "notears") {

    W_hat <- gnlearn::notears(X_dag)

    if (is.matrix(W_hat)) {
      W <- W_hat
    } else if (is.list(W_hat) && !is.null(W_hat$W)) {
      W <- as.matrix(W_hat$W)
    } else if (is.list(W_hat)) {
      W <- do.call(rbind, lapply(W_hat, as.numeric))
    } else {
      stop("Unsupported output format returned by gnlearn::notears().")
    }

    if (!all(dim(W) == c(p, p))) {
      stop("gnlearn::notears() returned an object with unexpected dimensions.")
    }

    A <- matrix(0L, p, p)
    A[abs(W) > threshold] <- 1L
    diag(A) <- 0L
    A <- t(A)  # convert to child x parent
  }

  colnames(A) <- rownames(A) <- var_names
  A
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

