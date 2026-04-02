#' Theoretical values of the R2/DY network under multiple methodologies
#'
#' @param A0 Structural matrix (p x p)
#' @param sigma_eps Structural shock standard deviations
#' @param var_names Variable names
#' @param method One of:
#'   "directed"          = DAG-based total/direct/indirect methodology
#'   "undirected_genizi" = undirected Genizi decomposition
#'   "dy"                = theoretical Diebold-Yilmaz network under VAR(0)
#' @return A list with theoretical network tables and measures
#' @export
theoretical_R2_network <- function(A0,
                                   sigma_eps = rep(1, nrow(A0)),
                                   var_names = paste0("X", 1:nrow(A0)),
                                   method = c("directed",
                                              "undirected_genizi",
                                              "dy")) {

  method <- match.arg(method)
  p <- nrow(A0)

  stopifnot(
    is.matrix(A0),
    ncol(A0) == p,
    length(sigma_eps) == p
  )

  # ============================================================
  # helpers
  # ============================================================

  .sym_sqrt <- function(M, tol = 1e-10) {
    ee <- eigen(M, symmetric = TRUE)
    vals <- pmax(ee$values, tol)
    ee$vectors %*% diag(sqrt(vals)) %*% t(ee$vectors)
  }

  .sym_inv_sqrt <- function(M, tol = 1e-10) {
    ee <- eigen(M, symmetric = TRUE)
    vals <- pmax(ee$values, tol)
    ee$vectors %*% diag(1 / sqrt(vals)) %*% t(ee$vectors)
  }

  .genizi_from_corr <- function(Rxx, rxy, pred_names = NULL, tol = 1e-10) {

    k <- length(rxy)

    if (is.null(pred_names)) {
      pred_names <- paste0("X", seq_len(k))
    }

    out <- rep(0, k)
    names(out) <- pred_names

    if (k == 0) return(out)

    if (k == 1) {
      out[1] <- as.numeric(rxy)^2
      return(out)
    }

    Rxx <- (Rxx + t(Rxx)) / 2
    diag(Rxx) <- 1

    P_half <- .sym_sqrt(Rxx, tol)
    P_inv_half <- .sym_inv_sqrt(Rxx, tol)

    z <- as.vector(P_inv_half %*% matrix(rxy, ncol = 1))
    G <- P_half %*% diag(z)
    gij <- rowSums(G^2)

    gij[abs(gij) < tol] <- 0
    out[] <- gij
    out
  }

  .get_ancestors <- function(A, node) {
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

  .all_directed_paths <- function(A, from, to) {
    children_list <- lapply(seq_len(nrow(A)), function(j) which(A[, j] != 0))
    paths <- list()

    dfs <- function(current, target, visited, path) {
      if (current == target) {
        paths[[length(paths) + 1]] <<- path
        return()
      }
      nxt <- setdiff(children_list[[current]], visited)
      for (v in nxt) {
        dfs(v, target, c(visited, v), c(path, v))
      }
    }

    dfs(from, to, from, from)
    paths
  }

  .path_product_B <- function(path, B) {
    out <- 1
    for (h in seq_len(length(path) - 1)) {
      out <- out * B[path[h + 1], path[h]]
    }
    out
  }

  # ============================================================
  # covariance + correlation
  # ============================================================

  A0_inv <- solve(A0)
  Sigma_eps <- diag(sigma_eps^2)

  cov_Y  <- A0_inv %*% Sigma_eps %*% t(A0_inv)
  corr_Y <- cov2cor(cov_Y)

  colnames(corr_Y) <- rownames(corr_Y) <- var_names
  colnames(cov_Y)  <- rownames(cov_Y)  <- var_names

  vars_raw <- diag(cov_Y)

  # ============================================================
  # adjacency
  # ============================================================

  if (method == "directed") {
    amat <- (abs(A0) > 0) * 1
    diag(amat) <- 0
  } else {
    amat <- matrix(1, p, p)
    diag(amat) <- 0
  }

  colnames(amat) <- rownames(amat) <- var_names

  # ============================================================
  # METHOD 1: DIRECTED (new methodology)
  # ============================================================

  if (method == "directed") {

    Bstd <- compute_standardized_betas(corr_Y, amat)
    Bstd[amat == 0] <- 0
    diag(Bstd) <- 0

    Total <- matrix(0, p, p, dimnames = list(var_names, var_names))

    for (i in seq_len(p)) {
      anc_idx <- .get_ancestors(amat, i)
      if (length(anc_idx) == 0) next

      Rxx <- corr_Y[anc_idx, anc_idx, drop = FALSE]
      rxy <- corr_Y[anc_idx, i]

      gij <- .genizi_from_corr(Rxx, rxy, var_names[anc_idx])
      Total[i, names(gij)] <- gij
    }

    diag(Total) <- 0

    q_direct   <- matrix(0, p, p, dimnames = list(var_names, var_names))
    q_indirect <- matrix(0, p, p, dimnames = list(var_names, var_names))

    for (j in seq_len(p)) {
      for (i in seq_len(p)) {
        if (i == j) next

        if (amat[i, j] != 0) {
          q_direct[i, j] <- (Bstd[i, j]^2) * (vars_raw[j] / vars_raw[i])
        }

        paths <- .all_directed_paths(amat, j, i)

        for (pth in paths) {
          if (length(pth) >= 3) {
            prod <- .path_product_B(pth, Bstd)
            q_indirect[i, j] <- q_indirect[i, j] +
              (prod^2) * (vars_raw[j] / vars_raw[i])
          }
        }
      }
    }

    Direct   <- matrix(0, p, p, dimnames = list(var_names, var_names))
    Indirect <- matrix(0, p, p, dimnames = list(var_names, var_names))

    for (i in seq_len(p)) {
      for (j in seq_len(p)) {

        gtot <- Total[i, j]
        qtot <- q_direct[i, j] + q_indirect[i, j]

        if (gtot > 0 && qtot > 0) {
          Direct[i, j]   <- gtot * q_direct[i, j] / qtot
          Indirect[i, j] <- gtot * q_indirect[i, j] / qtot
        }
      }
    }

    row_R2_total <- rowSums(Total)

    tci_direct   <- sum(rowSums(Direct))   / (p - 1)
    tci_indirect <- sum(rowSums(Indirect)) / (p - 1)
    tci_total    <- sum(row_R2_total)      / (p - 1)
  }

  # ============================================================
  # METHOD 2: UNDIRECTED GENIZI
  # ============================================================

  if (method == "undirected_genizi") {

    Bstd <- NULL
    Indirect <- NULL

    Total <- matrix(0, p, p, dimnames = list(var_names, var_names))

    for (i in seq_len(p)) {
      idx <- setdiff(seq_len(p), i)

      Rxx <- corr_Y[idx, idx, drop = FALSE]
      rxy <- corr_Y[idx, i]

      gij <- .genizi_from_corr(Rxx, rxy, var_names[idx])
      Total[i, names(gij)] <- gij
    }

    diag(Total) <- 0

    Direct <- Total

    row_R2_total <- rowSums(Total)

    tci_direct   <- NULL
    tci_indirect <- NULL
    tci_total    <- sum(row_R2_total) / p
  }

  # ============================================================
  # METHOD 3: DY under VAR(0)
  # generalized FEVD, row-normalized
  # ============================================================

  if (method == "dy") {

    Bstd <- NULL
    Indirect <- NULL

    sigma_diag <- diag(cov_Y)

    theta_raw <- matrix(0, p, p, dimnames = list(var_names, var_names))

    for (i in seq_len(p)) {
      for (j in seq_len(p)) {
        theta_raw[i, j] <- (cov_Y[i, j]^2) / (sigma_diag[i] * sigma_diag[j])
      }
    }

    theta_norm <- theta_raw

    for (i in seq_len(p)) {
      s <- sum(theta_raw[i, ])
      if (s > 0) {
        theta_norm[i, ] <- theta_raw[i, ] / s
      } else {
        theta_norm[i, ] <- 0
      }
    }

    Total <- theta_norm

    Direct <- Total

    row_R2_total <- rowSums(Total)

    tci_direct   <- NULL
    tci_indirect <- NULL
    tci_total    <- sum(row_R2_total) / p
  }

  # ============================================================
  # common outputs
  # ============================================================

  to_total   <- colSums(Total)
  from_total <- rowSums(Total)
  net_total  <- to_total - from_total
  npdc_total <- Total - t(Total)

  list(
    direct_table   = Direct,
    indirect_table = Indirect,
    total_table    = Total,
    to_total       = to_total,
    from_total     = from_total,
    net_total      = net_total,
    tci_direct     = tci_direct,
    tci_indirect   = tci_indirect,
    tci_total      = tci_total,
    npdc_total     = npdc_total,
    amat           = amat,
    Bstd           = Bstd,
    corr_Y         = corr_Y,
    cov_Y          = cov_Y
  )
}