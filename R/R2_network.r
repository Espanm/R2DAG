#' R2 network connectedness
#'
#' @param data data.frame or matrix with column names
#' @param directed logical
#' @param amat FALSE or adjacency matrix (child x parent)
#' @param dag_method "lingam", "notears", or "direct_lingam"
#' @param standardize_for_dag logical; if TRUE, variables are standardized
#'   columnwise before DAG learning. The rest of the function
#'   uses the original data scale.
#' @param mag integer seed used where relevant
#' @return A list with network tables and summary measures
#' @export
R2_network <- function(data,
                       directed = TRUE,
                       amat = FALSE,
                       dag_method = c("lingam", "notears", "direct_lingam"),
                       standardize_for_dag = FALSE,
                       mag = 123) {

  dag_method <- match.arg(dag_method)

  stopifnot(is.data.frame(data) || is.matrix(data))
  data <- as.data.frame(data)
  stopifnot(!is.null(colnames(data)))
  data <- stats::na.omit(data)

  p <- ncol(data)
  var_names <- colnames(data)
  X_raw <- as.matrix(data)

  # ============================================================
  # helpers
  # ============================================================

  .sym_sqrt <- function(M, tol = 1e-10) {
    ee <- eigen(M, symmetric = TRUE)
    vals <- pmax(ee$values, tol)
    ee$vectors %*% diag(sqrt(vals), length(vals)) %*% t(ee$vectors)
  }

  .sym_inv_sqrt <- function(M, tol = 1e-10) {
    ee <- eigen(M, symmetric = TRUE)
    vals <- pmax(ee$values, tol)
    ee$vectors %*% diag(1 / sqrt(vals), length(vals)) %*% t(ee$vectors)
  }

  .genizi_from_predictors <- function(y, X, tol = 1e-10) {
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

    P_half <- .sym_sqrt(Rxx, tol = tol)
    P_inv_half <- .sym_inv_sqrt(Rxx, tol = tol)

    z <- as.vector(P_inv_half %*% matrix(rxy, ncol = 1))
    G <- P_half %*% diag(z, nrow = length(z))
    gij <- rowSums(G^2)
    gij[abs(gij) < tol] <- 0

    out[keep_names] <- gij
    out
  }

  .get_ancestors <- function(A, node) {
    # A: child x parent
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
    # A: child x parent, so edges are parent -> child
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

  .path_product_B <- function(path, B) {
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
  # 0) Undirected default
  # ============================================================
  if (!directed && (is.logical(amat) && !amat)) {
    A <- matrix(1L, p, p)
    diag(A) <- 0L
    colnames(A) <- rownames(A) <- var_names
    amat <- A
  }

  # ============================================================
  # 1) Learn DAG if needed
  # ============================================================
  if ((is.logical(amat) && !amat) && directed) {

    X_dag <- if (standardize_for_dag) scale(X_raw) else X_raw

    if (dag_method == "lingam") {

      set.seed(mag)
      fit_lingam <- pcalg::lingam(X_dag, verbose = FALSE)

      Bm <- if (!is.null(fit_lingam$Bpruned)) {
        fit_lingam$Bpruned
      } else if (!is.null(fit_lingam$B)) {
        fit_lingam$B
      } else {
        stop("No 'B' or 'Bpruned' found in pcalg::lingam output.")
      }

      A <- matrix(0L, p, p)
      A[abs(Bm) > 1e-4] <- 1L
      diag(A) <- 0L

    } else if (dag_method == "direct_lingam") {

      if (!requireNamespace("reticulate", quietly = TRUE)) {
        stop("Package 'reticulate' is required for dag_method = 'direct_lingam'.")
      }

      lingam_available <- reticulate::py_module_available("lingam")

      if (!lingam_available) {
        message("Python package 'lingam' not found.")
        message("Installing now... this may take a few minutes.")
        reticulate::py_install("lingam")
        message("Installation finished. Continuing with DirectLiNGAM...")
      }

      lingam_py <- tryCatch(
        reticulate::import("lingam", delay_load = FALSE),
        error = function(e) {
          stop(
            "Could not import Python package 'lingam' even after installation attempt. ",
            "Original error: ", conditionMessage(e)
          )
        }
      )

      if (reticulate::py_module_available("numpy")) {
        np <- reticulate::import("numpy", delay_load = FALSE)
        np$random$seed(as.integer(mag))
      }

      model <- lingam_py$DirectLiNGAM()
      model$fit(X_dag)

      Bm <- reticulate::py_to_r(model$adjacency_matrix_)
      Bm <- as.matrix(Bm)

      if (!all(dim(Bm) == c(p, p))) {
        stop("DirectLiNGAM returned an adjacency matrix with unexpected dimensions.")
      }

      colnames(Bm) <- rownames(Bm) <- var_names

      A <- matrix(0L, p, p)
      A[abs(Bm) > 1e-4] <- 1L
      diag(A) <- 0L
      colnames(A) <- rownames(A) <- var_names

    } else if (dag_method == "notears") {

      W_hat <- gnlearn::notears(X_dag)

      W <- matrix(0, p, p)
      for (i in seq_along(W_hat)) {
        W[i, ] <- W_hat[i]
      }
      W <- t(W)

      A <- matrix(0L, p, p)
      A[abs(W) > 1e-4] <- 1L
      diag(A) <- 0L
      A <- t(A)  # child x parent
    }

    colnames(A) <- rownames(A) <- var_names
    amat <- A
  }

  if (is.logical(amat)) {
    stop("Provide an adjacency matrix 'amat' or set directed=TRUE.")
  }

  amat <- as.matrix(amat)
  colnames(amat) <- rownames(amat) <- var_names
  diag(amat) <- 0

  # ============================================================
  # 2) Correlation matrix
  # ============================================================
  corr_Y <- stats::cor(data, use = "pairwise.complete.obs")
  corr_Y <- (corr_Y + t(corr_Y)) / 2
  diag(corr_Y) <- 1

  # ============================================================
  # 3) Standardized structural betas from DAG
  # ============================================================
  if (directed) {
    Bstd <- compute_standardized_betas(corr_Y, amat)
    Bstd[amat == 0] <- 0
    diag(Bstd) <- 0
    colnames(Bstd) <- rownames(Bstd) <- var_names
  } else {
    Bstd <- NULL
  }

  # ============================================================
  # 4) TOTAL table:
  # target i is regressed only on its direct + indirect parents
  # ============================================================
  Total <- matrix(0, p, p, dimnames = list(var_names, var_names))

  if (directed) {
    for (i in seq_len(p)) {
      anc_idx <- .get_ancestors(amat, i)

      if (length(anc_idx) == 0) next

      y <- X_raw[, i]
      Xpred <- X_raw[, anc_idx, drop = FALSE]

      gij <- .genizi_from_predictors(y, Xpred)
      Total[i, names(gij)] <- gij
    }
  } else {
    for (i in seq_len(p)) {
      pred_idx <- setdiff(seq_len(p), i)
      y <- X_raw[, i]
      Xpred <- X_raw[, pred_idx, drop = FALSE]
      gij <- .genizi_from_predictors(y, Xpred)
      Total[i, pred_idx] <- gij
    }
  }

  diag(Total) <- 0

  # ============================================================
  # 5) DIRECT / INDIRECT allocation
  # ============================================================
  if (directed) {

    vars_raw <- apply(X_raw, 2, stats::var)
    names(vars_raw) <- var_names

    q_direct <- matrix(0, p, p, dimnames = list(var_names, var_names))
    q_indirect <- matrix(0, p, p, dimnames = list(var_names, var_names))

    for (j in seq_len(p)) {
      for (i in seq_len(p)) {
        if (i == j) next

        if (amat[i, j] != 0) {
          q_direct[i, j] <- (Bstd[i, j]^2) * (vars_raw[j] / vars_raw[i])
        }

        paths_ji <- .all_directed_paths(amat, from = j, to = i)

        if (length(paths_ji) > 0) {
          for (pth in paths_ji) {
            if (length(pth) >= 3) {
              prod_path <- .path_product_B(pth, Bstd)
              q_indirect[i, j] <- q_indirect[i, j] +
                (prod_path^2) * (vars_raw[j] / vars_raw[i])
            }
          }
        }
      }
    }

    Direct <- matrix(0, p, p, dimnames = list(var_names, var_names))
    Indirect <- matrix(0, p, p, dimnames = list(var_names, var_names))

    for (i in seq_len(p)) {
      for (j in seq_len(p)) {
        if (i == j) next

        gtot <- Total[i, j]
        qtot <- q_direct[i, j] + q_indirect[i, j]

        if (is.finite(gtot) && gtot > 0 && is.finite(qtot) && qtot > 0) {
          Direct[i, j]   <- gtot * q_direct[i, j] / qtot
          Indirect[i, j] <- gtot * q_indirect[i, j] / qtot
        }
      }
    }

    Direct[abs(Direct) < 1e-12] <- 0
    Indirect[abs(Indirect) < 1e-12] <- 0

    row_R2_direct   <- rowSums(Direct)
    row_R2_indirect <- rowSums(Indirect)
    row_R2_total    <- rowSums(Total)

    tci_direct   <- sum(row_R2_direct)   / (p - 1)
    tci_indirect <- sum(row_R2_indirect) / (p - 1)
    tci_total    <- sum(row_R2_total)    / (p - 1)

  } else {

    Direct <- Total
    Indirect <- NULL

    row_R2_total <- rowSums(Total)

    tci_direct   <- NULL
    tci_indirect <- NULL
    tci_total    <- sum(row_R2_total) / p
  }

  # ============================================================
  # 6) Network measures
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
    corr_Y         = corr_Y
  )
}