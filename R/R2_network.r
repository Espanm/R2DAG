#' R2 network connectedness
#'
#' @param data data.frame or matrix with column names
#' @param directed logical
#' @param amat FALSE or adjacency matrix (child x parent)
#' @param dag_method "lingam" or "notears"
#' @param standardize_for_dag logical; if TRUE, variables are standardized
#'   columnwise before DAG learning with LiNGAM. The rest of the function
#'   is still computed on the original data.
#' @return A list with network tables and summary measures
#' @export
R2_network <- function(data,
                       directed = TRUE,
                       amat = FALSE,
                       dag_method = c("lingam","notears"),
                       standardize_for_dag = FALSE) {

  dag_method <- match.arg(dag_method)

  stopifnot(is.data.frame(data) || is.matrix(data))
  data <- as.data.frame(data)
  stopifnot(!is.null(colnames(data)))

  # keep complete cases only
  data <- stats::na.omit(data)

  p <- ncol(data)
  var_names <- colnames(data)

  # ============================================================
  # helpers
  # ============================================================

  # symmetric matrix square root and inverse square root
  .sym_sqrt <- function(M, tol = 1e-10) {
    ee <- eigen(M, symmetric = TRUE)
    vals <- pmax(ee$values, tol)
    ee$vectors %*% diag(sqrt(vals), nrow = length(vals)) %*% t(ee$vectors)
  }

  .sym_inv_sqrt <- function(M, tol = 1e-10) {
    ee <- eigen(M, symmetric = TRUE)
    vals <- pmax(ee$values, tol)
    ee$vectors %*% diag(1 / sqrt(vals), nrow = length(vals)) %*% t(ee$vectors)
  }

  # Genizi decomposition for one target against a given predictor set
  .genizi_row_all_predictors <- function(target_index, X_df, tol = 1e-10) {
    y_name <- colnames(X_df)[target_index]
    x_names <- setdiff(colnames(X_df), y_name)

    out <- setNames(rep(0, ncol(X_df)), colnames(X_df))
    if (length(x_names) == 0) return(out)

    Xp <- as.matrix(X_df[, x_names, drop = FALSE])
    y  <- X_df[[y_name]]

    # standard deviations
    sdy <- stats::sd(y)
    sdx <- apply(Xp, 2, stats::sd)

    keep <- is.finite(sdx) & (sdx > tol)
    x_names_keep <- x_names[keep]

    if (length(x_names_keep) == 0 || !is.finite(sdy) || sdy <= tol) {
      return(out)
    }

    Xp <- as.matrix(X_df[, x_names_keep, drop = FALSE])

    # correlation objects
    Rxx <- stats::cor(Xp, use = "pairwise.complete.obs")
    rxy <- stats::cor(Xp, y, use = "pairwise.complete.obs")

    if (is.null(dim(Rxx))) {
      # one-predictor case
      gij <- as.numeric(rxy)^2
      out[x_names_keep] <- gij
      return(out)
    }

    # numerical stabilization
    Rxx <- (Rxx + t(Rxx)) / 2
    diag(Rxx) <- 1

    P_half     <- .sym_sqrt(Rxx, tol = tol)
    P_inv_half <- .sym_inv_sqrt(Rxx, tol = tol)

    z <- as.vector(P_inv_half %*% matrix(rxy, ncol = 1))
    G <- P_half %*% diag(z, nrow = length(z))
    gij <- rowSums(G^2)

    # tiny numerical negatives / drift
    gij[abs(gij) < tol] <- 0

    out[x_names_keep] <- gij
    out
  }

  # matrix of all directed paths of length >= 2 via powers of squared weights
  .indirect_from_squared_paths <- function(Bstd, tol = 1e-12) {
    S <- Bstd^2
    S[abs(S) < tol] <- 0
    p <- nrow(S)

    out <- matrix(0, p, p, dimnames = dimnames(S))
    if (p <= 2) return(out)

    current_power <- S
    for (k in 2:(p - 1)) {
      current_power <- current_power %*% S
      out <- out + current_power
    }
    out
  }

  # ============================================================
  # 0) If undirected and amat not provided: complete graph
  # ============================================================
  if (!directed && (is.logical(amat) && !amat)) {
    A <- matrix(1L, nrow = p, ncol = p)
    diag(A) <- 0L
    colnames(A) <- rownames(A) <- var_names
    amat <- A
  }

  # ============================================================
  # 1) Learn DAG if amat not provided
  # ============================================================
  if ((is.logical(amat) && !amat) && directed) {

    X_raw <- as.matrix(data)

    if (dag_method == "lingam") {

      X_dag <- if (standardize_for_dag) scale(X_raw) else X_raw

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

    } else if (dag_method == "notears") {

      X_dag <- X_raw
      W_hat <- gnlearn::notears(X_dag)

      W <- matrix(0, p, p)
      for (i in seq_along(W_hat)) {
        W[i, ] <- W_hat[i]
      }

      W <- t(W)

      A <- matrix(0L, p, p)
      A[abs(W) > 1e-4] <- 1L
      diag(A) <- 0L

      # convert to child x parent convention
      A <- t(A)
    }

    colnames(A) <- rownames(A) <- var_names
    amat <- A
  }

  # ============================================================
  # 1b) adjacency must exist
  # ============================================================
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
  # 3) TOTAL table = Genizi R² decomposition using all other vars
  # ============================================================
  Total <- matrix(0, p, p, dimnames = list(var_names, var_names))

  for (i in seq_len(p)) {
    Total[i, ] <- .genizi_row_all_predictors(i, data)
  }

  diag(Total) <- 0

  # ============================================================
  # 4) DIRECT / INDIRECT allocation
  # ============================================================
  if (directed) {

    # standardized structural coefficients implied by amat
    # convention: row = child, col = parent
    Bstd <- compute_standardized_betas(corr_Y, amat)
    Bstd[amat == 0] <- 0
    diag(Bstd) <- 0
    colnames(Bstd) <- rownames(Bstd) <- var_names

    # structural channel intensities
    q_direct <- Bstd^2
    q_indirect <- .indirect_from_squared_paths(Bstd)

    colnames(q_direct) <- rownames(q_direct) <- var_names
    colnames(q_indirect) <- rownames(q_indirect) <- var_names

    Direct   <- matrix(0, p, p, dimnames = list(var_names, var_names))
    Indirect <- matrix(0, p, p, dimnames = list(var_names, var_names))

    for (i in seq_len(p)) {
      for (j in seq_len(p)) {
        if (i == j) next

        g_total <- Total[i, j]
        q_total <- q_direct[i, j] + q_indirect[i, j]

        if (!is.finite(g_total) || g_total <= 0 || !is.finite(q_total) || q_total <= 0) {
          Direct[i, j]   <- 0
          Indirect[i, j] <- 0
        } else {
          direct_share   <- q_direct[i, j]   / q_total
          indirect_share <- q_indirect[i, j] / q_total

          Direct[i, j]   <- g_total * direct_share
          Indirect[i, j] <- g_total * indirect_share
        }
      }
    }

    # consistency-cleaning
    Direct[abs(Direct) < 1e-12] <- 0
    Indirect[abs(Indirect) < 1e-12] <- 0

    row_R2_direct   <- rowSums(Direct)
    row_R2_indirect <- rowSums(Indirect)
    row_R2_total    <- rowSums(Total)

    tci_direct   <- sum(row_R2_direct)   / (p - 1)
    tci_indirect <- sum(row_R2_indirect) / (p - 1)
    tci_total    <- sum(row_R2_total)    / (p - 1)

  } else {

    Bstd <- NULL
    q_direct <- NULL
    q_indirect <- NULL

    Direct <- Total
    Indirect <- NULL

    row_R2_direct <- NULL
    row_R2_total  <- rowSums(Total)

    tci_direct   <- NULL
    tci_indirect <- NULL
    tci_total    <- sum(row_R2_total) / p
  }

  # ============================================================
  # 5) Network measures
  # Convention kept the same as before:
  # rows = TO, cols = FROM  (a_ij = j -> i)
  # ============================================================
  to_total   <- colSums(Total)
  from_total <- rowSums(Total)
  net_total  <- to_total - from_total
  npdc_total <- Total - t(Total)

  list(
    direct_table        = Direct,
    indirect_table      = Indirect,
    total_table         = Total,
    to_total            = to_total,
    from_total          = from_total,
    net_total           = net_total,
    tci_direct          = tci_direct,
    tci_indirect        = tci_indirect,
    tci_total           = tci_total,
    npdc_total          = npdc_total,
    amat                = amat,
    Bstd                = Bstd,
    corr_Y              = corr_Y
  )
}