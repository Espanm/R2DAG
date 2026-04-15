#' R2 network connectedness
#'
#' @param data data.frame or matrix with column names
#' @param directed logical
#' @param amat FALSE or adjacency matrix (child x parent)
#' @param dag_method "lingam" or "notears"
#' @param standardize_for_dag logical; if TRUE, variables are standardized
#'   columnwise before DAG learning. The rest of the function
#'   uses the original data scale.
#' @return A list with network tables, summary measures, and structural residuals
#' @export
R2_network <- function(data,
                       directed = TRUE,
                       amat = FALSE,
                       dag_method = c("lingam", "notears"),
                       standardize_for_dag = FALSE) {

  dag_method <- match.arg(dag_method)

  stopifnot(is.data.frame(data) || is.matrix(data))
  data <- as.data.frame(data)
  stopifnot(!is.null(colnames(data)))
  data <- stats::na.omit(data)

  p <- ncol(data)
  var_names <- colnames(data)
  X_raw <- as.matrix(data)
  n <- nrow(X_raw)

  # ============================================================
  # 0) Undirected default
  # ============================================================
  if (!directed && is.logical(amat) && !amat) {
    amat <- matrix(1L, p, p, dimnames = list(var_names, var_names))
    diag(amat) <- 0L
  }

  # ============================================================
  # 1) Learn DAG if needed
  # ============================================================
  if (directed && is.logical(amat) && !amat) {
    amat <- estimate_amat(
      data = X_raw,
      dag_method = dag_method,
      standardize = standardize_for_dag
    )
  }

  if (is.logical(amat)) {
    stop("Provide an adjacency matrix 'amat' or set directed = TRUE.")
  }

  amat <- as.matrix(amat)

  if (!all(dim(amat) == c(p, p))) {
    stop("'amat' must have dimensions p x p, where p = ncol(data).")
  }

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
  # 3b) Raw-scale structural coefficients and residuals
  # ============================================================
  if (directed) {
    Braw <- matrix(0, p, p, dimnames = list(var_names, var_names))
    intercepts <- stats::setNames(numeric(p), var_names)

    residual_mat <- matrix(
      NA_real_,
      nrow = n,
      ncol = p,
      dimnames = list(rownames(data), var_names)
    )

    for (i in seq_len(p)) {
      parent_idx <- which(amat[i, ] != 0)
      y <- X_raw[, i]

      if (length(parent_idx) == 0) {
        intercepts[i] <- mean(y)
        residual_mat[, i] <- y - intercepts[i]
      } else {
        fit_df <- data.frame(y = y, X_raw[, parent_idx, drop = FALSE])
        fit <- stats::lm(y ~ ., data = fit_df)

        cf <- stats::coef(fit)
        intercepts[i] <- unname(cf[1])
        Braw[i, parent_idx] <- unname(cf[-1])

        residual_mat[, i] <- stats::residuals(fit)
      }
    }

    structural_residuals <- as.data.frame(residual_mat)

    if (dag_method == "lingam") {
      lingam_residuals <- structural_residuals
    } else {
      lingam_residuals <- NULL
    }

  } else {
    Braw <- NULL
    intercepts <- NULL
    structural_residuals <- NULL
    lingam_residuals <- NULL
  }

  # ============================================================
  # 4) TOTAL table
  # ============================================================
  Total <- matrix(0, p, p, dimnames = list(var_names, var_names))

  if (directed) {
    for (i in seq_len(p)) {
      anc_idx <- get_ancestors(amat, i)

      if (length(anc_idx) == 0) next

      y <- X_raw[, i]
      Xpred <- X_raw[, anc_idx, drop = FALSE]

      gij <- genizi_from_predictors(y, Xpred)
      Total[i, names(gij)] <- gij
    }
  } else {
    for (i in seq_len(p)) {
      pred_idx <- setdiff(seq_len(p), i)
      y <- X_raw[, i]
      Xpred <- X_raw[, pred_idx, drop = FALSE]

      gij <- genizi_from_predictors(y, Xpred)
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

        paths_ji <- all_directed_paths(amat, from = j, to = i)

        if (length(paths_ji) > 0) {
          for (pth in paths_ji) {
            if (length(pth) >= 3) {
              prod_path <- path_product_B(pth, Bstd)
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
    direct_table         = Direct,
    indirect_table       = Indirect,
    total_table          = Total,
    to_total             = to_total,
    from_total           = from_total,
    net_total            = net_total,
    tci_direct           = tci_direct,
    tci_indirect         = tci_indirect,
    tci_total            = tci_total,
    npdc_total           = npdc_total,
    amat                 = amat,
    Bstd                 = Bstd,
    Braw                 = Braw,
    intercepts           = intercepts,
    corr_Y               = corr_Y,
    structural_residuals = structural_residuals,
    lingam_residuals     = lingam_residuals
  )
}