#' R2 network connectedness
#'
#' @param data data.frame or matrix with column names
#' @param directed logical
#' @param amat FALSE or adjacency matrix (child x parent)
#' @param dag_method "lingam" or "notears"
#' @return A list with network tables and summary measures
#' @export
R2_network <- function(data,
                       directed = TRUE,
                       amat = FALSE,
                       dag_method = c("lingam","notears")) {

  dag_method <- match.arg(dag_method)

  stopifnot(is.data.frame(data) || is.matrix(data))
  data <- as.data.frame(data)
  stopifnot(!is.null(colnames(data)))
  p <- ncol(data)

  # 0) If undirected and amat not provided: complete graph
  if (!directed && (is.logical(amat) && !amat)) {
    A <- matrix(1L, nrow = p, ncol = p)
    diag(A) <- 0L
    colnames(A) <- rownames(A) <- colnames(data)
    amat <- A
  }

  # ------------------------------------------------------------
  # 1) Learn DAG if amat not provided
  # ------------------------------------------------------------
  if ((is.logical(amat) && !amat) && directed) {

    X <- as.matrix(data)

    if (dag_method == "lingam") {

      fit_lingam <- pcalg::lingam(X, verbose = FALSE)

      Bm <- if (!is.null(fit_lingam$Bpruned)) fit_lingam$Bpruned
      else if (!is.null(fit_lingam$B)) fit_lingam$B
      else stop("No 'B' or 'Bpruned' found in pcalg::lingam output.")

      A <- matrix(0L, p, p)
      A[abs(Bm) > 1e-4] <- 1L
      diag(A) <- 0L

    } else if (dag_method == "notears") {

      # run NOTEARS
      W_hat <- gnlearn::notears(X)

      # convert list -> matrix
      W <- matrix(0, p, p)

      for(i in seq_along(W_hat)){
        W[i, ] <- W_hat[i]
      }

      W <- t(W)

      # adjacency
      A <- matrix(0L, p, p)
      A[abs(W) > 1e-4] <- 1L
      diag(A) <- 0L

      # convert to child x parent convention
      A <- t(A)
    }

    colnames(A) <- rownames(A) <- colnames(data)
    amat <- A
  }

  # ------------------------------------------------------------
  # 1b) adjacency must exist
  # ------------------------------------------------------------
  if (is.logical(amat)) {
    stop("Provide an adjacency matrix 'amat' or set directed=TRUE.")
  }

  amat <- as.matrix(amat)
  colnames(amat) <- rownames(amat) <- colnames(data)

  # ------------------------------------------------------------
  # 2) Correlation matrix
  # ------------------------------------------------------------
  corr_Y <- stats::cor(data, use = "pairwise.complete.obs")

  # ------------------------------------------------------------
  # 3) Direct contributions
  # ------------------------------------------------------------
  Direct <- matrix(0, p, p, dimnames = list(colnames(data), colnames(data)))

  for (i in seq_len(p)) {
    Direct[i, ] <- .direct_row_genizi_fast(i, corr_Y, amat)
  }

  # ------------------------------------------------------------
  # 4) Indirect (only directed)
  # ------------------------------------------------------------
  if (directed) {

    Bstd <- compute_standardized_betas(corr_Y, amat)

    Indirect <- .indirect_from_powers_abs(Bstd) / 2
    colnames(Indirect) <- rownames(Indirect) <- colnames(data)

    Direct <- Direct / 2

    Total <- Direct + Indirect

    row_R2_direct   <- rowSums(Direct)
    row_R2_indirect <- rowSums(Indirect)
    row_R2_total    <- rowSums(Total)

    tci_direct   <- sum(row_R2_direct)   / (p - 1)
    tci_indirect <- sum(row_R2_indirect) / (p - 1)
    tci_total    <- sum(row_R2_total)    / (p - 1)

  } else {

    Bstd <- NULL
    Indirect <- NULL
    Total <- Direct

    row_R2_direct <- NULL
    row_R2_total  <- rowSums(Total)

    tci_direct   <- NULL
    tci_indirect <- NULL
    tci_total    <- sum(row_R2_total) / p
  }

  # ------------------------------------------------------------
  # 5) Network measures
  # ------------------------------------------------------------
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