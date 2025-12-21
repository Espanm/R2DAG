#' R2 network connectedness
#'
#' @param data data.frame or matrix with column names
#' @param method character, default "genizi"
#' @param directed logical
#' @param amat FALSE or adjacency matrix (child x parent)
#' @param prev_amat optional (currently unused)
#' @return A list with network tables and summary measures
#' @export
R2_network <- function(data,
                       method = "genizi",
                       directed = TRUE,
                       amat = FALSE,
                       prev_amat = NULL) {

  stopifnot(is.data.frame(data) || is.matrix(data))
  data <- as.data.frame(data)
  stopifnot(!is.null(colnames(data)))
  p <- ncol(data)

  # 1) LiNGAM -> amat[child, parent] = 1
  if ((is.logical(amat) && !amat) && directed) {
    X <- as.matrix(data)
    fit_lingam <- pcalg::lingam(X, verbose = FALSE)

    Bm <- if (!is.null(fit_lingam$Bpruned)) fit_lingam$Bpruned
    else if (!is.null(fit_lingam$B)) fit_lingam$B
    else stop("No 'B' or 'Bpruned' found in pcalg::lingam output.")

    A <- matrix(0L, nrow = p, ncol = p)
    A[abs(Bm) > 0] <- 1L
    diag(A) <- 0L
    colnames(A) <- rownames(A) <- colnames(data)
    amat <- A
  }

  if (is.logical(amat)) {
    stop("Provide a directed adjacency matrix 'amat' or set directed=TRUE to learn it via LiNGAM.")
  }
  colnames(amat) <- rownames(amat) <- colnames(data)

  # 2) Correlation matrix
  corr_Y <- stats::cor(data, use = "pairwise.complete.obs")

  # 3) Direct contributions
  Direct <- matrix(0, p, p, dimnames = list(colnames(data), colnames(data)))

  # Only Genizi is supported
  if (!identical(tolower(method), "genizi")) {
    stop("Only method = 'genizi' is supported in this package version.")
  }
  for (i in seq_len(p)) {
    Direct[i, ] <- .direct_row_genizi_fast(i, corr_Y, amat)
  }

  # 4) Standardized betas and indirect contributions
  Bstd <- compute_standardized_betas(corr_Y, amat)

  Indirect <- .indirect_from_powers_abs(Bstd)
  colnames(Indirect) <- rownames(Indirect) <- colnames(data)

  # 5) Aggregate network measures
  Total <- Direct + Indirect

  to_total   <- colSums(Total)
  from_total <- rowSums(Total)
  net_total  <- to_total - from_total
  npdc_total <- Total - t(Total)

  row_R2_direct   <- rowSums(Direct)
  row_R2_indirect <- rowSums(Indirect)
  row_R2_total    <- rowSums(Total)

  tci_direct   <- if (directed) sum(row_R2_direct)   / (p - 1) else sum(row_R2_direct)   / p
  tci_indirect <- if (directed) sum(row_R2_indirect) / (p - 1) else sum(row_R2_indirect) / p
  tci_total    <- if (directed) sum(row_R2_total)    / (p - 1) else sum(row_R2_total)    / p

  mult <- max(rowSums(Total))
  if (mult != 0) {
    tci_direct   <- tci_direct / mult
    tci_indirect <- tci_indirect / mult
    tci_total    <- tci_total / mult
  }

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
    mult           = mult
  )
}
