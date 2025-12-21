#' Theoretical values of the R2 network
#'
#' @param A0 Structural matrix (p x p)
#' @param var_names Variable names
#' @param directed logical
#' @return A list with theoretical network tables and measures
#' @export
theoretical_R2_network <- function(A0,
                                  var_names = paste0("X", 1:nrow(A0)),
                                  directed = TRUE) {
  p <- nrow(A0)
  stopifnot(is.matrix(A0), ncol(A0) == p)

  A0_inv <- tryCatch(solve(A0), error = function(e) NULL)
  if (is.null(A0_inv)) stop("`A0` is singular and cannot be inverted.")

  cov_Y  <- A0_inv %*% t(A0_inv)
  corr_Y <- stats::cov2cor(cov_Y)

  amat <- matrix(0L, p, p)
  for (i in seq_len(p)) {
    for (j in seq_len(p)) {
      if (i != j && A0[i, j] != 0) amat[i, j] <- 1L
    }
  }
  colnames(amat) <- rownames(amat) <- var_names

  Direct <- matrix(0, p, p, dimnames = list(var_names, var_names))
  for (i in seq_len(p)) {
    Direct[i, ] <- .direct_row_genizi_fast(i, corr_Y, amat)
  }

  Bstd <- compute_standardized_betas(corr_Y, amat)

  Indirect <- .indirect_from_powers_abs(Bstd)
  colnames(Indirect) <- rownames(Indirect) <- var_names

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
