#' Theoretical values of the R2 network
#'
#' @param A0 Structural matrix (p x p)
#' @param var_names Variable names
#' @param directed logical
#' @return A list with theoretical network tables and measures
#' @export
theoretical_R2_network <- function(A0,
                                   sigma_eps = rep(1, nrow(A0)),
                                   var_names = paste0("X", 1:nrow(A0)),
                                   directed = TRUE) {

  p <- nrow(A0)
  stopifnot(
    is.matrix(A0),
    ncol(A0) == p,
    length(sigma_eps) == p
  )

  # --- invert A0 ---
  A0_inv <- tryCatch(solve(A0), error = function(e) NULL)
  if (is.null(A0_inv)) stop("`A0` is singular and cannot be inverted.")

  # --- structural shock covariance ---
  Sigma_eps <- diag(sigma_eps^2)

  # --- implied covariance & correlation of Y ---
  cov_Y  <- A0_inv %*% Sigma_eps %*% t(A0_inv)
  corr_Y <- stats::cov2cor(cov_Y)

  # --- adjacency matrix ---
  if (directed) {
    amat <- matrix(0L, p, p)
    for (i in seq_len(p)) {
      for (j in seq_len(p)) {
        if (i != j && A0[i, j] != 0) amat[i, j] <- 1L
      }
    }
  } else {
    amat <- matrix(1L, p, p)
    diag(amat) <- 0L
  }
  colnames(amat) <- rownames(amat) <- var_names

  # --- Direct Genizi R2 ---
  Direct <- matrix(0, p, p, dimnames = list(var_names, var_names))
  for (i in seq_len(p)) {
    Direct[i, ] <- .direct_row_genizi_fast(i, corr_Y, amat)
  }

  # --- Indirect effects: only if directed=TRUE ---
  if (directed) {
    Bstd <- compute_standardized_betas(corr_Y, amat) 

    Indirect <- .indirect_from_powers_abs(Bstd) / 2
    colnames(Indirect) <- rownames(Indirect) <- var_names
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
    tci_total    <- sum(row_R2_total)  / p
  }

  # --- Aggregate measures (based on Total) ---
  to_total   <- colSums(Total)
  from_total <- rowSums(Total)
  net_total  <- to_total - from_total
  npdc_total <- Total - t(Total)

  mult <- max(rowSums(Total))

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
    cov_Y          = cov_Y,
    sigma_eps      = sigma_eps
  )
}
