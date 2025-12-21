#' IRF-based network decomposition from a VAR/SVAR
#'
#' Builds an SVAR identification (A matrix) consistent with an adjacency matrix
#' and computes an IRF-based connectedness matrix at horizon n.ahead.
#'
#' @param var_model A VAR model object (e.g., from vars::VAR)
#' @param n.ahead Integer horizon
#' @param cumsum Logical; if TRUE uses cumulative IRFs up to n.ahead
#' @param amat FALSE or an adjacency matrix (child x parent). If FALSE, learned from data.
#' @return A list with normalized IRF table and summary measures.
#' @export
irf_network <- function(var_model, n.ahead, cumsum = TRUE, amat = FALSE) {

  data <- var_model$y
  p <- ncol(data)

  if (is.logical(amat) && !amat) {
    lingam_amat <- data2amat(data)
    amat <- lingam_amat$amat
  }

  if (is.logical(amat)) {
    stop("`amat` must be a matrix (or set amat = FALSE to learn it).")
  }

  colnames(amat) <- rownames(amat) <- colnames(data)

  # SVAR template from adjacency
  Amat <- construct_svar_template(amat)

  # Use explicit namespace to satisfy R CMD check
  svar_model <- vars::SVAR(
    x = var_model,
    estmethod = "scoring",
    Amat = Amat,
    Bmat = NULL,
    max.iter = 1000,
    maxls = 1000,
    conv.crit = 1.0e-8
  )

  A <- svar_model[1]$A
  A_inv <- solve(A)

  IRF <- calculate_irf(var_model, n.ahead, ortho = FALSE, shock = A_inv)

  irf_matrix <- if (isTRUE(cumsum)) cumsum_irf(IRF, n.ahead) else IRF[, , (n.ahead + 1)]
  colnames(irf_matrix) <- rownames(irf_matrix) <- colnames(data)

  result <- list()

  mult <- max(rowSums(irf_matrix))
  result$mult <- mult
  result$table <- if (is.finite(mult) && mult > 0) irf_matrix / mult else irf_matrix

  result$from <- calculate_non_diag_percentage(irf_matrix, "row")
  names(result$from) <- colnames(data)

  result$to <- calculate_non_diag_percentage(irf_matrix, "col")
  names(result$to) <- colnames(data)

  result$tci <- mean(result$from, na.rm = TRUE)
  result$net <- result$to - result$from

  result$contamperanous$total <- abs(A_inv)
  result$contamperanous$direct <- abs(A)
  result$contamperanous$undirect <- abs(A_inv) - abs(A)

  if (n.ahead > 0) {
    result$lagged$total <- irf_matrix - IRF[, , 1]
    result$lagged$direct <- IRF[, , 2]
    result$lagged$undirect <- irf_matrix - IRF[, , 2] - IRF[, , 1]
  }

  svar_params <- list()
  svar_params$A0 <- A
  svar_params$var_params <- lapply(extract_lag_matrices(var_model), function(x) A %*% x)
  result$svar_params <- svar_params

  result
}
