theoretical_network <- function(A0, ..., eps_sd, n.ahead, var_names = NULL) {
  A_lags <- list(...)

  if (length(A_lags) == 1L && is.list(A_lags[[1L]]) && !is.matrix(A_lags[[1L]])) {
    A_lags <- A_lags[[1L]]
  }

  if (!is.matrix(A0) || nrow(A0) != ncol(A0)) {
    stop("`A0` must be a square matrix.")
  }

  n <- nrow(A0)
  p <- length(A_lags)

  if (!is.numeric(eps_sd) || length(eps_sd) != n) {
    stop("`eps_sd` must be a numeric vector of length ncol(A0).")
  }

  if (!is.numeric(n.ahead) || length(n.ahead) != 1L || n.ahead < 0 || n.ahead %% 1 != 0) {
    stop("`n.ahead` must be a single non-negative integer.")
  }

  for (j in seq_len(p)) {
    if (!is.matrix(A_lags[[j]]) || any(dim(A_lags[[j]]) != c(n, n))) {
      stop(sprintf("A%d must be a %d x %d matrix.", j, n, n))
    }
  }

  if (is.null(var_names)) {
    var_names <- colnames(A0)
    if (is.null(var_names)) {
      var_names <- paste0("Y", seq_len(n))
    }
  }

  colnames(A0) <- rownames(A0) <- var_names

  if (p > 0) {
    for (j in seq_len(p)) {
      colnames(A_lags[[j]]) <- rownames(A_lags[[j]]) <- var_names
    }
    names(A_lags) <- paste0("A", seq_len(p))
  }

  pct_non_diag <- function(M, margin = c("row", "col")) {
    margin <- match.arg(margin)

    if (margin == "row") {
      total <- rowSums(M)
      non_diag <- total - diag(M)
      out <- ifelse(total > 0, 100 * non_diag / total, NA_real_)
      names(out) <- rownames(M)
      return(out)
    }

    total <- colSums(M)
    non_diag <- total - diag(M)
    out <- ifelse(total > 0, 100 * non_diag / total, NA_real_)
    names(out) <- colnames(M)
    out
  }

  safe_row_normalize <- function(M) {
    rs <- rowSums(M)
    out <- M
    ok <- rs > 0

    if (any(ok)) {
      out[ok, ] <- M[ok, , drop = FALSE] / rs[ok]
    }
    if (any(!ok)) {
      out[!ok, ] <- 0
    }

    attr(out, "row_divisor") <- rs
    out
  }

  summarize_irf_network <- function(M) {
    mult <- max(rowSums(M))
    table <- if (is.finite(mult) && mult > 0) M / mult else M

    from <- pct_non_diag(M, "row")
    to <- pct_non_diag(M, "col")

    list(
      mult = mult,
      divisor = mult,
      table_raw = M,
      table = table,
      from = from,
      to = to,
      tci = mean(from, na.rm = TRUE),
      net = to - from
    )
  }

  summarize_dy_network <- function(M) {
    table <- safe_row_normalize(M)
    row_divisor <- attr(table, "row_divisor")

    from <- pct_non_diag(table, "row")
    to <- pct_non_diag(table, "col")

    list(
      mult = NA_real_,
      divisor = row_divisor,
      row_divisor = row_divisor,
      table_raw = M,
      table = table,
      from = from,
      to = to,
      tci = mean(from, na.rm = TRUE),
      net = to - from
    )
  }

  sum_matrices <- function(mats) {
    if (length(mats) == 0L) {
      return(matrix(0, n, n, dimnames = list(var_names, var_names)))
    }
    Reduce(`+`, mats)
  }

  A0_inv <- solve(A0)
  D_eps <- diag(eps_sd, nrow = n, ncol = n)
  Sigma_eps <- diag(eps_sd^2, nrow = n, ncol = n)

  Phi <- if (p > 0) {
    lapply(A_lags, function(Aj) A0_inv %*% Aj)
  } else {
    list()
  }

  IRF <- array(
    0,
    dim = c(n, n, n.ahead + 1L),
    dimnames = list(
      response = var_names,
      shock = var_names,
      horizon = paste0("h=", 0:n.ahead)
    )
  )

  IRF[, , 1L] <- A0_inv %*% D_eps

  if (p > 0 && n.ahead > 0) {
    for (h in 1:n.ahead) {
      tmp <- matrix(0, n, n, dimnames = list(var_names, var_names))
      for (j in 1:min(h, p)) {
        tmp <- tmp + Phi[[j]] %*% IRF[, , h - j + 1L]
      }
      IRF[, , h + 1L] <- tmp
    }
  }

  abs_irf_list <- lapply(seq_len(n.ahead + 1L), function(k) abs(IRF[, , k]))
  irf_cumulative_abs <- sum_matrices(abs_irf_list)

  result <- summarize_irf_network(irf_cumulative_abs)

  cont <- list(
    total = abs(A0_inv),
    direct = abs(A0),
    undirect = abs(A0_inv) - abs(A0),
    indirect = abs(A0_inv) - abs(A0),
    impact = abs(IRF[, , 1L]),
    impact_scaled = abs(IRF[, , 1L])
  )

  result$contamperanous <- cont
  result$contemporaneous <- cont

  zero_mat <- matrix(0, n, n, dimnames = list(var_names, var_names))

  if (n.ahead > 0) {
    lagged_total <- irf_cumulative_abs - abs(IRF[, , 1L])
    lagged_direct <- abs(IRF[, , 2L])
    lagged_undirect <- irf_cumulative_abs - abs(IRF[, , 1L]) - abs(IRF[, , 2L])
  } else {
    lagged_total <- zero_mat
    lagged_direct <- zero_mat
    lagged_undirect <- zero_mat
  }

  result$lagged <- list(
    total = lagged_total,
    direct = lagged_direct,
    undirect = lagged_undirect,
    indirect = lagged_undirect
  )

  dy_numerators_by_h <- lapply(seq_len(n.ahead + 1L), function(k) IRF[, , k]^2)
  dy_total_num <- sum_matrices(dy_numerators_by_h)

  dy_result <- summarize_dy_network(dy_total_num)

  dy_cont_num <- IRF[, , 1L]^2

  if (n.ahead > 0) {
    dy_lagged_total_num <- sum_matrices(lapply(2:(n.ahead + 1L), function(k) IRF[, , k]^2))
    dy_lagged_direct_num <- IRF[, , 2L]^2
  } else {
    dy_lagged_total_num <- zero_mat
    dy_lagged_direct_num <- zero_mat
  }

  if (n.ahead > 1) {
    dy_lagged_undirect_num <- sum_matrices(lapply(3:(n.ahead + 1L), function(k) IRF[, , k]^2))
  } else {
    dy_lagged_undirect_num <- zero_mat
  }

  dy_result$contamperanous <- list(
    total_raw = dy_cont_num,
    total = safe_row_normalize(dy_cont_num)
  )

  dy_result$contemporaneous <- dy_result$contamperanous

  dy_result$lagged <- list(
    total_raw = dy_lagged_total_num,
    total = safe_row_normalize(dy_lagged_total_num),
    direct_raw = dy_lagged_direct_num,
    direct = safe_row_normalize(dy_lagged_direct_num),
    undirect_raw = dy_lagged_undirect_num,
    undirect = safe_row_normalize(dy_lagged_undirect_num),
    indirect_raw = dy_lagged_undirect_num,
    indirect = safe_row_normalize(dy_lagged_undirect_num)
  )

  result$dy <- dy_result

  result$svar_params <- list(
    A0 = A0,
    var_params = A_lags,
    reduced_form_params = Phi,
    eps_sd = stats::setNames(eps_sd, var_names),
    Sigma_eps = stats::setNames(as.list(eps_sd^2), var_names)
  )

  result$IRF <- IRF
  result$IRF_cumulative_abs <- irf_cumulative_abs

  result
}