#' Generate data from an SVAR/VAR data-generating process
#'
#' @param n Integer, number of observations returned (after burn-in).
#' @param A0 Structural matrix (k x k). If NULL, identity is used.
#' @param list_A List of lag coefficient matrices A1..Ap (each k x k).
#' @param df Degrees of freedom for shocks. Use Inf for Gaussian.
#' @param burnin Integer burn-in length.
#' @param sd Scalar shock SD if sd_list is NULL.
#' @param sd_list Optional numeric vector of length k for shock SDs.
#' @return A time series object (ts) with n rows and k columns.
#' @export
gen_svar <- function(n,
                     A0 = NULL,
                     list_A = list(),
                     df = 10,
                     burnin = 100,
                     sd = 1,
                     sd_list = NULL) {

  # --- Basic checks
  if (!is.numeric(n) || length(n) != 1L || n < 1) stop("`n` must be a single integer >= 1.")
  if (!is.numeric(burnin) || length(burnin) != 1L || burnin < 0) stop("`burnin` must be a single integer >= 0.")
  if (!is.numeric(df) || length(df) != 1L) stop("`df` must be a single numeric value.")
  if (!is.infinite(df) && df <= 2) stop("For Student-t shocks, require df > 2 so variance is finite.")
  if (!is.numeric(sd) || length(sd) != 1L || sd < 0) stop("`sd` must be a single numeric value >= 0.")

  # --- Determine system dimension (k)
  if (length(list_A) > 0) {
    if (!is.matrix(list_A[[1]])) stop("`list_A` must be a list of matrices.")
    k <- nrow(list_A[[1]])
  } else if (!is.null(A0)) {
    if (!is.matrix(A0)) stop("`A0` must be a matrix.")
    k <- nrow(A0)
  } else {
    k <- 1L
  }

  # --- Default A0 if missing
  if (is.null(A0)) {
    A0 <- diag(k)
  }

  # --- Validate A0
  if (!is.matrix(A0) || nrow(A0) != k || ncol(A0) != k) {
    stop("`A0` must be a square (k x k) matrix with k consistent with `list_A` (if provided).")
  }
  if (isTRUE(any(!is.finite(A0)))) stop("`A0` must contain only finite values.")

  A0_inv <- tryCatch(solve(A0), error = function(e) NULL)
  if (is.null(A0_inv)) stop("`A0` is singular and cannot be inverted.")

  # --- Validate list_A (if provided)
  p <- length(list_A)
  if (p > 0) {
    ok_dims <- vapply(list_A, function(M) is.matrix(M) && all(dim(M) == c(k, k)), logical(1))
    if (!all(ok_dims)) stop("All elements of `list_A` must be (k x k) matrices.")
    if (isTRUE(any(vapply(list_A, function(M) any(!is.finite(M)), logical(1))))) {
      stop("All matrices in `list_A` must contain only finite values.")
    }
  }

  # --- Sanity checks for sd_list
  if (!is.null(sd_list)) {
    if (!is.numeric(sd_list)) stop("`sd_list` must be numeric.")
    if (length(sd_list) != k) stop("`sd_list` must have length equal to system dimension k.")
    if (any(sd_list < 0)) stop("`sd_list` must contain only non-negative values.")
    sd_vec <- as.numeric(sd_list)
  } else {
    sd_vec <- rep(sd, k)
  }

  total_n <- n + burnin

  # --- Generate structural shocks ε_t
  # df = Inf  -> Gaussian
  # df finite -> Student-t standardized to Var = 1
  if (is.infinite(df)) {
    eps_raw <- matrix(stats::rnorm(total_n * k), nrow = total_n, ncol = k)
  } else {
    t_raw <- matrix(stats::rt(total_n * k, df = df), nrow = total_n, ncol = k)
    scale_to_var1 <- sqrt(df / (df - 2))
    eps_raw <- t_raw / scale_to_var1
  }

  # --- Scale by standard deviations
  epsilon <- sweep(eps_raw, 2, sd_vec, `*`)

  # --- Initialize y
  y <- matrix(NA_real_, nrow = total_n, ncol = k)
  colnames(y) <- paste0("Y", seq_len(k))

  init_len <- max(1L, p)
  y[seq_len(init_len), ] <- matrix(stats::rnorm(init_len * k), nrow = init_len, ncol = k)

  # --- Simulate SVAR
  # Safe start index (avoids weird sequences when p=0 and keeps indices valid)
  start_i <- max(p + 1L, init_len + 1L)

  if (start_i <= total_n) {
    for (i in start_i:total_n) {
      y_sum <- rep(0, k)

      if (p > 0) {
        for (j in seq_len(p)) {
          y_sum <- y_sum + list_A[[j]] %*% y[i - j, ]
        }
      }

      y[i, ] <- A0_inv %*% (y_sum + epsilon[i, ])
    }
  }

  # --- Drop burn-in
  stats::ts(y[(burnin + 1L):total_n, , drop = FALSE])
}




