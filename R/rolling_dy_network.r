#' Rolling Diebold-Yilmaz connectedness (TCI)
#'
#' @param data data.frame or matrix
#' @param block_size Integer rolling window size
#' @param n.ahead Integer horizon (forecast steps used in FEVD)
#' @return A list containing a rolling TCI vector.
#' @export
rolling_dy_network <- function(data, block_size, n.ahead) {

  stopifnot(is.data.frame(data) || is.matrix(data))
  data <- as.data.frame(data)
  stopifnot(!is.null(colnames(data)))

  n <- nrow(data)
  p <- ncol(data)

  if (!is.numeric(block_size) || length(block_size) != 1L || block_size < 2L) {
    stop("`block_size` must be a single integer >= 2.")
  }
  if (block_size >= n) {
    stop("`block_size` must be smaller than the number of rows in `data`.")
  }
  if (!is.numeric(n.ahead) || length(n.ahead) != 1L || n.ahead < 1L) {
    stop("`n.ahead` must be a single integer >= 1.")
  }

  tci_vector <- numeric(0)

  for (i in 1:(n - block_size)) {
    rolling_data <- data[i:(i + block_size), , drop = FALSE]

    var_model <- vars::VAR(rolling_data, p = 1, type = "const")
    B_hat <- vars::Bcoef(var_model)

    R <- stats::resid(var_model)
    Sigma_hat <- stats::cov(R)

    TCi <- FEVD_TCI(Phi = B_hat, Sigma = Sigma_hat, nfore = n.ahead)
    tci_vector <- c(TCi, tci_vector)
  }

  list(tci_vector = tci_vector)
}


