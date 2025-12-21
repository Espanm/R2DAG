#' Rolling IRF network
#'
#' @param data data.frame or matrix
#' @param block_size Integer rolling window size
#' @param n.ahead Integer horizon
#' @param cumsum Logical; if TRUE uses cumulative IRFs
#' @param amat FALSE or adjacency matrix passed to irf_network
#' @return A list containing a rolling TCI vector and rolling TO/FROM data frames.
#' @export
rolling_irf_network <- function(data, block_size, n.ahead, cumsum = TRUE, amat = FALSE) {

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
  if (!is.numeric(n.ahead) || length(n.ahead) != 1L || n.ahead < 0L) {
    stop("`n.ahead` must be a single integer >= 0.")
  }

  tci_vector <- numeric(0)

  to_df <- data.frame(matrix(ncol = p, nrow = 0))
  colnames(to_df) <- colnames(data)

  from_df <- data.frame(matrix(ncol = p, nrow = 0))
  colnames(from_df) <- colnames(data)

  for (i in 1:(n - block_size)) {
    rolling_data <- data[i:(i + block_size), , drop = FALSE]

    var_model <- vars::VAR(rolling_data, p = 1, type = "const")

    network <- irf_network(var_model, n.ahead = n.ahead, cumsum = cumsum, amat = amat)

    tci_vector <- c(network$tci, tci_vector)
    to_df <- rbind(to_df, as.data.frame(as.list(network$to)))
    from_df <- rbind(from_df, as.data.frame(as.list(network$from)))
  }

  list(
    tci_vector = tci_vector,
    to_df = to_df,
    from_df = from_df
  )
}
