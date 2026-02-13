#' Rolling window R2 network analysis
#'
#' @param data data.frame or matrix
#' @param window_size integer window length
#' @param step integer step size
#' @param bootstrap number of bootstrap replications
#' @param amat_list FALSE or list of adjacency matrices
#' @param ... passed to R2_network (incl. directed = TRUE/FALSE)
#' @return A list of rolling network measures and (optionally) tables
#' @export
rolling_R2_network <- function(data,
                               window_size,
                               step = 1,
                               bootstrap = 1,
                               amat_list = FALSE,
                               ...) {

  stopifnot(is.data.frame(data) || is.matrix(data))
  data <- as.data.frame(data)
  stopifnot(!is.null(colnames(data)))

  n <- nrow(data)
  p <- ncol(data)

  if (window_size <= 1 || window_size > n) stop("Invalid window_size.")
  if (step < 1) stop("step must be >= 1.")
  if (bootstrap < 1) stop("bootstrap must be >= 1.")

  dots <- list(...)
  directed_flag <- if (!is.null(dots$directed)) isTRUE(dots$directed) else TRUE

  starts <- seq.int(1L, n - window_size + 1L, by = step)
  ends   <- starts + window_size - 1L
  nW     <- length(starts)

  tci_direct   <- rep(NA_real_, nW)
  tci_indirect <- rep(NA_real_, nW)
  tci_total    <- rep(NA_real_, nW)

  # --- NEW: window-by-node NET (TO-FROM) ---
  net_total <- matrix(NA_real_, nrow = nW, ncol = p,
                      dimnames = list(NULL, colnames(data)))

  amat_mode <- vector("list", nW)

  direct_tables   <- if (directed_flag) vector("list", nW) else NULL
  indirect_tables <- if (directed_flag) vector("list", nW) else NULL

  # --- Safe scalar extractor ---
  get_num1 <- function(x, name) {
    v <- x[[name]]
    if (is.null(v) || length(v) == 0) return(NA_real_)
    as.numeric(v[1])
  }

  # --- Safe vector extractor (length p) ---
  get_vecp <- function(x, name) {
    v <- x[[name]]
    if (is.null(v) || length(v) != p) return(rep(NA_real_, p))
    as.numeric(v)
  }

  # --- Single run ---
  one_run <- function(dat_sample, amat_override = FALSE) {
    tryCatch(
      if (isFALSE(amat_override)) {
        R2_network(dat_sample, ...)
      } else {
        R2_network(dat_sample, amat = amat_override, ...)
      },
      error = function(e) NULL
    )
  }

  # --- Align adjacency ---
  align_amat <- function(M) {
    if (is.null(M) || !is.matrix(M)) return(FALSE)
    if (!all(dim(M) == c(p, p))) return(FALSE)
    colnames(M) <- rownames(M) <- colnames(data)
    storage.mode(M) <- "integer"
    M
  }

  use_fixed_amats <- is.list(amat_list)
  if (use_fixed_amats) {
    if (length(amat_list) < nW) stop("amat_list shorter than number of windows.")
    amat_list <- lapply(amat_list[seq_len(nW)], align_amat)
  }

  key_of <- function(M) paste(c(M), collapse = ",")

  for (k in seq_len(nW)) {

    dat_k <- data[starts[k]:ends[k], , drop = FALSE]
    res_list <- vector("list", bootstrap)

    for (b in seq_len(bootstrap)) {

      dat_b <- if (bootstrap > 1) {
        dat_k[sample.int(nrow(dat_k), nrow(dat_k), replace = TRUE), , drop = FALSE]
      } else dat_k

      if (use_fixed_amats) {
        M <- amat_list[[k]]
        res_list[[b]] <- one_run(dat_b, amat_override = M)
      } else {
        res_list[[b]] <- one_run(dat_b, amat_override = FALSE)
      }
    }

    ok <- Filter(Negate(is.null), res_list)

    if (length(ok) == 0) {
      amat_mode[[k]] <- NULL
      if (directed_flag) {
        direct_tables[[k]] <- NULL
        indirect_tables[[k]] <- NULL
      }
      next
    }

    # --- Scalar averages ---
    tci_direct[k]   <- mean(vapply(ok, get_num1, numeric(1), "tci_direct"),   na.rm = TRUE)
    tci_indirect[k] <- mean(vapply(ok, get_num1, numeric(1), "tci_indirect"), na.rm = TRUE)
    tci_total[k]    <- mean(vapply(ok, get_num1, numeric(1), "tci_total"),    na.rm = TRUE)

    # --- NEW: NET averages (node-wise) ---
    net_mat <- vapply(ok, get_vecp, numeric(p), name = "net_total")  # p x B
    net_total[k, ] <- rowMeans(net_mat, na.rm = TRUE)

    # --- Direct / Indirect tables only if directed ---
    if (directed_flag) {

      ok_tab <- Filter(function(x)
        is.matrix(x$direct_table) && is.matrix(x$indirect_table), ok)

      if (length(ok_tab) > 0) {

        Direct_arr <- simplify2array(lapply(ok_tab, `[[`, "direct_table"))
        Indir_arr  <- simplify2array(lapply(ok_tab, `[[`, "indirect_table"))

        Direct_mean <- apply(Direct_arr, c(1,2), mean, na.rm = TRUE)
        Indir_mean  <- apply(Indir_arr,  c(1,2), mean, na.rm = TRUE)

        colnames(Direct_mean) <- rownames(Direct_mean) <- colnames(data)
        colnames(Indir_mean)  <- rownames(Indir_mean)  <- colnames(data)

        direct_tables[[k]]   <- Direct_mean
        indirect_tables[[k]] <- Indir_mean
      } else {
        direct_tables[[k]]   <- NULL
        indirect_tables[[k]] <- NULL
      }
    }

    # --- amat mode ---
    amats_ok <- lapply(ok, function(x) x$amat)
    amats_ok <- Filter(function(M)
      is.matrix(M) && all(dim(M) == c(p,p)), amats_ok)

    if (length(amats_ok) > 0) {
      keys <- vapply(amats_ok, key_of, character(1))
      key_mode <- names(sort(table(keys), decreasing = TRUE))[1]
      idx_mode <- match(key_mode, keys)
      amat_mode[[k]] <- amats_ok[[idx_mode]]
    } else {
      amat_mode[[k]] <- NULL
    }
  }

  out <- list(
    tci_direct   = tci_direct,
    tci_indirect = tci_indirect,
    tci_total    = tci_total,
    net_total    = net_total,   # <-- NEW: nW x p, minden ablakra minden idősor NET-je
    amat_mode    = amat_mode
  )

  if (directed_flag) {
    out$direct_tables   <- direct_tables
    out$indirect_tables <- indirect_tables
  }

  out
}
