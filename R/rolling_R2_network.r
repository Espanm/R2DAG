# =========================
# Rolling wrapper (no parallel) + optional amat_list
# Now returns Direct/Indirect tables per window
# =========================

#' Rolling window R2 network analysis
#'
#' @param data data.frame or matrix
#' @param window_size integer window length
#' @param step integer step size
#' @param bootstrap number of bootstrap replications
#' @param amat_list FALSE or list of adjacency matrices
#' @param ... passed to R2_network
#' @return A list of rolling network measures and tables
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

  if (!is.numeric(window_size) || length(window_size) != 1 || window_size <= 1)
    stop("`window_size` must be a single integer > 1.")
  if (window_size > n)
    stop("`window_size` is larger than the number of rows in `data`.")
  if (!is.numeric(step) || length(step) != 1 || step < 1)
    stop("`step` must be a single integer >= 1.")
  if (!is.numeric(bootstrap) || length(bootstrap) != 1 || bootstrap < 1)
    stop("`bootstrap` must be a single integer >= 1.")

  starts <- seq.int(1L, n - window_size + 1L, by = step)
  ends   <- starts + window_size - 1L
  nW     <- length(starts)

  tci_direct   <- rep(NA_real_, nW)
  tci_indirect <- rep(NA_real_, nW)
  tci_total    <- rep(NA_real_, nW)
  mult         <- rep(NA_real_, nW)

  amat_mode    <- vector("list", nW)

  # NEW: store per-window (averaged) Direct/Indirect tables
  direct_tables   <- vector("list", nW)
  indirect_tables <- vector("list", nW)

  # Single run (one bootstrap replicate)
  one_run <- function(dat_sample, amat_override = NULL) {
    tryCatch(
      if (is.null(amat_override)) R2_network(dat_sample, ...)
      else R2_network(dat_sample, amat = amat_override, ...),
      error = function(e) NULL
    )
  }

  # Align and validate adjacency matrix
  align_amat <- function(M) {
    if (is.null(M) || !is.matrix(M)) return(NULL)
    if (!all(dim(M) == c(p, p))) return(NULL)

    if (!is.null(colnames(M)) && !is.null(colnames(data))) {
      if (!setequal(colnames(M), colnames(data)))
        stop("One element of `amat_list` does not contain the same variables as `data`.")
      M <- M[colnames(data), colnames(data), drop = FALSE]
    } else {
      colnames(M) <- rownames(M) <- colnames(data)
    }

    storage.mode(M) <- "integer"
    M
  }

  # If amat_list is provided: map it to windows
  use_fixed_amats <- is.list(amat_list)
  if (use_fixed_amats) {
    if (!is.null(names(amat_list))) {
      wanted <- paste0("W", seq_len(nW))
      idx <- match(wanted, names(amat_list))
      if (any(is.na(idx))) stop("`amat_list` is incomplete (provide names W1..Wn or exactly nW elements).")
      amat_list <- amat_list[idx]
    } else {
      if (length(amat_list) < nW)
        stop("`amat_list` length is smaller than the number of required windows.")
      amat_list <- amat_list[seq_len(nW)]
    }
    amat_list <- lapply(amat_list, align_amat)
  }

  # Key for mode selection (when learning structure each time)
  key_of <- function(M) paste(c(M), collapse = ",")

  for (k in seq_len(nW)) {
    i1 <- starts[k]
    i2 <- ends[k]
    dat_k <- data[i1:i2, , drop = FALSE]

    # Run bootstrap replicates
    res_list <- vector("list", bootstrap)

    for (b in seq_len(bootstrap)) {
      dat_b <- if (bootstrap > 1) {
        dat_k[sample.int(nrow(dat_k), nrow(dat_k), replace = TRUE), , drop = FALSE]
      } else {
        dat_k
      }

      if (use_fixed_amats) {
        M <- amat_list[[k]]
        if (is.null(M)) { res_list[[b]] <- NULL; next }
        res_list[[b]] <- one_run(dat_b, amat_override = M)
      } else {
        res_list[[b]] <- one_run(dat_b, amat_override = NULL)
      }
    }

    ok <- Filter(Negate(is.null), res_list)

    if (length(ok) > 0) {
      # Average scalar metrics
      tci_direct[k]   <- mean(vapply(ok, function(x) x$tci_direct,   numeric(1)), na.rm = TRUE)
      tci_indirect[k] <- mean(vapply(ok, function(x) x$tci_indirect, numeric(1)), na.rm = TRUE)
      tci_total[k]    <- mean(vapply(ok, function(x) x$tci_total,    numeric(1)), na.rm = TRUE)
      mult[k]         <- mean(vapply(ok, function(x) x$mult,         numeric(1)), na.rm = TRUE)

      # NEW: average Direct/Indirect tables across successful replicates
      Direct_arr <- simplify2array(lapply(ok, `[[`, "direct_table"))
      Indir_arr  <- simplify2array(lapply(ok, `[[`, "indirect_table"))

      Direct_mean <- apply(Direct_arr, c(1, 2), mean, na.rm = TRUE)
      Indir_mean  <- apply(Indir_arr,  c(1, 2), mean, na.rm = TRUE)

      colnames(Direct_mean) <- rownames(Direct_mean) <- colnames(data)
      colnames(Indir_mean)  <- rownames(Indir_mean)  <- colnames(data)

      direct_tables[[k]]   <- Direct_mean
      indirect_tables[[k]] <- Indir_mean

      # amat_mode logic
      if (use_fixed_amats) {
        amat_mode[[k]] <- amat_list[[k]]
      } else {
        amats_ok <- lapply(ok, function(x) x$amat)
        amats_ok <- Filter(function(M) is.matrix(M) && all(dim(M) == c(p, p)), amats_ok)

        if (length(amats_ok) > 0) {
          keys <- vapply(amats_ok, key_of, character(1))
          tab  <- table(keys)
          key_mode <- names(tab)[which.max(tab)]
          idx_mode <- match(key_mode, keys)

          Mmode <- amats_ok[[idx_mode]]
          colnames(Mmode) <- rownames(Mmode) <- colnames(data)
          storage.mode(Mmode) <- "integer"
          amat_mode[[k]] <- Mmode
        } else {
          amat_mode[[k]] <- NULL
        }
      }
    } else {
      amat_mode[[k]] <- NULL
      direct_tables[[k]] <- NULL
      indirect_tables[[k]] <- NULL
    }
  }

  list(
    tci_direct      = tci_direct,
    tci_indirect    = tci_indirect,
    tci_total       = tci_total,
    mult            = mult,
    amat_mode       = amat_mode,
    direct_tables   = direct_tables,    # per window: averaged Direct table
    indirect_tables = indirect_tables   # per window: averaged Indirect table
  )
}