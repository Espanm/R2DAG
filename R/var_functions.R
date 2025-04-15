gen_svar <- function(n, A0 = NULL, list_A = list(), df = 10, burnin = 100) {
  # --- Determine system dimension (k)
  if (length(list_A) > 0) {
    k <- nrow(list_A[[1]])
  } else if (!is.null(A0)) {
    k <- nrow(A0)
  } else {
    k <- 1  # fallback
  }

  # --- Default A0 if missing
  if (is.null(A0)) {
    A0 <- diag(k)
  }

  p <- length(list_A)
  A0_inv <- solve(A0)

  total_n <- n + burnin  # simulate extra steps

  # --- Generate structural shocks ε_t
  if (is.infinite(df)) {
    epsilon <- matrix(rnorm(total_n * k), nrow = total_n)
  } else {
    t_raw <- matrix(rt(total_n * k, df = df), nrow = total_n)
    scale <- sqrt(df / (df - 2))  # standardize to variance 1
    epsilon <- t_raw / scale
  }

  # --- Initialize y
  y <- matrix(NA, nrow = total_n, ncol = k)
  colnames(y) <- paste0("X", seq_len(ncol(y)))
  y[1:max(1, p), ] <- matrix(rnorm(max(1, p) * k), nrow = max(1, p))

  # --- Simulate the SVAR process
  for (i in (p + 1):total_n) {
    y_sum <- rep(0, k)
    if (p > 0) {
      for (j in 1:p) {
        y_sum <- y_sum + list_A[[j]] %*% y[i - j, ]
      }
    }
    y[i, ] <- A0_inv %*% (y_sum + epsilon[i, ])
  }

  # Return only post-burnin sample as time series
  return(ts(y[(burnin + 1):total_n, ]))
}

extract_lag_matrices <- function(var_model) {
  K <- length(var_model$varresult)       # változók száma
  p <- var_model$p                       # késleltetések száma

  # Koefﬁciens-vektorok minden egyenlethez
  coef_vecs <- lapply(var_model$varresult, coef)

  varnames <- colnames(var_model$y)

  # Üres lista a lag mátrixokhoz
  A_list <- vector("list", length = p)

  # Minden lagmátrixhoz
  for (j in 1:p) {
    A_j <- matrix(NA, nrow = K, ncol = K)
    for (eq in 1:K) {
      # Összerakjuk az adott lag j mátrixát: változók × egyenletek
      coef_names <- paste0(varnames, ".", "l", j)
      A_j[, eq] <- coef_vecs[[eq]][coef_names]
    }
    A_list[[j]] <- t(A_j)
  }

  # Elnevezzük őket
  names(A_list) <- paste0("A", 1:p)

  return(A_list)
}

calculate_irf <- function(var_model, n.ahead, ortho = FALSE, shock = "None") {
  K <- length(var_model$varresult)
  p <- var_model$p
  varnames <- colnames(var_model$y)

  # Lag mátrixok lekérése
  A_list <- extract_lag_matrices(var_model)

  # Sokk mátrix: ortogonalizált (Cholesky) vagy identitás
  Sigma_u <- summary(var_model)$covres

  if (is.character(shock) && shock == "None") {
    shock_mat <- if (ortho) chol(Sigma_u) else diag(K)
  } else {
    shock_mat <- shock
  }

  # IRF tömb inicializálása: [válasz, sokk, idő]
  IRF <- array(0, dim = c(K, K, n.ahead + 1))
  IRF[, , 1] <- shock_mat  # kezdő sokkhatás

  if (n.ahead==0){
    return(abs(IRF))
  }

  # Rekurzív IRF számítás
  for (h in 1:n.ahead) {
    IRF_h <- matrix(0, nrow = K, ncol = K)
    for (j in 1:min(h, p)) {
      IRF_h <- IRF_h + A_list[[j]] %*% IRF[, , h - j + 1]
    }
    IRF[, , h + 1] <- IRF_h
  }

  # Dimenziónevek
  dimnames(IRF) <- list(
    response = varnames,
    shock = varnames,
    horizon = 0:n.ahead
  )

  return(abs(IRF))
}

library(igraph)

construct_svar_template <- function(adj) {
  k <- nrow(adj)

  # Step 1: Create graph and get causal (topological) ordering
  g <- graph_from_adjacency_matrix(adj, mode = "directed")
  causal_order <- as.numeric(topological.sort(g, mode = "in"))

  # Step 2: Construct lower-triangular matrix in causal order
  A_causal <- matrix(NA, k, k)
  diag(A_causal) <- 1
  A_causal[upper.tri(A_causal)] <- 0  # set upper triangle to 0

  # Step 3: Permute back to original order
  # Get inverse permutation
  inv_order <- order(causal_order)

  # Apply permutation to rows and columns
  A_original <- A_causal[inv_order, inv_order]

  return(A_original)
}

cumsum_irf <- function(IRF, h) {
  # Ellenőrzés: h nem lehet nagyobb, mint a dimenzió
  if (h > dim(IRF)[3]) {
    stop("h is larger than the number of horizons in IRF.")
  }

  # Méretek kinyerése
  n_row <- dim(IRF)[1]  # válaszadó változók száma
  n_col <- dim(IRF)[2]  # sokkoló változók száma

  # Üres mátrix inicializálása az eredményhez
  A <- matrix(0, nrow = n_row, ncol = n_col)

  # Összegzés ciklussal: IRF[,,1]..IRF[,,h+1]
  for (t in 1:(h + 1)) {
    A <- A + IRF[,,t]
  }

  return(A)
}



