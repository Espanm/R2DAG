gen_var <- function(n, list_A, sigma) {

  k <- nrow(list_A[[1]])    # Number of variables (k)
  p <- length(list_A)       # Number of lags (p)

  # Cholesky decomposition to generate correlated innovations
  chol_sigma <- chol(sigma)
  innovations <- matrix(rnorm(n * k), nrow = n) %*% chol_sigma

  # Initialize the time series matrix
  y <- matrix(NA, nrow = n, ncol = k)

  # Initial values: random for the first 'p' rows
  y[1:p, ] <- matrix(rnorm(p * k), nrow = p)

  # Simulate the VAR process
  for (i in (p+1):n) {
    y[i, ] <- 0  # Initialize the row
    for (j in 1:p) {
      y[i, ] <- y[i, ] + list_A[[j]] %*% y[i-j, ]
    }
    y[i, ] <- y[i, ] + innovations[i, ]  # Add the innovations (errors)
  }

  # Convert the result into a time series object
  y_ts <- ts(y)

  return(y_ts)
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
    A_list[[j]] <- A_j
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

  return(IRF)
}

cumsum_irf <- function(IRF, h) {
  # Ellenőrzés: h nem lehet nagyobb, mint a dimenzió
  if (h > dim(IRF)[3]) {
    stop("h is larger than the number of horizons in IRF.")
  }

  # Összegzés k = 1..h (azaz IRF[,,2] .. IRF[,,h+1])
  A <- apply(IRF[,,2:(h + 1), drop = FALSE], c(1, 2), sum)

  return(A)
}


