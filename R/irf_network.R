irf_network <- function(var_model, n.ahead, cumsum=TRUE, amat=FALSE){

  data <- var_model$y
  # Get the number of variables
  p <- ncol(data)

  if ((is.logical(amat) && !amat)) {
    lingam_amat <- data2amat(data)
    amat <- lingam_amat$amat
  }

  if (!is.logical(amat)){
    # Assign column names to adjacency matrix
    colnames(amat) <- rownames(amat) <- colnames(data)
  }
  paths <- compute_reachability(amat)

  # 1. Változók nevei
  vars <- names(paths)

  # 2. Inicializáljuk a k × k nullmátrixot
  paths_matrix <- matrix(0, nrow = p, ncol = p)
  rownames(paths_matrix) <- colnames(paths_matrix) <- vars

  # 3. Feltöltés: A[i, j] = 1, ha paths[[i]] tartalmazza j-t (i → j)
  for (i in vars) {
    targets <- paths[[i]]
    paths_matrix[i, targets] <- 1
  }


  # Kezdetben: 1-ek a diagonálban, 0 mindenhol máshol
  Amat <- construct_svar_template(amat)

  svar_model <- SVAR(x = var_model, estmethod = "scoring", Amat = Amat, Bmat = NULL, max.iter = 1000, maxls = 1000, conv.crit = 1.0e-8)
  A <- svar_model[1]$A
  A_inv <- solve(A)

  IRF <- calculate_irf(var_model, n.ahead, ortho = FALSE, shock = A_inv)

  if (cumsum){irf_matrix <- cumsum_irf(IRF, n.ahead)}
  else {irf_matrix <- IRF[,,(n.ahead+1)]}

  colnames(irf_matrix) <- rownames(irf_matrix) <- colnames(data)

  result <- list()

  result$mult <- max(rowSums(irf_matrix))

  result$table <- irf_matrix / max(rowSums(irf_matrix))

  result$from <- calculate_non_diag_percentage(irf_matrix, "col")
  names(result$from) <- colnames(data)

  result$to <- calculate_non_diag_percentage(irf_matrix, "row")
  names(result$to) <- colnames(data)

  result$tci <- mean(result$from)
  result$contamperanous$total <- abs(A_inv)
  result$contamperanous$direct <- abs(A)
  result$contamperanous$undirect <- abs(A_inv) - abs(A)
  if (n.ahead > 0){
    result$lagged$total <- irf_matrix - IRF[,,1]
    result$lagged$direct  <- IRF[,,2]
    result$lagged$undirect <- irf_matrix - IRF[,,2] - IRF[,,1]
  }
  svar_params <- list()
  svar_params$A0 <- A
  svar_params$var_params <- lapply(extract_lag_matrices(var_model), function(x) A %*% x)

  result$svar_params <- svar_params

  return(result)
}

off_diag_percentage <- function(mat) {
  if (!is.matrix(mat)) stop("Input must be a matrix.")

  total_sum <- sum(mat)
  diag_sum <- sum(diag(mat))

  if (total_sum == 0) return(NA)  # avoid division by zero

  percentage <- 1 - (diag_sum / total_sum)
  return(percentage)
}

calculate_non_diag_percentage <- function(mat, direction = "row") {
  # Ellenőrizzük, hogy a bemenet mátrix legyen
  if (!is.matrix(mat)) stop("Input must be a matrix.")

  # Ellenőrizzük, hogy a direction "row" vagy "col" legyen
  if (!(direction %in% c("row", "col"))) stop("Direction must be either 'row' or 'col'.")

  # Sorok és oszlopok számának lekérése
  nrow_mat <- nrow(mat)
  ncol_mat <- ncol(mat)

  # Készítünk egy vektort a nem diagonális elemek százalékos arányának
  result <- numeric()

  if (direction == "row") {
    # Sorok szerint számolunk
    for (i in 1:nrow_mat) {
      # Sor összegének kiszámítása
      row_sum <- sum(mat[i, ])

      # Nem diagonális elemek összegének kiszámítása
      non_diag_sum <- 0
      for (j in 1:ncol_mat) {
        if (i != j) {  # Ha nem diagonális elem
          non_diag_sum <- non_diag_sum + mat[i, j]
        }
      }

      # Százalékos arány: nem diagonális elemek / sor összegéhez
      result[i] <- (non_diag_sum / row_sum)
    }


  } else {
    # Oszlopok szerint számolunk
    for (j in 1:ncol_mat) {
      # Oszlop összegének kiszámítása
      col_sum <- sum(mat[, j])

      # Nem diagonális elemek összegének kiszámítása
      non_diag_sum <- 0
      for (i in 1:nrow_mat) {
        if (i != j) {  # Ha nem diagonális elem
          non_diag_sum <- non_diag_sum + mat[i, j]
        }
      }

      # Százalékos arány: nem diagonális elemek / oszlop összegéhez
      result[j] <- (non_diag_sum / col_sum)
    }

  }

  return(result)
}

