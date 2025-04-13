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
  A_inv <- solve(svar_model[1]$A)

  IRF <- calculate_irf(var_model, n.ahead, ortho = FALSE, shock = A_inv)

  if (cumsum){irf_matrix <- cumsum_irf(IRF, n.ahead)}
  else {irf_matrix <- IRF[,,(n.ahead+1)]}

  colnames(irf_matrix) <- rownames(irf_matrix) <- colnames(data)

  return(irf_matrix)
}
