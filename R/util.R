#' Zero-handling and filtering helper
#'
#' @param Z count matrix (samples x taxa)
#' @param imputation logical, whether to impute zeros
#' @param zero_ratio_threshold numeric, threshold for selecting taxa
#' @param zero_impute numeric, value to replace zeros
#' @return list with Z_processed and selected_taxa indices
zero_filter_impute <- function(Z, imputation = FALSE,
                               zero_ratio_threshold = 0.1,
                               zero_impute = 0.5) {
  zero_ratio <- colSums(Z == 0) / nrow(Z)
  if (imputation) {
    selected_taxa <- which(zero_ratio < zero_ratio_threshold)
    if (length(selected_taxa) == 0) {
      stop("No taxa meet the zero ratio threshold.")
    }
    Z[, selected_taxa] <- ifelse(Z[, selected_taxa] == 0,
                                 zero_impute,
                                 Z[, selected_taxa])
  } else {
    selected_taxa <- which(zero_ratio == 0)
  }
  list(Z = Z, selected_taxa = selected_taxa)
}

#' Filter taxa by minimum non‐zero counts per group
#'
#' @param Z numeric matrix (samples x taxa)
#' @param group integer or factor vector of length nrow(Z), coded 0/1 (two groups)
#' @param min_nonzero integer; minimum required non‐zero samples per group
#' @return filtered Z matrix
filter_by_prevalence <- function(Z, group, min_nonzero = 5) {
  if (length(group) != nrow(Z)) {
    stop("Length of 'group' must match nrow(Z)")
  }
  keep <- vapply(seq_len(ncol(Z)), function(j) {
    nz1 <- sum(Z[group == 0, j] != 0)
    nz2 <- sum(Z[group == 1, j] != 0)
    (nz1 >= min_nonzero) && (nz2 >= min_nonzero)
  }, logical(1))
  Z[, keep, drop = FALSE]
}

#' Generalized CLR transformation
#'
#' @param Z count matrix (samples x taxa)
#' @param selected_taxa integer vector, indices of taxa for balance
#' @param zero_impute numeric, replacement for zeros
#' @return log-ratio transformed matrix
gclr_transform <- function(Z, selected_taxa, zero_impute = 0.5) {
  p <- ncol(Z)
  D <- length(selected_taxa)
  # indicator and NA-mask
  U <- ifelse(Z == 0, 0, 1)
  U[, selected_taxa] <- 1
  NA_U <- ifelse(U == 0, NA, 1)

  # replace zeros and log
  M <- ifelse(Z == 0, zero_impute, Z)
  L <- log(M)

  # balance weights
  omega <- rep(0, p)
  omega[selected_taxa] <- 1 / D
  Gamma <- diag(p) - outer(rep(1, p), omega)

  # apply transformation and mask
  out <- (L %*% t(Gamma)) * NA_U
  return(out)
}

#' Generate covariance matrix based on phylogenetic distance
#'
#' @param rho A scalar controlling the decay rate in the exponential kernel.
#' @param tau2 Variance component for the phylogenetic signal.
#' @param sigma2 Variance component for the residual noise.
#' @param D_matrix A pairwise phylogenetic distance matrix.
#' @param Gamma_c A projection matrix to adjust for selected taxa.
#' @param omega A vector of weights for taxa.
#'
#' @return A p x p covariance matrix.
#' @noRd
generate_covariance_matrix <- function(rho, tau2, sigma2, D_matrix, Gamma_c) {
  G_rho <- exp(-rho * D_matrix)
  part1 <- tau2 * Gamma_c %*% G_rho %*% t(Gamma_c)
  part2 <- sigma2 * Gamma_c %*% t(Gamma_c)
  Sigma_c <- part1 + part2
  return(Sigma_c)
}

#' Solve constrained least squares problem
#'
#' @param H The Hessian matrix.
#' @param b The gradient vector.
#' @param C The constraint matrix.
#'
#' @return A vector of parameter estimates satisfying the constraints.
#' @noRd
solve_constrained_theta <- function(H, b, C) {
  aug_matrix <- rbind(
    cbind(H, t(C)),
    cbind(C, matrix(0, nrow = nrow(C), ncol = nrow(C)))
  )
  rhs <- rbind(b, matrix(0, nrow = nrow(C), ncol = 1))
  sol <- qr.solve(aug_matrix, rhs, tol = 1e-20)
  theta_hat <- sol[1:nrow(H), , drop = FALSE]
  return(theta_hat)
}

#' Compute the full Hessian matrix for model parameters
#'
#' @param alpha_hat Estimated alpha vector.
#' @param theta_hat Estimated constrained coefficients.
#' @param sigma_params A numeric vector of estimated sigma parameters.
#' @param X Design matrix.
#' @param Z Response matrix.
#' @param tree Phylogenetic tree.
#' @param selected_taxa Indices of selected taxa.
#' @param nonzero_index List of nonzero indices per sample.
#' @param PhiiXi_list List of design matrices for each sample.
#'
#' @return The full Hessian matrix.
#' @noRd
compute_full_hessian <- function(alpha_hat, theta_hat, sigma_params,
                                 X, Z, tree, selected_taxa,
                                 nonzero_index, PhiiXi_list) {
  n <- nrow(X)
  p <- ncol(Z)
  q <- ncol(X)

  D_matrix <- cophenetic(tree)[colnames(Z), colnames(Z)]
  weights <- rep(0, p)
  weights[selected_taxa] <- 100^2 / (1 + 100^2 * length(selected_taxa))
  Gamma_c <- diag(1, p) - outer(rep(1, p), weights)

  Sigma_c <- generate_covariance_matrix(
    sigma_params[1], sigma_params[2], sigma_params[3], D_matrix, Gamma_c
  )
  Sigma_list <- lapply(1:n, function(i) Sigma_c[nonzero_index[[i]], nonzero_index[[i]]])

  H_total <- matrix(0, nrow = p * q + 3, ncol = p * q + 3)
  H_theta_theta <- fast_H_theta_theta(PhiiXi_list, Sigma_list)
  H_total[1:(p * q), 1:(p * q)] <- -H_theta_theta

  PhiiZi_list <- lapply(1:n, function(i) Z[i, nonzero_index[[i]]])

  theta_gradient <- function(rho, tau2, sigma2) {
    Sigma_c <- generate_covariance_matrix(rho, tau2, sigma2, D_matrix, Gamma_c)
    Sigma_list <- lapply(1:n, function(i) Sigma_c[nonzero_index[[i]], nonzero_index[[i]]])
    fast_theta_gradient(PhiiXi_list, PhiiZi_list, Sigma_list, theta_hat)
  }

  grad_theta_sigma <- numDeriv::jacobian(
    func = function(sigvec) theta_gradient(sigvec[1], sigvec[2], sigvec[3]),
    x = sigma_params
  )

  H_total[1:(p * q), (p * q + 1):(p * q + 3)] <- grad_theta_sigma
  H_total[(p * q + 1):(p * q + 3), 1:(p * q)] <- t(grad_theta_sigma)

  sigma_obj <- function(sigparams) {
    alpha_matrix <- matrix(alpha_hat, nrow = q, ncol = p, byrow = TRUE)
    Sigma_c <- generate_covariance_matrix(sigparams[1], sigparams[2], sigparams[3], D_matrix, Gamma_c)
    fast_objective_c(Z, X, alpha_matrix, nonzero_index, Sigma_c)
  }

  sigma_hess <- numDeriv::hessian(sigma_obj, sigma_params)
  H_total[(p * q + 1):(p * q + 3), (p * q + 1):(p * q + 3)] <- sigma_hess

  return(H_total)
}
