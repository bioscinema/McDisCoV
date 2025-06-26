#' Core estimation function for differential abundance modeling
#'
#' This function performs constrained optimization with alternating updates of alpha and sigma parameters
#' under a compositional multivariate log-normal model.
#'
#' @param Z A matrix of microbial abundances (samples x taxa).
#' @param X A design matrix (samples x covariates), including intercept if needed.
#' @param selected_taxa A vector of indices indicating which taxa are selected.
#' @param tree A phylogenetic tree object (class "phylo").
#' @param max_iter Maximum number of iterations for alternating optimization.
#' @param tol Convergence tolerance.
#'
#' @return A list containing estimated alpha vector, sigma parameters, and covariance matrix.
#' @export
run_analysis_core <- function(Z, X, selected_taxa, tree, max_iter = 20, tol = 1e-5) {
  X <- as.matrix(X)
  n <- nrow(X)
  q <- ncol(X)
  p <- ncol(Z)

  init_params <- list(rho = 1, tau2 = 1, sigma2 = 1)
  nonzero_index <- lapply(1:n, function(i) which(!is.na(Z[i, ])))

  # Construct PhiiXi list
  PhiiXi_list <- lapply(1:n, function(i) {
    row_vector <- X[i, , drop = FALSE]
    kronecker(diag(p), row_vector)[nonzero_index[[i]], , drop = FALSE]
  })

  # Precompute constants
  c_val <- 100
  omega_val <- c_val^2 / (1 + c_val^2 * length(selected_taxa))
  weights <- rep(0, p)
  weights[selected_taxa] <- omega_val
  omega <- weights
  Gamma_c <- diag(1, p) - outer(rep(1, p), omega)
  D_matrix <- cophenetic(tree)[colnames(Z), colnames(Z)]

  # Constraint matrix C
  s <- rep(0, p)
  s[selected_taxa] <- 1
  C <- kronecker(t(s), diag(q))

  # Objective function
  objective_c <- function(params) {
    alpha_matrix <- matrix(params[1:(q * p)], nrow = q, byrow = TRUE)
    Sigma_c <- generate_covariance_matrix(
      rho = params[(q * p + 1)],
      tau2 = params[(q * p + 2)],
      sigma2 = params[(q * p + 3)],
      D_matrix, Gamma_c
    )
    fast_objective_c(Z, X, alpha_matrix, nonzero_index, Sigma_c)
  }

  # Alternating optimization
  optimize_alternating <- function(init_params, max_iter, tol) {
    sigma_params <- c(init_params$rho, init_params$tau2, init_params$sigma2)
    last_sigma <- sigma_params
    PhiiZi_list <- lapply(1:n, function(i) Z[i, nonzero_index[[i]]])

    Sigma_c <- generate_covariance_matrix(sigma_params[1], sigma_params[2], sigma_params[3], D_matrix, Gamma_c)
    Sigma_list <- lapply(1:n, function(i) Sigma_c[nonzero_index[[i]], nonzero_index[[i]]])
    result <- fast_half1_half2(PhiiXi_list, PhiiZi_list, Sigma_list)
    theta_hat <- solve_constrained_theta(result$half1, matrix(result$half2, ncol = 1), C)

    alpha_vector <- as.vector(t(matrix(theta_hat, nrow = q)))
    last_alpha <- alpha_vector
    last_objective <- objective_c(c(alpha_vector, sigma_params))

    for (iter in 1:max_iter) {
      sigma_opt <- nloptr::nloptr(
        x0 = sigma_params,
        eval_f = function(sig) {
          obj <- objective_c(c(alpha_vector, sig))
          if (is.infinite(obj) || is.nan(obj)) return(1e10)
          -obj
        },
        lb = c(0, 0, 0),
        ub = c(Inf, Inf, Inf),
        opts = list(algorithm = "NLOPT_LN_BOBYQA", xtol_rel = tol, maxeval = 500)
      )
      sigma_params <- sigma_opt$solution

      Sigma_c <- generate_covariance_matrix(sigma_params[1], sigma_params[2], sigma_params[3], D_matrix, Gamma_c)
      Sigma_list <- lapply(1:n, function(i) Sigma_c[nonzero_index[[i]], nonzero_index[[i]]])
      result <- fast_half1_half2(PhiiXi_list, PhiiZi_list, Sigma_list)
      theta_hat <- solve_constrained_theta(result$half1, matrix(result$half2, ncol = 1), C)

      alpha_vector <- as.vector(t(matrix(theta_hat, nrow = q)))
      current_objective <- objective_c(c(alpha_vector, sigma_params))
      rel_change <- abs(current_objective - last_objective) / max(abs(last_objective), 1e-10)

      if (rel_change < tol &&
          norm(alpha_vector - last_alpha, type = "2") < tol &&
          norm(Sigma_c - generate_covariance_matrix(last_sigma[1], last_sigma[2], last_sigma[3], D_matrix, Gamma_c), type = "F") < tol) {
        break
      }

      last_objective <- current_objective
      last_alpha <- alpha_vector
      last_sigma <- sigma_params
    }

    list(theta_hat = theta_hat, alpha_vector = alpha_vector,
         sigma_params = sigma_params, objective = current_objective)
  }

  result <- optimize_alternating(init_params, max_iter, tol)

  # Compute full Hessian
  H_full <- compute_full_hessian(
    result$alpha_vector, result$theta_hat, result$sigma_params,
    X, Z, tree, selected_taxa, nonzero_index, PhiiXi_list
  )

  H_full <- H_full + diag(1e-8, nrow(H_full))  # Stabilize
  JC <- t(cbind(C, 0, 0, 0))
  inv_H <- -qr.solve(H_full, tol = 1e-20)
  cov_matrix <- inv_H - inv_H %*% JC %*%
    qr.solve(t(JC) %*% inv_H %*% JC, tol = 1e-20) %*% t(JC) %*% inv_H

  list(
    estimated_alpha = result$alpha_vector,
    estimated_sigma_params = result$sigma_params,
    cov_matrix = cov_matrix
  )
}
