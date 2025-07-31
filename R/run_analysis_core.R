#’ Core alternating‐optimization routine for McDisCoV
#’
#’ Given a response matrix Z (n×p) with missing entries and a design matrix X (n×q),
#’ this function builds block‐diagonal design submatrices, estimates residuals via
#’ separate linear models, shrinks the empirical covariance toward its diagonal,
#’ selects the best shrinkage parameter by maximizing a log‐likelihood, and then
#’ solves for constrained parameters (using fast C++ backends).  Finally it computes
#’ the asymptotic covariance for the estimated coefficients.
#’
#’ @param Z A numeric matrix/data.frame (n × p), possibly containing NA’s.
#’ @param X A numeric matrix/data.frame (n × q) of covariates.
#’ @param selected_taxa Integer indices (or logical vector) of length p indicating
#’   which columns of Z are subject to the linear constraint.
#’ @return A list with components:
#’   \describe{
#’     \item{estimated_alpha}{Numeric vector (length p·q) of fitted coefficients.}
#’     \item{estimated_sigma_c}{Estimated shrinkage covariance matrix (p×p).}
#’     \item{cov_matrix}{Asymptotic covariance matrix of the alpha estimates.}
#’   }
#’ @details
#’  Internally this calls two C++ routines:
#’  \itemize{
#’    \item \code{fast_half1_half2(PhiiXi_list, PhiiZi_list, Sigma_list)} to form and solve
#’      the normal‐equations blocks,
#’    \item \code{fast_H_theta_theta(PhiiXi_list, Sigma_list)} to build the Hessian.
#’  }
#’
#’ @useDynLib McDisCoV, .registration=TRUE
#’ @importFrom Matrix nearPD
#’ @importFrom stats lm residuals pchisq kappa
#’ @export
run_analysis_core <- function(Z, X, selected_taxa) {
  X <- as.matrix(X)
  n <- nrow(X)
  q <- ncol(X)
  p <- ncol(Z)

  nonzero_index <- lapply(1:n, function(i) which(!is.na(Z[i, ])))

  PhiiXi_list <- lapply(1:n, function(i) {
    row_vector <- X[i, , drop = FALSE]
    identity_p <- diag(p)
    kronecker(identity_p, row_vector)[nonzero_index[[i]], , drop = FALSE]
  })

  # initial estimation
  get_residual_matrix <- function(Z, X) {
    n <- nrow(Z)
    p <- ncol(Z)
    resid_mat <- matrix(NA, nrow = n, ncol = p)

    for (j in 1:p) {
      y_j <- Z[, j]
      non_missing <- !is.na(y_j)

      # 仅用非缺失的位置拟合 lm
      fit <- lm(y_j[non_missing] ~ X[non_missing, , drop = FALSE])
      resid_mat[non_missing, j] <- residuals(fit)
    }

    return(resid_mat)
  }

  # Sigma估计目标函数
  log_lik <- function(Sigma_c) {
    objective_value <- 0

    for (i in 1:n) {
      nonzero_i <- nonzero_index[[i]]
      PhiiZi <- Z[i, nonzero_i]
      PhiiSigma_cPhiiT <- Sigma_c[nonzero_i, nonzero_i]

      Ri <- chol(PhiiSigma_cPhiiT)

      log_det_i <- 2 * sum(log(diag(Ri)))
      v <- forwardsolve(t(Ri), PhiiZi)
      quad <- sum(v^2)

      # Update total objective
      objective_value <- objective_value + (-0.5) * log_det_i + (-0.5) * quad
    }
    objective_value
  }

  # Sigma估计
  generate_covariance_matrix <- function(lambda_seq = seq(0, 1, by = 0.1)) {
    delta_mat <- ifelse(is.na(Z), 0, 1)
    residual <- get_residual_matrix(Z, X)

    # 构造残差协方差估计
    Sigma_emp <- matrix(0, p, p)
    count_mat <- matrix(0, p, p)

    for (i in 1:n) {
      delta_i <- delta_mat[i, ]
      r_i <- residual[i, ]
      if (sum(delta_i) == 0) next

      r_i[delta_i == 0] <- 0
      Sigma_emp <- Sigma_emp + tcrossprod(r_i)
      count_mat <- count_mat + tcrossprod(delta_i)
    }

    Sigma_emp <- ifelse(count_mat > 0, Sigma_emp / count_mat, 0)

    # 正定性判断函数
    is_posdef <- function(Sigma, tol = 1e-8) {
      if (!isSymmetric(Sigma)) return(FALSE)
      eigvals <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
      all(eigvals > tol)
    }

    best_loglik <- -Inf
    best_Sigma <- NULL
    loglik_values <- numeric(length(lambda_seq))

    for (i in seq_along(lambda_seq)) {
      lambda <- lambda_seq[i]

      # shrink to diagonal
      Tm <- diag(diag(Sigma_emp))
      Sigma_c <- (1 - lambda) * Sigma_emp + lambda * Tm

      # 若不正定强行正定
      Sigma_c <- if (is_posdef(Sigma_c)) Sigma_c else as.matrix(Matrix::nearPD(Sigma_c)$mat)

      # 计算 log-likelihood损失函数
      ll <- log_lik(Sigma_c)
      loglik_values[i] <- ll

      if (ll > best_loglik) {
        best_loglik <- ll
        best_Sigma <- Sigma_c
      }
    }
    return(list(
      best_lambda = lambda_seq[which.max(loglik_values)],
      best_loglik = best_loglik,
      best_cov = best_Sigma,
      all_loglik = loglik_values,
      lambda_seq = lambda_seq
    ))
  }

  # theta约束方程求解构造
  s <- rep(0, p)
  s[selected_taxa] <- 1
  q <- ncol(X)
  C <- kronecker(t(s), diag(q))
  solve_constrained_theta <- function(H, b, C) {
    aug_matrix <- rbind(
      cbind(H, t(C)),
      cbind(C, matrix(0, nrow = q, ncol = q))
    )
    rhs <- rbind(b, matrix(0, nrow = q, ncol = 1))
    sol <- qr.solve(aug_matrix, rhs, tol = 1e-20)
    theta_hat <- sol[1:nrow(H), , drop = FALSE]
    return(theta_hat)
  }

  optimize_alternating <- function(max_iter, tol) {
    print("sigma......")
    Sigma_c <- generate_covariance_matrix()$best_cov
    Sigma_list <- lapply(1:n, function(i) Sigma_c[nonzero_index[[i]], nonzero_index[[i]]])

    print("alpha......")
    PhiiZi_list <- lapply(1:n, function(i) Z[i, nonzero_index[[i]]])
    result <- fast_half1_half2(PhiiXi_list, PhiiZi_list, Sigma_list)
    half1 <- result$half1
    half2 <- matrix(result$half2, ncol=1)
    theta_hat <- solve_constrained_theta(half1, half2, C)
    alpha_vector <- as.vector(t(matrix(theta_hat, nrow = q, ncol = p)))

    list(alpha_vector = alpha_vector, Sigma_c = Sigma_c)
  }

  result <- optimize_alternating(max_iter, tol)

  print("Calculating asymptotic covariance...")
  compute_full_hessian <- function() {
    alpha_hat <- result$alpha_vector
    Sigma_c_hat <- result$Sigma_c

    Sigma_list <- lapply(1:n, function(i) Sigma_c_hat[nonzero_index[[i]], nonzero_index[[i]]])

    # 2. Rcpp加速求累加
    H_theta_theta <- fast_H_theta_theta(PhiiXi_list, Sigma_list)
    H_total <- -H_theta_theta

    return(H_total)
  }
  H_full <- compute_full_hessian()
  # 如果病态则加扰动
  if (kappa(H_full) > 1e12) {
    H_full <- H_full + diag(1e-8, nrow(H_full))
    print("add")
  }
  JC <- t(C)
  inv_Hessian <- -qr.solve(H_full, tol = 1e-20)

  cov_matrix <- inv_Hessian - inv_Hessian %*% JC %*%
    qr.solve(t(JC) %*% inv_Hessian %*% JC, tol = 1e-20) %*% t(JC) %*% inv_Hessian
  # 如果数值原因小于零则替换为零
  diag(cov_matrix)[diag(cov_matrix) < 0] <- 0

  estimated_alpha <- result$alpha_vector
  estimated_sigma_c <- result$Sigma_c

  return(list(
    estimated_alpha = estimated_alpha,
    estimated_sigma_c = estimated_sigma_c,
    cov_matrix = cov_matrix
  ))
}
