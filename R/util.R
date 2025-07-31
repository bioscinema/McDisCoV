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

#' Generalized Centered Log-Ratio (gCLR) Transformation
#'
#' Apply a generalized CLR transform to a count matrix by balancing a subset
#' of taxa. Zeros are first replaced by a small constant, then log-transformed,
#' and finally centered on the geometric mean of the selected taxa. Any taxa
#' that were zero in the original data (and not part of the balance) are set to NA.
#'
#' @param matrix1 Numeric matrix of counts (rows = samples, columns = taxa).
#' @param selected_taxa Integer vector of column indices specifying which taxa
#'   define the geometric mean (the “balance”).
#' @param zero_impute Numeric scalar value used to replace zero counts before
#'   the log transformation.
#'
#' @return A numeric matrix of the same dimensions as \code{matrix1}, containing
#'   the gCLR-transformed values. Entries corresponding to non-balanced taxa
#'   that were originally zero are returned as \code{NA}.
#'
#' @noRd
GCLR_transform <- function(matrix1, selected_taxa, zero_impute) {
  U <- ifelse(matrix1 == 0, 0, 1)
  U[, selected_taxa] <- 1
  NA_U <- ifelse(U == 0, NA, 1)
  p <- ncol(matrix1)
  D <- length(selected_taxa)
  matrix2 <- ifelse(matrix1 == 0, zero_impute, matrix1)
  L <- log(matrix2)
  weights <- t(rep(0, p))
  omega <- 100^2 / (1 + 100^2 * D) # 1/D
  weights[selected_taxa] <- omega
  Gamma <- diag(1, p) - rep(1, p) %*% weights
  matrix2 <- L %*% t(Gamma)
  matrix2 <- NA_U * matrix2
  return(matrix2)
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


#' Preprocess phyloseq object for differential abundance analysis
#'
#' @param physeq A phyloseq object containing count data and sample metadata
#' @param group_var Character; column name in sample_data defining groups
#' @param level Character or NULL; taxonomic rank to aggregate (e.g. "Family"); default NULL
#' @param imputation Logical; whether to perform zero-value imputation
#' @param zero_ratio_threshold Numeric; threshold for zero-rate above which taxa are imputed
#' @param zero_impute Numeric; value to add to zero counts when imputing
#' @param min_nonzero Integer; minimum number of non-zero samples per group to retain a taxon
#' @return A list with elements: otu (filtered OTU matrix), metadata, design (model matrix), tree, taxa_kept
#' @export
preprocess_inputs <- function(
    physeq,
    group_var,
    level                = NULL,
    imputation           = TRUE,
    zero_ratio_threshold = 0.1,
    zero_impute          = 0.5,
    min_nonzero          = 5
) {
  library(phyloseq)

  # 1) Optional aggregation at specified taxonomic rank
  if (!is.null(level)) {
    physeq <- tax_glom(physeq, taxrank = level)
  }

  # 2) Extract OTU count matrix and metadata
  otu_mat <- t(as(otu_table(physeq), "matrix"))
  meta_df <- as(sample_data(physeq), "data.frame")
  otu_mat <- otu_mat[rownames(meta_df), , drop = FALSE]

  # 3) Build dummy-coded design matrix for group_var
  group      <- factor(meta_df[[group_var]])
  design_mat <- model.matrix(~ 0 + group)
  colnames(design_mat) <- levels(group)
  rownames(design_mat) <- rownames(meta_df)

  # 4) Filter taxa by prevalence per group
  taxa_names <- colnames(otu_mat)
  keep_taxa  <- vapply(taxa_names, function(taxon) {
    nz_count <- tapply(otu_mat[, taxon] > 0, group, sum)
    all(nz_count >= min_nonzero)
  }, logical(1))
  otu_mat_f     <- otu_mat[, keep_taxa, drop = FALSE]
  filtered_taxa <- taxa_names[keep_taxa]

  # 5) Prune phyloseq object and extract its phylogenetic tree
  phy_filtered <- prune_taxa(filtered_taxa, physeq)
  tree         <- phy_tree(phy_filtered)

  # 6) Zero-value imputation if requested
  if (imputation) {
    zero_rate     <- colMeans(otu_mat_f == 0)
    to_impute_idx <- zero_rate > zero_ratio_threshold
    otu_mat_f[, to_impute_idx] <- otu_mat_f[, to_impute_idx] + zero_impute
  }

  # Return preprocessed inputs
  list(
    otu        = otu_mat_f,
    metadata   = meta_df,
    design     = design_mat,
    tree       = tree,
    taxa_kept  = filtered_taxa
  )
}
