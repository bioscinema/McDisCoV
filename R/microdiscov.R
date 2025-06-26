#' Full differential‐abundance pipeline with a single grouping variable
#'
#' @param physeq      A phyloseq object with OTU table, sample_data, and phy_tree
#' @param group_var   Character; name of a sample_data column (numeric or factor)
#' @param imputation  Logical; if TRUE, impute zeros, else drop taxa with any zeros
#' @param zero_ratio_threshold Numeric; maximum zero ratio for imputation
#' @param zero_impute Numeric; value to replace zeros when imputing
#' @param Wc          Numeric in (0,1); threshold for clustering effect strengths
#' @param max_iter    Integer; max iterations for the core optimizer
#' @param tol         Numeric; convergence tolerance
#' @param min_nonzero Integer; minimum non-zero counts per group (only for factor group)
#' @return A list with elements:
#'   * wald_results: data.frame of per-taxon Wald tests
#'   * group_results: list of taxon groupings by effect size
#'   * selected_taxa: integer vector of taxa kept after zero-handling
#'   * core_fit: raw output from run_analysis_core
#' @useDynLib McDisCoV, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @export
differential_abundance <- function(physeq,
                                   group_var,
                                   imputation           = TRUE,
                                   zero_ratio_threshold = 0.1,
                                   zero_impute          = 0.5,
                                   Wc                   = 0.7,
                                   max_iter             = 30,
                                   tol                  = 1e-10,
                                   min_nonzero          = 5) {
  # dependencies
  requireNamespace("ape")
  requireNamespace("nloptr")

  # 1) Extract count matrix Z
  Z <- as(otu_table(physeq), "matrix")
  if (taxa_are_rows(physeq)) Z <- t(Z)

  # 2) Extract grouping variable and build design matrix X
  sd <- data.frame(sample_data(physeq))
  if (!group_var %in% colnames(sd))
    stop("`group_var` not found in sample_data")
  group <- sd[[group_var]]
  X <- model.matrix(~ group, data = sd)

  # 3) Optionally filter low-prevalence taxa for two-level factors
  if (is.factor(group) && nlevels(group) == 2) {
    gb <- as.integer(group) - 1L  # convert to 0/1
    Z <- filter_by_prevalence(Z, gb, min_nonzero)
  }

  # 4) Extract phylogenetic tree
  tree <- phy_tree(physeq)

  # 5) Zero-handling
  zres         <- zero_filter_impute(Z, imputation,
                                     zero_ratio_threshold,
                                     zero_impute)
  Z_filt       <- zres$Z
  selected_taxa <- zres$selected_taxa

  # 6) GCLR transform
  Z_trans <- gclr_transform(Z_filt, selected_taxa, zero_impute)

  # 7) Core model fitting
  core_fit <- run_analysis_core(Z_trans, X, selected_taxa,
                                tree, max_iter, tol)

  # 8) Wald tests on the “group” coefficient (column 2 of X)
  p         <- ncol(Z)
  q         <- ncol(X)
  est_alpha <- core_fit$estimated_alpha[(q + 1):(2 * q)]
  cov_mat   <- core_fit$cov_matrix[(q + 1):(2 * q),
                                   (q + 1):(2 * q)]
  se        <- sqrt(diag(cov_mat))
  wald_stat <- (est_alpha / se)^2
  p_vals    <- 1 - pchisq(wald_stat, df = 1)
  wald_results <- data.frame(
    taxon      = colnames(Z),
    estimate   = est_alpha,
    std_error  = se,
    wald_stat  = wald_stat,
    p_value    = p_vals,
    significant = p_vals < 0.05
  )

  # 9) Cluster taxa by effect‐size differences
  group_results <- cluster_taxa_by_wald(est_alpha, cov_mat, Wc)

  # 10) Return all results
  list(
    wald_results   = wald_results,
    group_results  = group_results,
    selected_taxa  = selected_taxa,
    core_fit       = core_fit
  )
}
