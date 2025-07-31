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
McDisCov <- function(
    Z, X,
    imputation = FALSE,
    selected_taxa = NULL,
    zero_ratio_threshold = 0.2,
    zero_impute = 0.5,
    Wc = 0.7,
    which_X = NULL)
{
  #-------------------------------------------------------------------------------
  # 1. Data preprocessing
  print("Data preprocessing...")
  print("Selecting reference...")

  if (!all(which_X %in% colnames(X))) {
    stop("Elements of which_X must be column names of X!")
  }
  X <- cbind(X[,which_X,drop=FALSE], Intercept = 1)

  zero_ratio <- colSums(Z == 0) / nrow(Z)
  if (!is.null(selected_taxa)) {
    for (j in selected_taxa) {
      Z[, j] <- ifelse(Z[, j] == 0, zero_impute, Z[, j])
    }
  } else {
    if (imputation) {
      selected_taxa <- which(zero_ratio < zero_ratio_threshold)
      if (length(selected_taxa)==0) {
        stop("No taxon is selected. Choose another zero ratio threshold!")
      }
      for (j in selected_taxa) {
        Z[, j] <- ifelse(Z[, j] == 0, zero_impute, Z[, j])
      }
    } else {
      if (length(which(zero_ratio == 0)) == 0) {
        # 假设 Z 和 X 已经定义好了
        zero_ratio <- colSums(Z == 0) / nrow(Z)
        sorted_indices <- order(zero_ratio)

        i <- 1
        repeat {
          print("CHOOSE")
          # 当前考虑的微生物是第 sorted_indices[i] 个
          least0_index <- sorted_indices[i]

          # 找到在这个微生物上为 0 的样本索引
          remove_sample_index <- which(Z[, least0_index] == 0)

          # 删除这些样本
          Z_new <- Z[-remove_sample_index, , drop = FALSE]
          X_new <- X[-remove_sample_index, , drop = FALSE]

          # 去掉 X 的最后一列
          X_tmp <- X_new[, -ncol(X_new), drop = FALSE]

          # 检查是否有某一列全为 0 或全为 1
          is_all_0_or_1 <- apply(X_tmp, 2, function(col) all(col == 0) || all(col == 1))

          if (any(is_all_0_or_1)) {
            # 如果有问题，继续检查下一个 zero_ratio 排名
            i <- i + 1
            if (i > length(sorted_indices)) {
              stop("找不到满足条件的微生物，所有处理结果都会导致 X 的协变量列全为0或1。")
            }
          } else {
            # 满足条件，停止
            break
          }
        }

        # 成功选择的微生物 index 和处理后的 Z 与 X
        selected_taxa <- sorted_indices[i]
        Z <- Z_new
        X <- X_new

      } else {
        selected_taxa <- which(zero_ratio == 0)
      }
    }
  }

  print("GCLR transforming...")

  Z_trans <- GCLR_transform(Z, selected_taxa, zero_impute)

  #-------------------------------------------------------------------------------
  # 2. Core analysis: model fitting and parameter estimation

  analysis_result <- run_analysis_core(Z_trans, X, selected_taxa)
  #-------------------------------------------------------------------------------
  # 3. Wald test and grouping for all specified X columns
  p <- ncol(Z)
  q <- ncol(X)
  which_X_index <- match(which_X, colnames(X))

  wald_results_list <- list()
  group_results_list <- list()

  for (k in which_X_index) {
    idx_range <- seq(k, p * q, by = q)
    alpha_extraction_index <- (k*p-p+1):(k*p)
    est_alpha_this <- analysis_result$estimated_alpha[alpha_extraction_index]
    cov_this <- analysis_result$cov_matrix[idx_range, idx_range]
    std_error <- sqrt(diag(cov_this))
    # Wald statistics and p-values
    wald_stat <- (est_alpha_this / std_error)^2
    pval <- 1 - pchisq(wald_stat, df=1)
    p_star <- ifelse(pval < 0.05, 1, 0)
    # Grouping by pairwise difference
    wald_test_matrix <- matrix(0, nrow=p, ncol=p)
    for (i in 1:p) for (j in 1:p) {
      if (i != j) {
        var_i <- cov_this[i, i]
        var_j <- cov_this[j, j]
        cov_ij <- cov_this[i, j]
        wald_diff <- (est_alpha_this[i] - est_alpha_this[j]) / sqrt(var_i + var_j - 2 * cov_ij)
        p_diff <- 1 - pchisq(wald_diff^2, df=1)
        if (is.na(p_diff)) {wald_test_matrix[i, j] <- 0} else {
          if (p_diff < 0.05) {wald_test_matrix[i, j] <- 1}
        }
      }
    }
    # 1) extract variances
    var_vec <- diag(cov_this)               # length-p

    D <- outer(est_alpha_this, est_alpha_this, "-")  # p×p
    Var_mat <- outer(var_vec, var_vec, "+") - 2 * cov_this

    # 4) compute Wald Z‐scores
    W <- D / sqrt(Var_mat)
    P <- 1 - pchisq(W^2, df = 1)

    # 6) build your indicator matrix (0/1), setting NAs to 0
    wald_test_matrix <- (P < 0.05) * 1
    wald_test_matrix[is.na(wald_test_matrix)] <- 0
    row_sums <- rowSums(wald_test_matrix)
    hca_count <- hclust(dist(row_sums), "complete")
    hca_count_index <- cutree(hca_count, k=2)
    min_group_index <- hca_count_index[which.min(row_sums)]
    max_group_index <- hca_count_index[which.max(row_sums)]
    Group_A <- which(hca_count_index == min_group_index)
    Group_BC <- which(hca_count_index == max_group_index)
    Group_B <- Group_BC[which(row_sums[Group_BC] < (p*Wc))]
    Group_C <- Group_BC[which(row_sums[Group_BC] >= (p*Wc))]
    group_star <- rep(0, p)
    group_star[Group_C] <- 2
    group_star[Group_B] <- 1
    # Effect Size
    effect_size <- est_alpha_this - mean(est_alpha_this[Group_A])

    wald_results_list[[colnames(X)[k]]] <- data.frame(
      Response = colnames(Z),
      Estimate = est_alpha_this,
      StdError = std_error,
      WaldStatistic = wald_stat,
      p_value = pval,
      star = p_star
    )
    group_results_list[[colnames(X)[k]]] <- list(
      effect_size = effect_size,
      p_star = group_star,
      Group_A = Group_A,
      Group_B = Group_B,
      Group_C = Group_C,
      row_sums = row_sums
    )
  }
  print("Testing done.")
  return(list(
    wald_results = wald_results_list,
    group_results = group_results_list,
    selected_taxa = selected_taxa,
    analysis_result = analysis_result
  ))
}
