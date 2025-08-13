#' Differential-abundance pipeline for microbiome data (multiple covariates)
#'
#' Runs the McDisCoV analysis on a count table with covariates, performing
#' zero handling, GCLR transform, model fitting, Wald tests, and grouping.
#' Factors in `meta_data` are automatically expanded to 0/1 dummy variables.
#' An intercept term is appended internally as the last column of the design matrix.
#'
#' @param taxa_data matrix or data.frame. Observed counts with
#'   **rows = samples** and **columns = taxa**.
#' @param meta_data data.frame. Covariates with the **same rows (order) as `taxa_data`**.
#' @param variables character vector. Names of covariates in `meta_data` to include in
#'   the model (all will be tested while adjusting for the others). Categorical variables
#'   may be factors or characters; they will be dummy-encoded (0/1). Continuous variables
#'   should be numeric.
#' @param imputation logical. If `TRUE`, impute zeros for a selected set of taxa
#'   (see `zero_ratio_threshold` and `zero_impute`). If `FALSE`, a reference taxon is
#'   chosen; samples with zero for that reference may be dropped to avoid degenerate design.
#' @param selected_taxa integer vector or `NULL`. Indices of taxa to impute (if
#'   `imputation = TRUE`). If `NULL`, taxa are chosen by `zero_ratio_threshold`.
#' @param zero_ratio_threshold numeric in [0,1]. When `imputation = TRUE`, taxa with
#'   zero proportion `<` this threshold are selected for imputation.
#' @param zero_impute numeric. Value used to replace zeros during imputation.
#' @param Wc numeric in (0,1). Threshold used when clustering pairwise effect differences
#'   to define weak/strong association groups.
#' @param max_iter integer. Maximum iterations for the core optimizer (passed to
#'   internal routines if applicable).
#' @param tol numeric. Convergence tolerance for the core optimizer (passed through).
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Builds a design matrix from `meta_data[ , variables]`, dummy-encodes factors,
#'         and appends an intercept column.
#'   \item Handles zeros via imputation or reference-taxon selection.
#'   \item Applies a GCLR transform and fits the core model.
#'   \item Computes Wald tests per taxon and per variable and groups taxa by pairwise
#'         effect differences into \emph{not} / \emph{weakly} / \emph{strongly} associated.
#' }
#'
#' @return A list with components:
#' \describe{
#'   \item{AA}{A named list (one element per tested variable). Each element is a
#'     `data.frame` with columns:
#'       \describe{
#'         \item{Taxa}{Taxon name.}
#'         \item{effect_size}{Estimated effect size relative to the Group-A mean.}
#'         \item{DA_group}{Association strength per taxon: \code{0}=not, \code{1}=weak, \code{2}=strong.}
#'         \item{W}{Total number of significant pairwise differences for that taxon.}
#'       }}
#'   \item{RA}{A named list (one element per tested variable). Each element is a
#'     `data.frame` with columns:
#'       \describe{
#'         \item{Taxa}{Taxon name.}
#'         \item{Estimate}{Per-taxon coefficient estimate.}
#'         \item{StdError}{Standard error of the estimate.}
#'         \item{TestStatistic}{Wald test statistic (chi-squared, 1 df).}
#'         \item{p_val}{Wald p-value.}
#'         \item{DA}{Differential-abundance flag: \code{0}=not associated, \code{1}=associated.}
#'       }}
#'   \item{selected_taxa}{Indices of taxa retained/used after zero handling.}
#'   \item{numerical_res}{A list with numerical outputs:
#'     \code{estimated_theta}, \code{estimated_sigma_c}, \code{theta_cov_matrix}.}
#' }
#'
#'
#' @useDynLib McDisCoV, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @export
McDisCov <- function(
    taxa_data,                      # samples x taxa counts (matrix or data.frame)
    meta_data,                      # samples x covariates (data.frame; same rows as taxa_data)
    variables,                      # character vector of covariate names to include & test
    imputation = FALSE,
    selected_taxa = NULL,
    zero_ratio_threshold = 0.2,
    zero_impute = 0.5,
    Wc = 0.7
){
  ## --------------------------- 1) Data checks & design matrix ---------------------------
  message("Data preprocessing...")

  # coerce Z (counts) to numeric matrix
  Z <- as.matrix(taxa_data)
  mode(Z) <- "numeric"

  # row alignment
  if (nrow(Z) != nrow(meta_data)) {
    stop("Rows of taxa_data and meta_data must match (same samples, same order).")
  }

  # check variables
  if (!all(variables %in% colnames(meta_data))) {
    miss <- setdiff(variables, colnames(meta_data))
    stop("variables not found in meta_data: ", paste(miss, collapse = ", "))
  }

  message("Building design matrix X from meta_data[variables] ...")
  X_raw <- meta_data[, variables, drop = FALSE]

  # Expand factors to 0/1 dummies (no intercept); keep numeric as is
  if (any(!vapply(X_raw, is.numeric, logical(1)))) {
    X_vars <- model.matrix(~ 0 + ., data = X_raw)          # dummy expansion, no intercept
  } else {
    X_vars <- as.matrix(X_raw)
  }
  storage.mode(X_vars) <- "numeric"

  # Add intercept as the LAST column (your code expects that)
  X <- cbind(X_vars, Intercept = 1)

  # indices of the variables we want to TEST (exclude the intercept)
  which_X_index <- seq_len(ncol(X_vars))                    # test all provided variables

  ## --------------------------- 2) Zero handling / reference selection -------------------
  message("Selecting reference / zero handling ...")
  zero_ratio <- colSums(Z == 0) / nrow(Z)

  if (!is.null(selected_taxa)) {
    for (j in selected_taxa) Z[, j] <- ifelse(Z[, j] == 0, zero_impute, Z[, j])
  } else {
    if (imputation) {
      selected_taxa <- which(zero_ratio < zero_ratio_threshold)
      if (length(selected_taxa) == 0) {
        stop("No taxon selected. Choose another zero_ratio_threshold.")
      }
      for (j in selected_taxa) Z[, j] <- ifelse(Z[, j] == 0, zero_impute, Z[, j])
    } else {
      if (length(which(zero_ratio == 0)) == 0) {
        # pick a taxon with fewest zeros, drop samples where it is zero,
        # and ensure no covariate column becomes all-0 or all-1
        sorted_indices <- order(zero_ratio)
        i <- 1
        repeat {
          least0_index <- sorted_indices[i]
          remove_idx <- which(Z[, least0_index] == 0)
          Z_new <- Z[-remove_idx, , drop = FALSE]
          X_new <- X[-remove_idx, , drop = FALSE]

          X_tmp <- X_new[, -ncol(X_new), drop = FALSE]  # drop intercept
          is_all_0_or_1 <- apply(X_tmp, 2, function(col) all(col == 0) || all(col == 1))

          if (any(is_all_0_or_1)) {
            i <- i + 1
            if (i > length(sorted_indices)) {
              stop("No suitable taxon; all possibilities make some X column all 0 or all 1.")
            }
          } else {
            break
          }
        }
        selected_taxa <- sorted_indices[i]
        Z <- Z_new
        X <- X_new
      } else {
        selected_taxa <- which(zero_ratio == 0)
      }
    }
  }

  ## --------------------------- 3) Transform & core analysis -----------------------------
  message("GCLR transforming ...")
  Z_trans <- GCLR_transform(Z, selected_taxa, zero_impute)

  analysis_result <- run_analysis_core(Z_trans, X, selected_taxa)

  ## --------------------------- 4) Wald tests & grouping --------------------------------
  p <- ncol(Z)   # #taxa
  q <- ncol(X)   # #covariates incl. intercept

  which_X_index <- setdiff(seq_len(q), q)   # by default: test all except intercept
  # If you already built which_X_index earlier, keep that instead.

  AA <- list()   # <-- new
  RA <- list()   # <-- new

  taxa_names <- colnames(Z)

  for (k in which_X_index) {
    idx_range <- seq(k, p * q, by = q)
    alpha_extraction_index <- (k*p - p + 1):(k*p)

    est_alpha_this <- analysis_result$estimated_alpha[alpha_extraction_index]
    cov_this       <- analysis_result$cov_matrix[idx_range, idx_range]
    std_error      <- sqrt(diag(cov_this))

    ## --- Wald stats (you already had these) ---
    wald_stat <- (est_alpha_this / std_error)^2
    pval      <- 1 - pchisq(wald_stat, df = 1)
    p_star    <- ifelse(pval < 0.05, 1, 0)

    ## --- Pairwise differences (you already had these) ---
    var_vec  <- diag(cov_this)
    D        <- outer(est_alpha_this, est_alpha_this, "-")
    Var_mat  <- outer(var_vec, var_vec, "+") - 2 * cov_this
    Wmat     <- D / sqrt(Var_mat)
    Pmat     <- 1 - pchisq(Wmat^2, df = 1)

    wald_test_matrix <- (Pmat < 0.05) * 1
    wald_test_matrix[is.na(wald_test_matrix)] <- 0

    row_sums <- rowSums(wald_test_matrix)   # total significant pairwise differences per taxon

    ## --- Your grouping logic (A/B/C) gives group_star = 0/1/2 ---
    hca_count       <- hclust(dist(row_sums), "complete")
    hca_count_index <- cutree(hca_count, k = 2)

    min_group_index <- hca_count_index[which.min(row_sums)]
    max_group_index <- hca_count_index[which.max(row_sums)]
    Group_A  <- which(hca_count_index == min_group_index)
    Group_BC <- which(hca_count_index == max_group_index)
    Group_B  <- Group_BC[which(row_sums[Group_BC] <  (p * Wc))]
    Group_C  <- Group_BC[which(row_sums[Group_BC] >= (p * Wc))]

    group_star <- rep(0, p)  # 0 = not, 1 = weak, 2 = strong
    group_star[Group_B] <- 1
    group_star[Group_C] <- 2

    ## --- Effect size relative to Group A mean (as in your code) ---
    effect_size <- est_alpha_this - mean(est_alpha_this[Group_A])

    var_name <- colnames(X)[k]

    ## ===================== PACK NEW OUTPUTS =====================

    # 1) AA: per-variable summary for association strength
    #    columns: Taxa, effect_size, DA_group (0/1/2), W (row_sums)
    AA[[var_name]] <- data.frame(
      Taxa        = taxa_names,
      effect_size = effect_size,
      DA_group    = as.integer(group_star),
      W           = as.numeric(row_sums),
      row.names   = NULL
    )

    # 2) RA: per-variable regression table
    #    columns: Taxa, Estimate, StdError, TestStatistic, p_val, DA (0/1)
    RA[[var_name]] <- data.frame(
      Taxa          = taxa_names,
      Estimate      = est_alpha_this,
      StdError      = std_error,
      TestStatistic = wald_stat,          # (χ^2 with 1 df). If you prefer Z, use est/std instead.
      p_val         = pval,
      DA            = as.integer(p_star), # 0/1
      row.names     = NULL
    )
  }

  ## --------------------------- 5) Numerical results bundle ------------------------------
  # Your core currently exposes estimated_alpha and cov_matrix.
  # If it also provides theta/sigma_c/theta_cov_matrix, pass those through.
  estimated_theta    <- if ("estimated_theta"    %in% names(analysis_result)) analysis_result$estimated_theta    else analysis_result$estimated_alpha
  estimated_sigma_c  <- if ("estimated_sigma_c"  %in% names(analysis_result)) analysis_result$estimated_sigma_c  else NA
  theta_cov_matrix   <- if ("theta_cov_matrix"   %in% names(analysis_result)) analysis_result$theta_cov_matrix   else analysis_result$cov_matrix

  res <- list(
    AA             = AA,                         # per-variable: effect_size, DA_group, W
    RA             = RA,                         # per-variable regression table
    selected_taxa  = selected_taxa,
    numerical_res  = list(
      estimated_theta   = estimated_theta,
      estimated_sigma_c = estimated_sigma_c,
      theta_cov_matrix  = theta_cov_matrix
    )
  )

  return(res)
}

# McDisCov <- function(
#     Z, X,
#     imputation = FALSE,
#     selected_taxa = NULL,
#     zero_ratio_threshold = 0.2,
#     zero_impute = 0.5,
#     Wc = 0.7,
#     which_X = NULL)
# {
#   #-------------------------------------------------------------------------------
#   # 1. Data preprocessing
#   print("Data preprocessing...")
#   print("Selecting reference...")
#
#   if (!all(which_X %in% colnames(X))) {
#     stop("Elements of which_X must be column names of X!")
#   }
#   X <- cbind(X[,which_X,drop=FALSE], Intercept = 1)
#
#   zero_ratio <- colSums(Z == 0) / nrow(Z)
#   if (!is.null(selected_taxa)) {
#     for (j in selected_taxa) {
#       Z[, j] <- ifelse(Z[, j] == 0, zero_impute, Z[, j])
#     }
#   } else {
#     if (imputation) {
#       selected_taxa <- which(zero_ratio < zero_ratio_threshold)
#       if (length(selected_taxa)==0) {
#         stop("No taxon is selected. Choose another zero ratio threshold!")
#       }
#       for (j in selected_taxa) {
#         Z[, j] <- ifelse(Z[, j] == 0, zero_impute, Z[, j])
#       }
#     } else {
#       if (length(which(zero_ratio == 0)) == 0) {
#         # 假设 Z 和 X 已经定义好了
#         zero_ratio <- colSums(Z == 0) / nrow(Z)
#         sorted_indices <- order(zero_ratio)
#
#         i <- 1
#         repeat {
#           print("CHOOSE")
#           # 当前考虑的微生物是第 sorted_indices[i] 个
#           least0_index <- sorted_indices[i]
#
#           # 找到在这个微生物上为 0 的样本索引
#           remove_sample_index <- which(Z[, least0_index] == 0)
#
#           # 删除这些样本
#           Z_new <- Z[-remove_sample_index, , drop = FALSE]
#           X_new <- X[-remove_sample_index, , drop = FALSE]
#
#           # 去掉 X 的最后一列
#           X_tmp <- X_new[, -ncol(X_new), drop = FALSE]
#
#           # 检查是否有某一列全为 0 或全为 1
#           is_all_0_or_1 <- apply(X_tmp, 2, function(col) all(col == 0) || all(col == 1))
#
#           if (any(is_all_0_or_1)) {
#             # 如果有问题，继续检查下一个 zero_ratio 排名
#             i <- i + 1
#             if (i > length(sorted_indices)) {
#               stop("找不到满足条件的微生物，所有处理结果都会导致 X 的协变量列全为0或1。")
#             }
#           } else {
#             # 满足条件，停止
#             break
#           }
#         }
#
#         # 成功选择的微生物 index 和处理后的 Z 与 X
#         selected_taxa <- sorted_indices[i]
#         Z <- Z_new
#         X <- X_new
#
#       } else {
#         selected_taxa <- which(zero_ratio == 0)
#       }
#     }
#   }
#
#   print("GCLR transforming...")
#
#   Z_trans <- GCLR_transform(Z, selected_taxa, zero_impute)
#
#   #-------------------------------------------------------------------------------
#   # 2. Core analysis: model fitting and parameter estimation
#
#   analysis_result <- run_analysis_core(Z_trans, X, selected_taxa)
#   #-------------------------------------------------------------------------------
#   # 3. Wald test and grouping for all specified X columns
#   p <- ncol(Z)
#   q <- ncol(X)
#   which_X_index <- match(which_X, colnames(X))
#
#   wald_results_list <- list()
#   group_results_list <- list()
#
#   for (k in which_X_index) {
#     idx_range <- seq(k, p * q, by = q)
#     alpha_extraction_index <- (k*p-p+1):(k*p)
#     est_alpha_this <- analysis_result$estimated_alpha[alpha_extraction_index]
#     cov_this <- analysis_result$cov_matrix[idx_range, idx_range]
#     std_error <- sqrt(diag(cov_this))
#     # Wald statistics and p-values
#     wald_stat <- (est_alpha_this / std_error)^2
#     pval <- 1 - pchisq(wald_stat, df=1)
#     p_star <- ifelse(pval < 0.05, 1, 0)
#     # Grouping by pairwise difference
#     wald_test_matrix <- matrix(0, nrow=p, ncol=p)
#     for (i in 1:p) for (j in 1:p) {
#       if (i != j) {
#         var_i <- cov_this[i, i]
#         var_j <- cov_this[j, j]
#         cov_ij <- cov_this[i, j]
#         wald_diff <- (est_alpha_this[i] - est_alpha_this[j]) / sqrt(var_i + var_j - 2 * cov_ij)
#         p_diff <- 1 - pchisq(wald_diff^2, df=1)
#         if (is.na(p_diff)) {wald_test_matrix[i, j] <- 0} else {
#           if (p_diff < 0.05) {wald_test_matrix[i, j] <- 1}
#         }
#       }
#     }
#     # 1) extract variances
#     var_vec <- diag(cov_this)               # length-p
#
#     D <- outer(est_alpha_this, est_alpha_this, "-")  # p×p
#     Var_mat <- outer(var_vec, var_vec, "+") - 2 * cov_this
#
#     # 4) compute Wald Z‐scores
#     W <- D / sqrt(Var_mat)
#     P <- 1 - pchisq(W^2, df = 1)
#
#     # 6) build your indicator matrix (0/1), setting NAs to 0
#     wald_test_matrix <- (P < 0.05) * 1
#     wald_test_matrix[is.na(wald_test_matrix)] <- 0
#     row_sums <- rowSums(wald_test_matrix)
#     hca_count <- hclust(dist(row_sums), "complete")
#     hca_count_index <- cutree(hca_count, k=2)
#     min_group_index <- hca_count_index[which.min(row_sums)]
#     max_group_index <- hca_count_index[which.max(row_sums)]
#     Group_A <- which(hca_count_index == min_group_index)
#     Group_BC <- which(hca_count_index == max_group_index)
#     Group_B <- Group_BC[which(row_sums[Group_BC] < (p*Wc))]
#     Group_C <- Group_BC[which(row_sums[Group_BC] >= (p*Wc))]
#     group_star <- rep(0, p)
#     group_star[Group_C] <- 2
#     group_star[Group_B] <- 1
#     # Effect Size
#     effect_size <- est_alpha_this - mean(est_alpha_this[Group_A])
#
#     wald_results_list[[colnames(X)[k]]] <- data.frame(
#       Response = colnames(Z),
#       Estimate = est_alpha_this,
#       StdError = std_error,
#       WaldStatistic = wald_stat,
#       p_value = pval,
#       star = p_star
#     )
#     group_results_list[[colnames(X)[k]]] <- list(
#       effect_size = effect_size,
#       p_star = group_star,
#       Group_A = Group_A,
#       Group_B = Group_B,
#       Group_C = Group_C,
#       row_sums = row_sums
#     )
#   }
#   print("Testing done.")
#   return(list(
#     wald_results = wald_results_list,
#     group_results = group_results_list,
#     selected_taxa = selected_taxa,
#     analysis_result = analysis_result
#   ))
# }
