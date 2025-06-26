// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
using namespace Rcpp;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// [[Rcpp::export]]
double fast_objective_c(const Eigen::MatrixXd& Z,
                        const Eigen::MatrixXd& X,
                        const Eigen::MatrixXd& alpha_matrix,
                        const List& nonzero_index,
                        const Eigen::MatrixXd& Sigma_c) {
  int n = Z.rows();
  int p = Z.cols();
  int q = X.cols();
  double obj = 0.0;

  Eigen::MatrixXd mu_matrix = X * alpha_matrix;

  for (int i = 0; i < n; ++i) {
    IntegerVector nz = nonzero_index[i];
    int m = nz.size();
    if (m == 0) continue;
    VectorXd PhiiZi(m), Phiimui(m);
    MatrixXd PhiiSigma_cPhiiT(m, m);
    for (int j = 0; j < m; ++j) {
      PhiiZi(j) = Z(i, nz[j] - 1);
      Phiimui(j) = mu_matrix(i, nz[j] - 1);
      for (int k = 0; k < m; ++k) {
        PhiiSigma_cPhiiT(j, k) = Sigma_c(nz[j] - 1, nz[k] - 1);
      }
    }

    Eigen::LLT<Eigen::MatrixXd> lltOfA(PhiiSigma_cPhiiT);
    if (lltOfA.info() != Eigen::Success) {
      ::Rf_warning("Cholesky decomposition failed for sample %d", i+1);
      continue;
    }
    Eigen::MatrixXd Ri = lltOfA.matrixU();
    double log_det_i = 2 * Ri.diagonal().array().log().sum();

    VectorXd mu_difference = PhiiZi - Phiimui;
    VectorXd v = Ri.transpose().triangularView<Eigen::Lower>().solve(mu_difference);

    double quad = v.squaredNorm();

    obj += (-0.5) * log_det_i + (-0.5) * quad;
  }
  return obj;
}

// [[Rcpp::export]]
List fast_half1_half2(const List& PhiiXi_list,
                      const List& PhiiZi_list,
                      const List& Sigma_list) {
  int n = PhiiXi_list.size();
  int pq = as<MatrixXd>(PhiiXi_list[0]).cols();
  MatrixXd half1 = MatrixXd::Zero(pq, pq);
  VectorXd half2 = VectorXd::Zero(pq);

  for (int i = 0; i < n; ++i) {
    MatrixXd PhiiXi = as<MatrixXd>(PhiiXi_list[i]);
    VectorXd PhiiZi = as<VectorXd>(PhiiZi_list[i]);
    MatrixXd M = as<MatrixXd>(Sigma_list[i]);
    if(M.rows() == 0 || M.cols() == 0) {
      stop("There is at least one pair of intercept and variable linearly dependent after zero separation!");
    }
    Eigen::LLT<Eigen::MatrixXd> lltOfA(M);
    if(lltOfA.info() != Eigen::Success) {
      ::Rf_warning("Cholesky decomposition failed at sample %d", i+1);
      continue;
    }
    MatrixXd inv_Mi = lltOfA.solve(MatrixXd::Identity(M.rows(), M.cols()));
    MatrixXd tPhiiXi_invMi = PhiiXi.transpose() * inv_Mi;
    half1.noalias() += tPhiiXi_invMi * PhiiXi;
    half2.noalias() += tPhiiXi_invMi * PhiiZi;
  }
  return List::create(Named("half1") = half1, Named("half2") = half2);
}

// [[Rcpp::export]]
Eigen::MatrixXd fast_H_theta_theta(const List& PhiiXi_list,
                                   const List& Sigma_list) {
  int n = PhiiXi_list.size();
  int pq = as<Eigen::MatrixXd>(PhiiXi_list[0]).cols();
  Eigen::MatrixXd H_theta_theta = Eigen::MatrixXd::Zero(pq, pq);

  for (int i = 0; i < n; ++i) {
    Eigen::MatrixXd PhiiXi = as<Eigen::MatrixXd>(PhiiXi_list[i]);
    Eigen::MatrixXd M = as<Eigen::MatrixXd>(Sigma_list[i]);
    if (M.rows() == 0 || M.cols() == 0) {
      stop("Zero-dimension covariance in Hessian block!");
    }
    Eigen::LLT<Eigen::MatrixXd> lltOfA(M);
    if (lltOfA.info() != Eigen::Success) {
      ::Rf_warning("Cholesky decomposition failed at sample %d (Hessian)", i+1);
      continue;
    }
    Eigen::MatrixXd inv_Mi = lltOfA.solve(Eigen::MatrixXd::Identity(M.rows(), M.cols()));
    H_theta_theta.noalias() += PhiiXi.transpose() * inv_Mi * PhiiXi;
  }
  return H_theta_theta;
}

// [[Rcpp::export]]
Eigen::VectorXd fast_theta_gradient(const List& PhiiXi_list,
                                    const List& PhiiZi_list,
                                    const List& Sigma_list,
                                    const Eigen::VectorXd& theta_hat) {
  int n = PhiiXi_list.size();
  int pq = as<Eigen::MatrixXd>(PhiiXi_list[0]).cols();
  Eigen::VectorXd grad_theta = Eigen::VectorXd::Zero(pq);

  for (int i = 0; i < n; ++i) {
    Eigen::MatrixXd PhiiXi = as<Eigen::MatrixXd>(PhiiXi_list[i]);
    Eigen::VectorXd PhiiZi = as<Eigen::VectorXd>(PhiiZi_list[i]);
    Eigen::MatrixXd M = as<Eigen::MatrixXd>(Sigma_list[i]);
    if (M.rows() == 0 || M.cols() == 0) {
      stop("Zero-dimension covariance in theta-gradient block!");
    }
    Eigen::LLT<Eigen::MatrixXd> lltOfA(M);
    if (lltOfA.info() != Eigen::Success) {
      ::Rf_warning("Cholesky decomposition failed at sample %d (grad)", i+1);
      continue;
    }
    Eigen::MatrixXd inv_Mi = lltOfA.solve(Eigen::MatrixXd::Identity(M.rows(), M.cols()));
    Eigen::VectorXd mu_i = PhiiXi * theta_hat;
    Eigen::VectorXd diff = PhiiZi - mu_i;
    grad_theta.noalias() += PhiiXi.transpose() * inv_Mi * diff;
  }
  return grad_theta;
}
