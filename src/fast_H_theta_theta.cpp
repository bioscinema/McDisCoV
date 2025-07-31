// fast_H_theta_theta.cpp

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
using namespace Rcpp;

// [[Rcpp::export]]
Eigen::MatrixXd fast_H_theta_theta(const List& PhiiXi_list,
                                   const List& Sigma_list) {
  int n = PhiiXi_list.size();
  int pq = as<Eigen::MatrixXd>(PhiiXi_list[0]).cols();
  Eigen::MatrixXd H_theta_theta = Eigen::MatrixXd::Zero(pq, pq);

  for (int i = 0; i < n; ++i) {
    Eigen::MatrixXd PhiiXi = as<Eigen::MatrixXd>(PhiiXi_list[i]);
    Eigen::MatrixXd M     = as<Eigen::MatrixXd>(Sigma_list[i]);
    if (M.rows() == 0 || M.cols() == 0) {
      stop("Zero-dimension covariance in Hessian block!");
    }
    Eigen::LLT<Eigen::MatrixXd> lltOfA(M);
    if (lltOfA.info() != Eigen::Success) {
      ::Rf_warning("Cholesky decomposition failed at sample %d (Hessian)", i + 1);
      continue;
    }
    Eigen::MatrixXd inv_Mi = lltOfA.solve(Eigen::MatrixXd::Identity(M.rows(), M.cols()));
    H_theta_theta.noalias() += PhiiXi.transpose() * inv_Mi * PhiiXi;
  }
  return H_theta_theta;
}
