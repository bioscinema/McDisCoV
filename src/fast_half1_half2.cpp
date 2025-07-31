// fast_half1_half2.cpp

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;

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
    MatrixXd M     = as<MatrixXd>(Sigma_list[i]);

    if (M.rows() == 0 || M.cols() == 0) {
      stop("There is at least one pair of intercept and variable linearly dependent after zero separation!");
    }
    Eigen::LLT<MatrixXd> lltOfA(M);
    if (lltOfA.info() != Eigen::Success) {
      ::Rf_warning("Cholesky decomposition failed at sample %d", i + 1);
      continue;
    }
    MatrixXd inv_Mi         = lltOfA.solve(MatrixXd::Identity(M.rows(), M.cols()));
    MatrixXd tPhiiXi_invMi  = PhiiXi.transpose() * inv_Mi;
    half1.noalias()        += tPhiiXi_invMi * PhiiXi;
    half2.noalias()        += tPhiiXi_invMi * PhiiZi;
  }

  return List::create(
    Named("half1") = half1,
    Named("half2") = half2
  );
}
