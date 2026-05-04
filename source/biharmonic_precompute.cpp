#include "biharmonic_precompute.h"
#include "min_quad_with_fixed_helpers.h"
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/min_quad_with_fixed.h>

void biharmonic_precompute(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                           const Eigen::VectorXi &b,
                           igl::min_quad_with_fixed_data<double> &data) {
  Eigen::SparseMatrix<double> M, L, Q;
  igl::cotmatrix(V, F, L);
  igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);

  // Assert M is diagonal
  for (int k = 0; k < M.outerSize(); k++)
    for (Eigen::SparseMatrix<double>::InnerIterator it(M, k); it; ++it)
      assert(it.row() == it.col() && "M is not diagonal");

  // Invert diagonal
  Eigen::SparseMatrix<double> M_inv(M.rows(), M.cols());
  for (int i = 0; i < M.rows(); i++) {
    M_inv.coeffRef(i, i) = 1.0 / M.coeff(i, i);
  }

  Q = L.transpose() * M_inv * L;

  precompute_min_quad_with_fixed(Q, b, data);
}
