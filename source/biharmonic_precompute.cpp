#include "biharmonic_precompute.h"
#include <igl/cotmatrix.h>
#include <igl/invert_diag.h>
#include <igl/massmatrix.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/slice.h>
#include <stdexcept>

void biharmonic_precompute(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                           const Eigen::VectorXi &b,
                           igl::min_quad_with_fixed_data<double> &data) {
  // See min_quad_with_fixed_data struct structure at
  // https://github.com/libigl/libigl/blob/main/include/igl/min_quad_with_fixed.h#L165
  data.n = V.rows();
  data.Auu_pd = true;
  data.Auu_sym = true;
  data.known = b;

  data.unknown.resize(V.rows() - b.size());
  std::vector<bool> is_known(V.rows(), false);
  for (int i = 0; i < b.size(); i++) {
    is_known[b(i)] = true;
  }

  int idx = 0;
  for (int i = 0; i < V.rows(); i++) {
    if (!is_known[i]) {
      data.unknown(idx) = i;
      idx++;
    }
  }

  Eigen::SparseMatrix<double> M, L, Q;
  igl::cotmatrix(V, F, L);
  igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);

  // Assert M is diagonal
  for (int k = 0; k < M.outerSize(); k++)
    for (Eigen::SparseMatrix<double>::InnerIterator it(M, k); it; ++it)
      assert(it.row() == it.col() && "M is not diagonal");

  // Invert diagonal
  Eigen::SparseMatrix<double> M_inv(M.rows(), M.cols());
  for (int i = 0; i < M.rows(); i++)
    M_inv.coeffRef(i, i) = 1.0 / M.coeff(i, i);

  Q = L.transpose() * M_inv * L;

  igl::slice(Q, data.unknown, data.unknown, data.Auu);
  igl::slice(Q, data.unknown, data.known, data.preY);
  data.preY = -data.preY;
  data.llt.compute(data.Auu);
  if (data.llt.info() != Eigen::Success) {
    throw std::runtime_error(
        "LLT factorization failed - matrix may not be positive definite");
  }
}
