#include "min_quad_with_fixed_helpers.h"
#include <igl/min_quad_with_fixed.h>

void precompute_min_quad_with_fixed(
    const Eigen::SparseMatrix<double> &M, const Eigen::VectorXi &b,
    igl::min_quad_with_fixed_data<double> &data) {
  // See min_quad_with_fixed_data struct structure at
  // https://github.com/libigl/libigl/blob/main/include/igl/min_quad_with_fixed.h#L165

  data.n = M.rows();
  data.known = b;

  std::vector<bool> is_known(M.rows(), false);
  for (int i = 0; i < b.size(); i++) {
    is_known[b(i)] = true;
  }

  data.unknown.resize(M.rows() - b.size());
  int idx = 0;
  for (int i = 0; i < M.rows(); i++) {
    if (!is_known[i]) {
      data.unknown(idx++) = i;
    }
  }

  igl::slice(M, data.unknown, data.unknown, data.Auu);
  igl::slice(M, data.unknown, data.known, data.preY);
  data.preY = -data.preY;

  data.ldlt.compute(data.Auu);
  if (data.ldlt.info() != Eigen::Success) {
    throw std::runtime_error("LDLT factorization failed");
  }
}