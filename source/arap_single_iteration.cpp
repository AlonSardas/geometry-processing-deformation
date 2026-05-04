#include "../include/arap_single_iteration.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/polar_svd3x3.h>
#include <igl/slice_into.h>

void arap_single_iteration(const igl::min_quad_with_fixed_data<double> &data,
                           const Eigen::SparseMatrix<double> &K,
                           const Eigen::MatrixXd &bc, Eigen::MatrixXd &U) {
  // Local step, finding best R
  Eigen::MatrixXd C = U.transpose() * K;
  Eigen::MatrixXd R(3 * data.n, 3);
  for (int i = 0; i < data.n; ++i) {
    // Not sure why the extra transpose here...
    Eigen::Matrix3d Ci = C.block<3, 3>(0, i * 3).transpose();
    Eigen::Matrix3d Ri;
    igl::polar_svd3x3(Ci, Ri);
    R.block<3, 3>(i * 3, 0) = Ri;
  }

  // Global step, vertices
  Eigen::MatrixXd B = K * R;

  Eigen::MatrixXd Y = data.preY * bc - 2.0 * B(data.unknown, Eigen::all);

  Eigen::MatrixXd U_unknown = U(data.unknown, Eigen::all);
  Eigen::MatrixXd LHS = data.Auu * U_unknown;
  Eigen::MatrixXd solved = data.ldlt.solve(Y);
  U(data.unknown, Eigen::all) = solved;
  igl::slice_into(bc, data.known, 1, U);
}
