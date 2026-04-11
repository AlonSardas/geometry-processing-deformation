#include "biharmonic_solve.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/slice_into.h>
#include <iostream>

void biharmonic_solve(const igl::min_quad_with_fixed_data<double> &data,
                      const Eigen::MatrixXd &bc, Eigen::MatrixXd &D) {
  Eigen::MatrixXd Y = data.preY * bc;
  Eigen::MatrixXd U = data.llt.solve(Y);

  D.resize(data.n, 3);
  igl::slice_into(bc, data.known, 1, D);
  igl::slice_into(U, data.unknown, 1, D);
}
