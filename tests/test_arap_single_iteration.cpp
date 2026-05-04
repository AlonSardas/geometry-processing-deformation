#include "arap_precompute.h"
#include "arap_single_iteration.h"
#include <Eigen/Dense>
#include <gtest/gtest.h>
#include <igl/min_quad_with_fixed.h>

TEST(TestArapStep, TestNoDeformation) {
  // Grid mesh: 3x3 vertices, 8 triangles
  Eigen::MatrixXd V(9, 3);
  V << 0, 0, 0, //
      1, 0, 0,  //
      2, 0, 0,  //
      0, 1, 0,  //
      1, 1, 0,  //
      2, 1, 0,  //
      0, 2, 0,  //
      1, 2, 0,  //
      2, 2, 0;

  Eigen::MatrixXi F(8, 3);
  F << 0, 1, 4, 0, 4, 3, 1, 2, 5, 1, 5, 4, 3, 4, 7, 3, 7, 6, 4, 5, 8, 4, 8, 7;

  Eigen::VectorXi b(4);
  b << 0, 2, 5, 7;
  igl::min_quad_with_fixed_data<double> data;
  Eigen::SparseMatrix<double> K;
  arap_precompute(V, F, b, data, K);

  Eigen::MatrixXd bc(b.size(), 3);
  bc << V.row(b(0)), V.row(b(1)), V.row(b(2)), V.row(b(3));

  Eigen::MatrixXd U = V;
  arap_single_iteration(data, K, bc, U);

  //   EXPECT_TRUE(U.isApprox(V, 1e-6));
  for (int i = 0; i < U.rows(); i++)
    for (int j = 0; j < U.cols(); j++)
      EXPECT_NEAR(U(i, j), V(i, j), 1e-6) << "at (" << i << "," << j << ")";
}

TEST(TestArapStep, TestRigidRotation) {
  Eigen::MatrixXd V(9, 3);
  V << 0, 0, 0, 1, 0, 0, 2, 0, 0, 0, 1, 0, 1, 1, 0, 2, 1, 0, 0, 2, 0, 1, 2, 0,
      2, 2, 0;
  Eigen::MatrixXi F(8, 3);
  F << 0, 1, 4, 0, 4, 3, 1, 2, 5, 1, 5, 4, 3, 4, 7, 3, 7, 6, 4, 5, 8, 4, 8, 7;

  Eigen::VectorXi b(3);
  b << 0, 2, 5;
  igl::min_quad_with_fixed_data<double> data;
  Eigen::SparseMatrix<double> K;
  arap_precompute(V, F, b, data, K);

  // 90 degree rotation around Z axis
  double angle = M_PI / 2.0;
  Eigen::Matrix3d R;
  R << cos(angle), -sin(angle), 0, sin(angle), cos(angle), 0, 0, 0, 1;

  // Only boundary conditions are rotated
  Eigen::MatrixXd bc(3, 3);
  bc << (R * V.row(0).transpose()).transpose(),
      (R * V.row(2).transpose()).transpose(),
      (R * V.row(5).transpose()).transpose();

  // Initialize U from original V
  Eigen::MatrixXd U = V;
  for (int i = 0; i < 1000; ++i) {
    arap_single_iteration(data, K, bc, U);
  }

  Eigen::MatrixXd V_rot = (R * V.transpose()).transpose();
  for (int i = 0; i < U.rows(); i++)
    for (int j = 0; j < U.cols(); j++)
      EXPECT_NEAR(U(i, j), V_rot(i, j), 1e-6) << "at (" << i << "," << j << ")";
}