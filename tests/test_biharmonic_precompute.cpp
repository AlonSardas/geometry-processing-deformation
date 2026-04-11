#include "biharmonic_precompute.h"
#include <gtest/gtest.h>
#include <igl/min_quad_with_fixed.h>

TEST(TestBiharmonicPrecompute, TestSanity) {
  // Grid mesh: 3x3 vertices, 8 triangles
  Eigen::MatrixXd V(9, 3);
  V << 0, 0, 0, 1, 0, 0, 2, 0, 0, 0, 1, 0, 1, 1, 0, 2, 1, 0, 0, 2, 0, 1, 2, 0,
      2, 2, 0;

  Eigen::MatrixXi F(8, 3);
  F << 0, 1, 4, 0, 4, 3, 1, 2, 5, 1, 5, 4, 3, 4, 7, 3, 7, 6, 4, 5, 8, 4, 8, 7;

  // 4 handle vertices at corners
  Eigen::VectorXi b(2);
  b << 0, 2;
  igl::min_quad_with_fixed_data<double> data;
  biharmonic_precompute(V, F, b, data);

  ASSERT_EQ(data.Auu.rows(), 7);
  ASSERT_EQ(data.Auu.cols(), 7);
  ASSERT_EQ(data.preY.rows(), 7);
  ASSERT_EQ(data.preY.cols(), 2);
}
