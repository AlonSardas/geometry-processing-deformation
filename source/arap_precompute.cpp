#include "arap_precompute.h"
#include "min_quad_with_fixed_helpers.h"
#include <igl/arap_linear_block.h>
#include <igl/cotmatrix.h>
#include <igl/min_quad_with_fixed.h>

typedef Eigen::Triplet<double> T;

void compute_K_matrix(int n, const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                      Eigen::SparseMatrix<double> &K);

void push_half_vertex(const Eigen::MatrixXd &V, int i, int j, int other,
                      double c_ij, std::vector<T> &triplets);

void arap_precompute(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                     const Eigen::VectorXi &b,
                     igl::min_quad_with_fixed_data<double> &data,
                     Eigen::SparseMatrix<double> &K) {
  Eigen::SparseMatrix<double> L;
  igl::cotmatrix(V, F, L);

  precompute_min_quad_with_fixed(L, b, data);
  compute_K_matrix(data.n, V, F, K);
}

void compute_K_matrix(int n, const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                      Eigen::SparseMatrix<double> &K) {
  // Create K
  K.resize(n, n * 3);

  Eigen::MatrixXd C;
  igl::cotmatrix_entries(V, F, C);

  std::vector<T> triplets;
  triplets.reserve(F.rows() * 3 * 3 * 6);
  for (int f = 0; f < F.rows(); ++f) {
    push_half_vertex(V, F(f, 0), F(f, 1), F(f, 2), C(f, 2), triplets);
    push_half_vertex(V, F(f, 1), F(f, 2), F(f, 0), C(f, 0), triplets);
    push_half_vertex(V, F(f, 2), F(f, 0), F(f, 1), C(f, 1), triplets);
  }

  assert(triplets.size() == F.rows() * 3 * 3 * 6);
  K.setFromTriplets(triplets.begin(), triplets.end());
}

void push_half_vertex(const Eigen::MatrixXd &V, int i, int j, int other,
                      double c_ij, std::vector<T> &triplets) {
  Eigen::RowVector3d e_ij = c_ij * (V.row(i) - V.row(j)) / 6.0;

  int ks[] = {i, j, other};
  for (int k_idx = 0; k_idx < 3; ++k_idx) {
    triplets.emplace_back(i, 3 * ks[k_idx] + 0, e_ij(0));
    triplets.emplace_back(i, 3 * ks[k_idx] + 1, e_ij(1));
    triplets.emplace_back(i, 3 * ks[k_idx] + 2, e_ij(2));
    triplets.emplace_back(j, 3 * ks[k_idx] + 0, -e_ij(0));
    triplets.emplace_back(j, 3 * ks[k_idx] + 1, -e_ij(1));
    triplets.emplace_back(j, 3 * ks[k_idx] + 2, -e_ij(2));
  }
}
