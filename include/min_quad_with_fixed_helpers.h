#ifndef MIN_QUAD_WITH_FIXED_HELPERS_H
#define MIN_QUAD_WITH_FIXED_HELPERS_H

#include <Eigen/Sparse>

namespace igl {
template <typename T> struct min_quad_with_fixed_data;
}

void precompute_min_quad_with_fixed(
    const Eigen::SparseMatrix<double> &M, const Eigen::VectorXi &b,
    igl::min_quad_with_fixed_data<double> &data);

#endif