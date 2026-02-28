#pragma once

#include <Eigen/Dense>
#include <stdexcept>

namespace volesti {
namespace io {

inline void validate(
    const Eigen::MatrixXd& A,
    const Eigen::VectorXd& b)
{
    if (A.rows() != b.size())
        throw std::runtime_error("Dimension mismatch");

    if (A.cols() == 0)
        throw std::runtime_error("Zero-dimensional polytope");

    if (!A.allFinite() || !b.allFinite())
        throw std::runtime_error("Non-finite values detected");
}

} // namespace io
} // namespace volesti