// VolEsti (volume computation and sampling library)
// Converts .mm sparse format to volesti HPolytope
#ifndef MATRIX_MARKET_READER_H
#define MATRIX_MARKET_READER_H

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include "convex_bodies/hpolytope.h"

template <typename NT>
Eigen::SparseMatrix<NT>
read_matrix_market_sparse(const std::string& filename)
{
    std::ifstream file(filename);
    if (!file.is_open())
        throw std::runtime_error("Cannot open file: " + filename);

    std::string line;

    while (std::getline(file, line))
        if (!line.empty() && line[0] != '%') break;

    int rows, cols, nnz;
    std::istringstream header(line);
    if (!(header >> rows >> cols >> nnz))
        throw std::runtime_error("Bad header in: " + filename);

    std::vector<Eigen::Triplet<NT>> triplets;
    triplets.reserve(nnz);

    for (int i = 0; i < nnz; i++) {
        int r, c; NT val;
        if (!(file >> r >> c >> val))
            throw std::runtime_error("Unexpected end of file: " + filename);
        triplets.emplace_back(r-1, c-1, val); 
    }

    Eigen::SparseMatrix<NT> M(rows, cols);
    M.setFromTriplets(triplets.begin(), triplets.end());
    return M;
}

template <typename NT>
Eigen::Matrix<NT, Eigen::Dynamic, 1>
read_matrix_market_vector(const std::string& filename)
{
    std::ifstream file(filename);
    if (!file.is_open())
        throw std::runtime_error("Cannot open file: " + filename);

    std::string line;
    while (std::getline(file, line))
        if (!line.empty() && line[0] != '%') break;

    int rows, cols, nnz;
    std::istringstream header(line);
    if (!(header >> rows >> cols >> nnz))
        throw std::runtime_error("Bad header in: " + filename);

    Eigen::Matrix<NT, Eigen::Dynamic, 1> v =
        Eigen::Matrix<NT, Eigen::Dynamic, 1>::Zero(rows);

    for (int i = 0; i < nnz; i++) {
        int r, c; NT val;
        if (!(file >> r >> c >> val))
            throw std::runtime_error("Unexpected end of file: " + filename);
        v(r-1) = val; 
    }
    return v;
}

// A_file      : constraint matrix (.mm)  — sparse
// bounds_file : variable upper bounds (.mm) — one per variable
template <typename Point>
HPolytope<Point> matrix_market_to_hpolytope(
    const std::string& A_file,
    const std::string& bounds_file)
{
    typedef typename Point::FT NT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;

    // Read sparse constraint matrix A
    Eigen::SparseMatrix<NT> A_sparse =
        read_matrix_market_sparse<NT>(A_file);

    int num_constraints = A_sparse.rows();
    int num_vars        = A_sparse.cols();

    // b = zero vector for standard LP feasibility (Ax <= 0)
    VT b = VT::Zero(num_constraints);

    // Read bounds as sparse matrix (num_vars × 2)
    // col 1 = lower bounds, col 2 = upper bounds
    Eigen::SparseMatrix<NT> bounds_sparse =
        read_matrix_market_sparse<NT>(bounds_file);

    // Convert to dense for easy access
    MT bounds_dense = MT(bounds_sparse);

    // bounds_dense has shape (num_bounds_rows × 2)
    // Some variables may have no entry = unbounded
    // We only add rows for finite upper bounds (< 1e8)
    NT INF_THRESHOLD = NT(1e8);

    // Count finite upper bounds
    std::vector<std::pair<int,NT>> finite_bounds;
    for (int i = 0; i < bounds_dense.rows(); i++) {
        NT ub = (bounds_dense.cols() >= 2) ? bounds_dense(i, 1) : bounds_dense(i, 0);
        if (ub < INF_THRESHOLD) {
            finite_bounds.push_back({i, ub});
        }
    }

    // Build extended system
    int total_rows = num_constraints + (int)finite_bounds.size();
    MT A_extended  = MT::Zero(total_rows, num_vars);
    VT b_extended  = VT::Zero(total_rows);

    // Original constraints
    A_extended.topRows(num_constraints) = MT(A_sparse);
    b_extended.head(num_constraints)    = b;

    // Add finite upper bound rows: x_i <= ub
    for (int k = 0; k < (int)finite_bounds.size(); k++) {
        int var_idx = finite_bounds[k].first;
        NT  ub      = finite_bounds[k].second;
        if (var_idx < num_vars) {
            A_extended(num_constraints + k, var_idx) = NT(1);
            b_extended(num_constraints + k)          = ub;
        }
    }

    unsigned int dim = (unsigned int) num_vars;
    return HPolytope<Point>(dim, A_extended, b_extended);
}
#endif
