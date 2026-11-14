// VolEsti (volume computation and sampling library)
// Matrix Market format reader for Netlib LP benchmarks
// Converts .mm sparse format to volesti HPolytope
#ifndef MATRIX_MARKET_READER_H
#define MATRIX_MARKET_READER_H

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>
#include <Eigen/Eigen>
#include "convex_bodies/hpolytope.h"
#include "cartesian_geom/cartesian_kernel.h"

// Step 1: Read a Matrix Market .mm file into a dense Eigen matrix
template <typename NT>
Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>
read_matrix_market(const std::string& filename) {

    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: cannot open file " << filename << std::endl;
        return Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>();
    }

    std::string line;

    // Skip comment lines starting with %
    while (std::getline(file, line)) {
        if (line[0] != '%') break;
    }

    // Read dimensions: rows cols non_zeros
    int rows, cols, nnz;
    std::istringstream header(line);
    header >> rows >> cols >> nnz;

    // Initialize matrix with zeros
    Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> M =
        Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>::Zero(rows, cols);

    // Read each non-zero entry: row col value (1-indexed)
    for (int i = 0; i < nnz; i++) {
        int r, c;
        NT val;
        file >> r >> c >> val;
        M(r-1, c-1) = val;  // convert to 0-indexed
    }

    file.close();
    return M;
}

// Step 2: Convert Matrix Market files to HPolytope
// A_file: constraint matrix (.mm)
// b_file: bounds file (.mm)
template <typename Point>
HPolytope<Point> matrix_market_to_hpolytope(
    const std::string& A_file,
    const std::string& b_file)
{
    typedef typename Point::FT NT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;

    // Read constraint matrix A
    MT A = read_matrix_market<NT>(A_file);
    int num_constraints = A.rows();
    int num_vars = A.cols();

    // Read bounds matrix (rows x 2: col1=lower, col2=upper)
    MT bounds = read_matrix_market<NT>(b_file);

    // Build b vector from upper bounds (col index 1 = second column)
    VT b(num_constraints);
    for (int i = 0; i < num_constraints; i++) {
        // Use 0 as default RHS (standard LP feasibility region)
        b(i) = NT(0);
    }

    // Add variable bound constraints: x_i <= ub_i
    // bounds matrix: row = variable index, col 1 = upper bound
    int num_bound_rows = bounds.rows();
    MT A_extended(num_constraints + num_bound_rows, num_vars);
    VT b_extended(num_constraints + num_bound_rows);

    A_extended.topRows(num_constraints) = A;
    b_extended.head(num_constraints) = b;

    // Each variable upper bound becomes: x_i <= ub
    for (int i = 0; i < num_bound_rows; i++) {
        A_extended.row(num_constraints + i).setZero();
        A_extended(num_constraints + i, i) = NT(1);
        b_extended(num_constraints + i) = bounds(i, 1); // upper bound
    }

    unsigned int dim = num_vars;
    return HPolytope<Point>(dim, A_extended, b_extended);
}

#endif // MATRIX_MARKET_READER_H

