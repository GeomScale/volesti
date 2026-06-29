#pragma once

#include <Eigen/Dense>
#include <fstream>
#include <vector>
#include <stdexcept>

namespace volesti {
namespace io {

struct HPolytopeData {
    Eigen::MatrixXd A;
    Eigen::VectorXd b;
};

HPolytopeData import_lp(const std::string& filename)
{
    std::ifstream in(filename);
    if (!in)
        throw std::runtime_error("Cannot open file");

    int n, m;
    in >> n >> m;   // n = variables, m = equality constraints

    // Read Aeq matrix (m x n)
    Eigen::MatrixXd Aeq(m, n);
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < n; ++j)
            in >> Aeq(i, j);

    // Read beq vector
    Eigen::VectorXd beq(m);
    for (int i = 0; i < m; ++i)
        in >> beq(i);

    // Read lower bounds
    Eigen::VectorXd lower(n);
    for (int i = 0; i < n; ++i)
        in >> lower(i);

    // Read upper bounds
    Eigen::VectorXd upper(n);
    for (int i = 0; i < n; ++i)
        in >> upper(i);

    // Build inequality system Ax ≤ b
    std::vector<Eigen::VectorXd> rows;
    std::vector<double> bounds;

    // --- Convert equalities Ax = b ---
    for (int i = 0; i < m; ++i) {
        rows.push_back(Aeq.row(i));
        bounds.push_back(beq(i));

        rows.push_back(-Aeq.row(i));
        bounds.push_back(-beq(i));
    }

    // --- Convert bounds l ≤ x ≤ u ---
    for (int i = 0; i < n; ++i) {

        // x_i ≤ upper
        Eigen::VectorXd r = Eigen::VectorXd::Zero(n);
        r(i) = 1.0;
        rows.push_back(r);
        bounds.push_back(upper(i));

        // -x_i ≤ -lower
        r(i) = -1.0;
        rows.push_back(r);
        bounds.push_back(-lower(i));
    }

    // Convert vectors → matrices
    int k = rows.size();
    Eigen::MatrixXd A(k, n);
    Eigen::VectorXd b(k);

    for (int i = 0; i < k; ++i) {
        A.row(i) = rows[i];
        b(i) = bounds[i];
    }

    return {A, b};
}

} // namespace io
} // namespace volesti