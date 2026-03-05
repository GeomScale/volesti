#pragma once

#include <Eigen/Dense>
#include "Highs.h"

namespace volesti {

template <typename NT>
bool compute_chebyshev_highs(
    const Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>& A,
    const Eigen::Matrix<NT, Eigen::Dynamic, 1>& b,
    Eigen::Matrix<NT, Eigen::Dynamic, 1>& center,
    NT& radius
) {
    int m = A.rows();
    int d = A.cols();

    Highs highs;
    highs.silent();

    // Variables: x1..xd, r
    int nvar = d + 1;

    HighsLp lp;
    lp.num_col_ = nvar;
    lp.num_row_ = m;

    lp.col_cost_.assign(nvar, 0.0);
    lp.col_cost_[d] = -1.0;  // maximize r

    lp.col_lower_.assign(nvar, -kHighsInf);
    lp.col_upper_.assign(nvar, kHighsInf);

    lp.row_lower_.assign(m, -kHighsInf);
    lp.row_upper_.assign(b.data(), b.data() + m);

    std::vector<int> astart(nvar + 1, 0);
    std::vector<int> aindex;
    std::vector<double> avalue;

    for (int j = 0; j < d; ++j) {
        astart[j] = aindex.size();
        for (int i = 0; i < m; ++i) {
            if (A(i, j) != 0) {
                aindex.push_back(i);
                avalue.push_back(A(i, j));
            }
        }
    }

    // r column
    astart[d] = aindex.size();
    for (int i = 0; i < m; ++i) {
        NT normAi = A.row(i).norm();
        aindex.push_back(i);
        avalue.push_back(normAi);
    }

    astart[nvar] = aindex.size();

    lp.a_matrix_.start_ = astart;
    lp.a_matrix_.index_ = aindex;
    lp.a_matrix_.value_ = avalue;

    highs.passModel(lp);

    if (highs.run() != HighsStatus::kOk) return false;

    auto sol = highs.getSolution();
    if (!sol.value_valid) return false;

    center = Eigen::Map<const Eigen::VectorXd>(sol.col_value.data(), d);
    radius = sol.col_value[d];

    return true;
}

}