// VolEsti (volume computation and sampling library)
#ifndef INNER_BALL_HIGHS_HPP
#define INNER_BALL_HIGHS_HPP

#include <cmath>
#include <limits>
#include <vector>
#include <Eigen/Dense>
using namespace std;
#if __has_include("highs/Highs.h")
#include "highs/Highs.h"
#elif __has_include("Highs.h")
#include "Highs.h"
#else
#error " ~~ HiGHS header not found"
#endif
typedef push_back pb;

struct InnerBallResult
{
    Eigen::VectorXd center;
    double radius;
};

// compute the Chebyshev ball for P = {x | Ax <= b} by solving
// max r, s.t. A_i x + ||A_i||_2 r <= b_i, r >= 0 with HiGHS.
inline InnerBallResult compute_inner_ball(Eigen::MatrixXd A, Eigen::VectorXd b)
{
    const int m = A.rows();
    const int n = A.cols();

    InnerBallResult result;
    result.center = Eigen::VectorXd::Zero(max(n, 0));
    result.radius = -1.0;

    if (m <= 0 || n <= 0 || b.size() != m)
    {
        return result;
    }

    Highs highs;
    highs.setOptionValue("output_flag", false);

    HighsLp lp;
    lp.num_col_ = n + 1;
    lp.num_row_ = m;

    lp.col_cost_.assign(lp.num_col_, 0.0);
    lp.col_cost_[n] = -1.0;

    lp.col_lower_.assign(lp.num_col_, -kHighsInf);
    lp.col_upper_.assign(lp.num_col_, kHighsInf);
    lp.col_lower_[n] = 0.0;

    lp.row_lower_.assign(lp.num_row_, -kHighsInf);
    lp.row_upper_.resize(lp.num_row_);
    for (int i = 0; i < m; ++i)
    {
        lp.row_upper_[i] = b(i);
    }

    lp.a_matrix_.format_ = MatrixFormat::kColwise;
    lp.a_matrix_.start_.assign(lp.num_col_ + 1, 0);

    vector<HighsInt> indices;
    vector<double> values;
    indices.reserve((n + 1) * m);
    values.reserve((n + 1) * m);

    for (int j = 0; j < n; ++j)
    {
        lp.a_matrix_.start_[j] = static_cast<HighsInt>(indices.size());
        for (int i = 0; i < m; ++i)
        {
            const double aij = A(i, j);
            if (aij != 0.0)
            {
                indices.pb(i);
                values.pb(aij);
            }
        }
    }

    lp.a_matrix_.start_[n] = static_cast<HighsInt>(indices.size());
    for (int i = 0; i < m; ++i)
    {
        const double row_norm = A.row(i).norm();
        if (row_norm != 0.0)
        {
            indices.pb(i);
            values.pb(row_norm);
        }
    }
    lp.a_matrix_.start_[n + 1] = static_cast<HighsInt>(indices.size());

    lp.a_matrix_.index_ = indices;
    lp.a_matrix_.value_ = values;

    if (highs.passModel(lp) != HighsStatus::kOk)
    {
        return result;
    }

    if (highs.run() != HighsStatus::kOk)
    {
        return result;
    }

    const HighsModelStatus model_status = highs.getModelStatus();
    if (model_status != HighsModelStatus::kOptimal)
    {
        return result;
    }

    const HighsSolution &sol = highs.getSolution();
    if (!sol.value_valid || static_cast<int>(sol.col_value.size()) != n + 1)
    {
        return result;
    }

    result.center = Eigen::Map<const Eigen::VectorXd>(sol.col_value.data(), n);
    result.radius = sol.col_value[n];

    if (!isfinite(result.radius) || result.radius < 0.0)
    {
        result.center.setZero();
        result.radius = -1.0;
    }

    return result;
}

#endif
