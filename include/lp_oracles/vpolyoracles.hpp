// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2026 Vissarion Fisikopoulos
// Copyright (c) 2018-2026 Apostolos Chalkis
// Copyright (c) 2026      Dimitrios Pavlou

// Contributed and/or modified by Dimitrios Pavlou, as part of Google Summer of Code 2026 program

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef VPOLYORACLES_HPP
#define VPOLYORACLES_HPP

#include <vector>
#include <utility>
#include <cmath>
#include "Highs.h"
#include "lp_oracle_options.hpp"

// Builds the LP of the line p + l v against a V-polytope or zonotope with vertices in V:
//
//      V^T w + l v = p, sum(w) = 1 (V-polytope only)
//      w >= 0 (V-polytope) or -1 <= w <= 1 (zonotope), l free, objective l
//
// The matrix is stored in CSC format, so the model is passed to HiGHS in one call.
template <typename MT, typename Point>
HighsLp line_intersection_lp(MT const& V, Point const& p, Point const& v, bool zonotope)
{
    unsigned const d = v.dimension();
    unsigned const m = V.rows();
    unsigned const num_rows = zonotope ? d : d+1;

    HighsLp lp;
    lp.num_col_ = m+1;
    lp.num_row_ = num_rows;

    lp.col_cost_.assign(m+1, 0.0);
    lp.col_cost_[m] = 1.0;
    lp.col_lower_.assign(m+1, zonotope ? -1.0 : 0.0);
    lp.col_upper_.assign(m+1, 1.0);
    lp.col_lower_[m] = -kHighsInf;
    lp.col_upper_[m] = kHighsInf;

    for (unsigned i = 0; i < d; ++i) {
        lp.row_lower_.push_back((double)p[i]);
        lp.row_upper_.push_back((double)p[i]);
    }

    if (!zonotope) {
        lp.row_lower_.push_back(1.0);
        lp.row_upper_.push_back(1.0);
    }

    HighsSparseMatrix& A = lp.a_matrix_;
    A.format_ = MatrixFormat::kColwise;
    A.num_col_ = m+1;
    A.num_row_ = num_rows;
    A.start_.assign(1, 0);

    for (unsigned j = 0; j < m; ++j) {
        for (unsigned i = 0; i < d; ++i) {
            double const a = (double)V(j, i);
            if (std::abs(a) >= LPOracleZeroTol) {
                A.index_.push_back((HighsInt)i);
                A.value_.push_back(a);
            }
        }
        if (!zonotope) {
            A.index_.push_back((HighsInt)d);
            A.value_.push_back(1.0);
        }
        A.start_.push_back((HighsInt)A.index_.size());
    }

    for (unsigned i = 0; i < d; ++i) {
        double const a = (double)v[i];
        if (std::abs(a) >= LPOracleZeroTol) {
            A.index_.push_back((HighsInt)i);
            A.value_.push_back(a);
        }
    }
    A.start_.push_back((HighsInt)A.index_.size());

    return lp;
}

// Decides whether q belongs in the V-Polytope given by the vertex matrix V.
//
// The LP looks for a hyperplane seperating q from the vertices, so it
// maximizes q k - l s.t. V k <= l. Such a hyperplane exists when the optimum
// is positive, and q is inside the hull when it is not.
// @tparam MT the matrix type of V
// @tparam Point the point type
// @param V the vertex matrix, one vertex per row
// @param q the point to test
// @param opts optional callback to configure highs, see LPOracleOptions
// @return whether q belongs to V, and whether the lp was solved successfully
template <typename MT, typename Point>
LPOracleResult<bool> memLP_Vpoly(const MT& V, const Point& q,
                                 LPOracleOptions const& opts = nullptr)
{
    unsigned d = q.dimension();
    unsigned m = V.rows();

    Highs highs;
    lp_oracles_configure_highs(highs, opts);

    for (unsigned i = 0; i <= d; ++i)
        highs.addVar(-kHighsInf, kHighsInf);

    // Adds the variables for k and l. 
    std::vector<HighsInt> indices(d+1);
    std::vector<double> values(d+1);
    for (unsigned i = 0; i <= d; ++i)
        indices[i] = (HighsInt)i;

    // Forces V k <= l.
    for (unsigned i = 0; i < m; ++i) {
        for (unsigned j = 0; j < d; ++j) {
            values[j] = (double)V(i, j);
        }
        values[d] = -1.0;
        highs.addRow(-kHighsInf, 0.0, d+1, indices.data(), values.data());
    }

    // Bounds the objective so the LP stays bounded.
    for (unsigned i = 0; i < d; ++i) {
        values[i] = (double)q[i];
    }
    values[d] = -1.0;
    highs.addRow(-kHighsInf, 1.0, d+1, indices.data(), values.data());

    // Maximizes q k - l
    for (unsigned i = 0; i < d; ++i) 
        highs.changeColCost((HighsInt)i, (double)q[i]);

    highs.changeColCost((HighsInt)d, -1.0);
    highs.changeObjectiveSense(ObjSense::kMaximize);
    highs.run();
    
    if (highs.getModelStatus() != HighsModelStatus::kOptimal) {
        #ifdef VOLESTI_DEBUG
            std::cout << "Could not solve the Linear Program for membership"
                      << ", highs returned code "
                      << (int)highs.getModelStatus()
                      << std::endl;
        #endif
        return {false, false};
    }

    return {highs.getObjectiveValue() <= 0.0, true};
}

// Computes the intersection of the ray p + l v with a V-Polytope, or with
// a zonotope when zonotope is true.
//
// The point on the ray is written as a combination of the vertices, so the LP
// has the m combination weights and l as variables, with the d rows V^T w + l v = p. 
// Additionally Z-Polytope and V-Polytope conditions are applied appropriately.
// @tparam NT the number type
// @tparam MT the matrix type of V
// @tparam Point the point type
// @param V the vertex matrix, one vertex per row
// @param p the line origin
// @param v the line direction
// @param conv_comb set to the combination of weight of the intersection point
// @param maxi if true the largest lambda is computed, otherwise the smallest
// @param zonotope true when V describes a zonotope
// @param opts optional callback to configure highs, see LPOracleOptions
// @return lambda, and whether the LP was solved successfully
template <typename NT, typename MT, typename Point>
LPOracleResult<NT> intersect_line_Vpoly(MT const& V, Point const& p, Point const& v,
                                        double *conv_comb, bool maxi, bool zonotope,
                                        LPOracleOptions const& opts = nullptr)
{   
    unsigned const m = V.rows();

    Highs highs;
    lp_oracles_configure_highs(highs, opts);
    highs.passModel(line_intersection_lp(V, p, v, zonotope));
    highs.changeObjectiveSense(maxi ? ObjSense::kMaximize : ObjSense::kMinimize);
    highs.run();

    if (highs.getModelStatus() != HighsModelStatus::kOptimal) {
        #ifdef VOLESTI_DEBUG
            std::cout << "Could not solve the Linear Program for V-polytope line intersection"
                      << ", highs returned status code "
                      << (int)highs.getModelStatus()
                      << std::endl;
        #endif
        return {NT(0), false};
    }

    const auto& sol = highs.getSolution().col_value;
    for (unsigned i = 0; i < m; ++i) {
        conv_comb[i] = (std::abs(sol[i]) < LPOracleZeroTol) ? 0.0 : sol[i];
    }

    return {NT(-highs.getObjectiveValue()), true};
}

// Computes both intersections of the line p + l v with a V-polytope or zonotope.
//
// The model is built once, it is solved minimizing, then the objective sense is
// flipped and it is re-solved from the previous basis.
// @tparam NT the number type
// @tparam MT the matrix type of V
// @tparam Point the point type
// @param V the vertex matrix, one vertex per row
// @param p the line origin
// @param v the line direction
// @param zonotope true if V represents a zonotope
// @param opts optional callback to configure highs, see LPOracleOptions
// @return the pair (first, second): first from minimizing l, second from maximizing l,
// both as -l, and whether both solves succeeded.
template <typename NT, typename MT, typename Point>
LPOracleResult<std::pair<NT, NT>> intersect_double_line(MT const& V, Point const& p, Point const& v,
                                                        bool zonotope,
                                                        LPOracleOptions const& opts = nullptr)
{
    Highs highs;
    lp_oracles_configure_highs(highs, opts);
    highs.passModel(line_intersection_lp(V, p, v, zonotope));

    highs.changeObjectiveSense(ObjSense::kMinimize);
    highs.run();
    if (highs.getModelStatus() != HighsModelStatus::kOptimal)
        return {{NT(0), NT(0)}, false};
    NT const first = NT(-highs.getObjectiveValue());

    highs.changeObjectiveSense(ObjSense::kMaximize);
    highs.run();
    if (highs.getModelStatus() != HighsModelStatus::kOptimal)
        return {{NT(0), NT(0)}, false};
    NT const second = NT(-highs.getObjectiveValue());

    return {{first, second}, true};
}

// Computes both intersections of the line p + l v with a V-Polytope.
//
// The model is built once and solved twice, see intersect_double_line.
// @tparam MT the matrix type of V
// @tparam Point the point type
// @param V the vertex matrix, one vertex per row
// @param p the line origin
// @param v the line direction
// @param opts optional callback to configure highs, see LPOracleOptions
// @return the pair (lambda_min, lambda_max), and whether the LP was solved successfully
template <typename NT, typename MT, typename Point>
LPOracleResult<std::pair<NT, NT>> intersect_double_line_Vpoly(MT const& V, Point const& p, Point const& v,
                                                              LPOracleOptions const& opts = nullptr)
{
    return intersect_double_line<NT>(V, p, v, false, opts);
}
#endif