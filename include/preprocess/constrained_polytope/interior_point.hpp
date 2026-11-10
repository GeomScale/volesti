// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2026 Vissarion Fisikopoulos
// Copyright (c) 2018-2026 Apostolos Chalkis
// Copyright (c) 2026      Dimitrios Pavlou

// Contributed and/or modified by Dimitrios Pavlou, as part of Google Summer of Code 2026 program

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef INTERIOR_POINT_HPP
#define INTERIOR_POINT_HPP

#include <cmath>
#include <vector>
#include "Highs.h"
#include "convex_bodies/constrained_polytope.hpp"

// The result of an interior point computation.
// @tparam VT the vector type.
template <typename VT>
struct InteriorPoint {
    // The interior point found.
    VT point;

    // The uniform slack of the point against every inequality.
    double slack = 0.0;

    // True if it's a valid point
    bool valid = false;
};

// Finds a point in the strict interior of the polytope by maximizing a uniform slack y.
//
// The LP solved is:
//
//      max y
//      s.t. A_eq x = b_eq
//           A_in x + ||a_i|| y <= b_in(i) for every row i
//           lb(j) + y <= x_j
//           ub(j) + y >= x_j
//           0 <= y <= 1
//
// @tparam Point the point type
// @param P the polytope
// @param tol the minimum slack for the point to count
// @return the result, meaningful only if true is returned
template <typename Point>
InteriorPoint<typename ConstrainedPolytope<Point>::VT>
find_interior_point(ConstrainedPolytope<Point> const& P,
                    double tol = 1e-9)
{
    using NT = typename ConstrainedPolytope<Point>::NT;
    using VT = typename ConstrainedPolytope<Point>::VT;
    using MT = typename ConstrainedPolytope<Point>::MT;

    MT const& A_eq = P.getEqualities();
    VT const& b_eq = P.getEqualityRHS();
    MT const& A_in = P.getInequalities();
    VT const& b_in = P.getInequalityRHS();
    VT const& lb = P.getLowerBounds();
    VT const& ub = P.getUpperBounds();
    unsigned const d = P.getDimension();
    unsigned m_eq = P.getNumEqualities();
    unsigned m_in = P.getNumInequalities();

    InteriorPoint<VT> ip;

    Highs highs;
    highs.setOptionValue("output_flag", false);
    highs.setOptionValue("presolve", "on");

    // The original variables.
    for (unsigned j = 0; j < d; ++j)
        highs.addVar(-kHighsInf, kHighsInf);

    // The slack variable
    highs.addVar(0.0, 1.0);
    HighsInt const y = (HighsInt)d;

    // Adds the equality rows (no slack needed)
    for (unsigned i = 0; i < m_eq; ++i) {
        std::vector<HighsInt> ids;
        std::vector<double> vals;

        for (typename MT::InnerIterator it(A_eq, i); it; ++it) {
            ids.push_back((HighsInt)it.col());
            vals.push_back((double)it.value());
        }

        highs.addRow((double)b_eq(i), (double)b_eq(i), (HighsInt)ids.size(), 
                     ids.data(), vals.data());
    }

    // Adds the inequality rows a_i x + || a_i || y <= b_in(i)
    for (unsigned i = 0; i < m_in; ++i) {
        std::vector<HighsInt> ids;
        std::vector<double> vals;

        for (typename MT::InnerIterator it(A_in, i); it; ++it) {
            ids.push_back((HighsInt)it.col());
            vals.push_back((double)it.value());
        }

        double norm = A_in.row(i).norm();

        if (norm > 0.0) {
            ids.push_back((HighsInt)d);
            vals.push_back(norm);
        }

        highs.addRow(-kHighsInf, (double)b_in(i), (HighsInt)ids.size(),
                     ids.data(), vals.data());

    }

    // Adds the box bounds, each with norm 1.
    for (unsigned j = 0; j < d; ++j) {
        HighsInt ids[2] = {(HighsInt)j, y};

        if (std::isfinite((double)lb(j))) {
            double vals[2] = {1.0, -1.0};
            highs.addRow((double)lb(j), kHighsInf, 2, ids, vals);
        }

        if (std::isfinite((double)ub(j))) {
            double vals[2] = {1.0, 1.0};
            highs.addRow(-kHighsInf, (double)ub(j), 2, ids, vals);
        }
    }

    // Maximizes the slack
    highs.changeColCost((HighsInt)d, -1.0);
    highs.run();

    if (highs.getModelStatus() != HighsModelStatus::kOptimal)
        return ip;

    ip.slack = -highs.getObjectiveValue();

    if (ip.slack < tol)
        return ip;

    auto const& sol = highs.getSolution().col_value;
    ip.point.resize(d);

    for (unsigned j = 0; j < d; ++j)
        ip.point(j) = (NT)sol[j];

    ip.valid = true;
    return ip;
}
#endif