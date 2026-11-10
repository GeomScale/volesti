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
#include <memory>
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
            A.index_.push_back((HighsInt)i);
            A.value_.push_back((double)V(j,i));
        }
        if (!zonotope) {
            A.index_.push_back((HighsInt)d);
            A.value_.push_back(1.0);
        }
        A.start_.push_back((HighsInt)A.index_.size());
    }

    for (unsigned i = 0; i < d; ++i) {
        A.index_.push_back((HighsInt)i);
        A.value_.push_back((double)v[i]);
    }
    A.start_.push_back((HighsInt)A.index_.size());

    return lp;
}

// Keeps one HiGHS model of the line intersection LP across calls.
//
// Between calls only the line changes, so the first load builds the model and
// later loads only update p (the row bounds) and v (the column of l). HiGHS then
// re-solves from the previous basis, which is close to optimal for nearby points.
// A copy starts empty and builds its own model, reset must be called when V changes.
class LineIntersectionSolver {
    private:
        LPOracleOptions _opts;
        std::unique_ptr<Highs> _highs;

    public:
        // Creates an empty solver, the model is built on the first load.
        // @param opts optional callback to configure highs, see LPOracleOptions
        explicit LineIntersectionSolver(LPOracleOptions opts = nullptr) : _opts(std::move(opts)) {}

        // Copies only the options, the copy builds its own model on the first load.
        // @param solver the solver to copy the options from
        LineIntersectionSolver(LineIntersectionSolver const& solver) : _opts(solver._opts) {}

        // Copies only the options and drops the current model.
        // @param solver the solver to copy options from
        // @return this solver
        LineIntersectionSolver& operator=(LineIntersectionSolver const& solver)
        {
            _opts = solver._opts;
            reset();
            return *this;
        }

        // Drops the model, so the next load builds it again from V.
        void reset() { _highs.reset(); }

        // Loads the line p + l v, building the model on the first call.
        //
        // V and zonotope are only read when the model is built, later calls
        // change the d row bounds to p and the column of l to v.
        // @tparam MT the matrix type of V
        // @tparam Point the point type
        // @param V the vertex matrix, one vertex (or generator) per row
        // @param p the line origin
        // @param v the line direction
        // @param zonotope true if V holds the generators of a zonotope
        template <typename MT, typename Point>
        void load(MT const& V, Point const& p, Point const& v, bool zonotope)
        {
            if (!_highs) {
                _highs.reset(new Highs());
                lp_oracles_configure_highs(*_highs, _opts);
                _highs->passModel(line_intersection_lp(V, p, v, zonotope));
                return;
            }
            HighsInt const l = (HighsInt)V.rows();
            for (unsigned i = 0; i < v.dimension(); ++i) {
                _highs->changeRowBounds((HighsInt)i, (double)p[i], (double)p[i]);
                _highs->changeCoeff((HighsInt)i, l, (double)v[i]);
            }
        }

        // Solves the loaded model from the previous basis.
        //
        // A failed solve drops the model, so a bad basis is never reused.
        // @param maxi if true l is maximized, otherwise minimized
        // @return the HiGHS model status, kOptimal when the solve succeeded
        HighsModelStatus run(bool maxi)
        {
            _highs->changeObjectiveSense(maxi ? ObjSense::kMaximize : ObjSense::kMinimize);
            _highs->run();

            HighsModelStatus const status = _highs->getModelStatus();
            if (status != HighsModelStatus::kOptimal)
                reset();
            return status;
        }

        // @return the model, for reading the solution after a successful run
        Highs const& highs() const { return *_highs; }
};

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
// @param solver the model of the polytope, reused across calls, see LineIntersectionSolver
// @param V the vertex matrix, one vertex per row
// @param p the line origin
// @param v the line direction
// @param conv_comb set to the combination of weight of the intersection point
// @param maxi if true the largest lambda is computed, otherwise the smallest
// @param zonotope true when V describes a zonotope
// @return lambda, and whether the LP was solved successfully
template <typename NT, typename MT, typename Point>
LPOracleResult<NT> intersect_line_Vpoly(LineIntersectionSolver& solver,
                                        MT const& V, Point const& p, Point const& v,
                                        double *conv_comb, bool maxi, bool zonotope)
{   
    unsigned const m = V.rows();

    solver.load(V, p, v, zonotope);
    HighsModelStatus const status = solver.run(maxi);

    if (status != HighsModelStatus::kOptimal) {
        #ifdef VOLESTI_DEBUG
            std::cout << "Could not solve the Linear Program for V-polytope line intersection"
                      << ", highs returned status code "
                      << (int)status
                      << std::endl;
        #endif
        return {NT(0), false};
    }

    const auto& sol = solver.highs().getSolution().col_value;
    for (unsigned i = 0; i < m; ++i) {
        conv_comb[i] = (std::abs(sol[i]) < LPOracleZeroTol) ? 0.0 : sol[i];
    }

    return {NT(-solver.highs().getObjectiveValue()), true};
}

// Computes the intersection of the ray p + l v with a V-polytope or zonotope,
// with a new model for this call only, see the overload above.
// @param opts optional callback to configure highs, see LPOracleOptions
// @return lambda, and whether the LP was solved successfully
template <typename NT, typename MT, typename Point>
LPOracleResult<NT> intersect_line_Vpoly(MT const& V, Point const& p, Point const& v,
                                        double *conv_comb, bool maxi, bool zonotope,
                                        LPOracleOptions const& opts = nullptr)
{
    LineIntersectionSolver solver(opts);
    return intersect_line_Vpoly<NT>(solver, V, p, v, conv_comb, maxi, zonotope);
}

// Computes both intersections of the line p + l v with a V-polytope or zonotope.
//
// The model is built once, it is solved minimizing, then the objective sense is
// flipped and it is re-solved from the previous basis.
// @tparam NT the number type
// @tparam MT the matrix type of V
// @tparam Point the point type
// @param solver the model of the polytope, reused across calls, see LineIntersectionSolver
// @param V the vertex matrix, one vertex per row
// @param p the line origin
// @param v the line direction
// @param zonotope true if V represents a zonotope
// @return the pair (first, second): first from minimizing l, second from maximizing l,
// both as -l, and whether both solves succeeded.
template <typename NT, typename MT, typename Point>
LPOracleResult<std::pair<NT, NT>> intersect_double_line(LineIntersectionSolver& solver,
                                                        MT const& V, Point const& p, Point const& v,
                                                        bool zonotope)
{
    solver.load(V, p, v, zonotope);

    if (solver.run(false) != HighsModelStatus::kOptimal)
        return {{NT(0), NT(0)}, false};
    NT const first = NT(-solver.highs().getObjectiveValue());

    if (solver.run(true) != HighsModelStatus::kOptimal)
        return {{NT(0), NT(0)}, false};
    NT const second = NT(-solver.highs().getObjectiveValue());

    return {{first, second}, true};
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
    LineIntersectionSolver solver(opts);
    return intersect_double_line<NT>(solver, V, p, v, zonotope);
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

// Computes both intersections of the line p + l v with a V-Polytope.
//
// The model is built once and solved twice, see intersect_double_line.
// @tparam MT the matrix type of V
// @tparam Point the point type
// @param solver the model of the polytope, reused across calls, see LineIntersectionSolver
// @param V the vertex matrix, one vertex per row
// @param p the line origin
// @param v the line direction
// @return the pair (lambda_min, lambda_max), and whether the LP was solved successfully
template <typename NT, typename MT, typename Point>
LPOracleResult<std::pair<NT, NT>> intersect_double_line_Vpoly(LineIntersectionSolver& solver,
                                                              MT const& V, Point const& p, Point const& v)
{
    return intersect_double_line<NT>(solver, V, p, v, false);
}
#endif