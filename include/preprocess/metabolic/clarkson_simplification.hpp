// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2026 Vissarion Fisikopoulos
// Copyright (c) 2018-2026 Apostolos Chalkis
// Copyright (c) 2026      Dimitrios Pavlou

// Contributed and/or modified by Dimitrios Pavlou, as part of Google Summer of Code 2026 program

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef METABOLIC_CLARKSON_SIMPLIFICATION_HPP
#define METABOLIC_CLARKSON_SIMPLIFICATION_HPP

#include <iostream>
#include <array>
#include <vector>
#include <set>
#include <queue>
#include <utility>
#include <random>
#include <iterator>
#include <cmath>
#include <limits>
#include <Eigen/Eigen>
#include "convex_bodies/metabolic_polytope.hpp"
#include "preprocess/metabolic/exhaustive_simplification.hpp"
#include "Highs.h"

namespace clarkson_simplification {
    // Reuses the verbosity level and logging function from exhaustive simplification.
    using exhaustive_simplification::VerbosityLevel;
    using exhaustive_simplification::log_diagnostics;

    // Configuration parameters controlling the simplification process.
    struct Config {
        // Tolerance for marking a bound as redundant, a bound is relaxed if
        // moving it does not change the optimum by more than this quantity.
        double facet_tolerance = 1e-7;

        // Tolerance for marking a dimension as degenerate. A dimension is
        // fixed when its max and min differ less than this quantity.
        double dim_tolerance = 1e-7;

        // The error tolerance for the interior point.
        double interior_tolerance = 1e-9;
        
        // The error tolerance for the ray shooting stage of clarkson.
        double ray_tolerance = 1e-9;

        // The gap by which the bound is relaxed in clarkson's lp test.
        double relaxation_gap = 1.0;

        // The bound on the number of failed iteration's in clarkson.
        unsigned failed_iter_count = 50;
        
        // The seed used by clarkson to select inequalities.
        unsigned clarkson_seed = 0;

        // If true, degenerate dimensions are fixed before the redundancy removal.
        // Clarkson needs an interior point to work, which is difficult to compute
        // for a polytope with degenerate dimensions, so having this turned off
        // will usually send the run down the exhaustive fallback.
        bool fix_dimensions = false;
        
        // The primal feasibility tolerance of the LP solver.
        double primal_feasibility_tol = 1e-7;

        // The dual feasibility tolerance of the LP solver.
        double dual_feasibility_tol = 1e-7;

        // Iteration limit for the simplex solves.
        double simplex_iter_limit = 1000;

        // Solver time limit.
        double time_limit = 200;

        // If true, it prints diagnostic information about the simplification
        // process.
        VerbosityLevel verbosity = VerbosityLevel::Silent;

        // Where the diagnostic information is printed. Defaults to std::cerr.
        std::ostream* log_stream = &std::cerr;
        
        // Returns the stream where diagnostic information should be printed.
        // @param level the level the message belongs to
        // @return the stream, or nullptr
        std::ostream* log_at(VerbosityLevel level) const {
            return verbosity >= level ? log_stream : nullptr;
        }
    };

    // Result of the simplification process, containing the simplified polytope and
    // statistics.
    // @tparam Point the point type of the polytope
    template <typename Point>
    struct Result {
        // The simplified polytope.
        MetabolicPolytope<Point> P;

        // Number of finite bounds relaxed to +-infinity.
        unsigned bounds_relaxed = 0;

        // Number of dimensions fixed, i.e. tight box constraints converted
        // to equalities.
        unsigned dims_fixed = 0;      

        // Tracks if simplification was successful.
        bool success = true;
    };

    // A single side of the box bound, treated as a single row of the equivalent inequality
    // system A x <= b. An upper bound is of the form x_k <= b_u(k) and a lower bound is of
    // the form b_l(k) <= x_k.
    struct Ineq {
        // The index of the variable.
        unsigned k;

        // True for an upper bound, false for a lower bound.
        bool is_upper;

        // Orders inequalities so they can be added in an std::set.
        // @param other_ineq the inequality to compare against
        // @return true if this inequality precedes other_ineq
        bool operator<(Ineq const& other_ineq) const {
            if (k != other_ineq.k) return k < other_ineq.k;
            return is_upper < other_ineq.is_upper;
        }
    };

    // Evaluates the left hand side of the constraint a x <= b. Trivially
    // returns x_k or -x_k depending on the side of the inequality.
    // @param VT the vector type of x
    // @param a the inequality
    // @param the point to evaluate at
    // @return the inner dot product <a,x>
    template <typename VT>
    inline double get_row_value(Ineq const& a, VT const& x) {
        double xk = (double)x(a.k);
        return a.is_upper ? xk : -xk;
    }

    // Returns the right hand side of the box bound c written
    // as a row a x <= b. Notice that the lower bound b_l(k) <= x_k becomes -x_k
    // <= -b_l(k).
    // @tparam VT the vector type of the bounds
    // @param a the inequality
    // @param b the bound vector
    // @return the right hand side of the row
    template <typename VT>
    inline double get_row_rhs(Ineq const& a, VT const& b) {
        return a.is_upper ? (double)b(a.k) : -(double)b(a.k);
    }

    // Applies one bound of P to the highs model.
    // @tparam Point the point type of the polytope
    // @param highs the model
    // @param P the polytope
    // @param ineq the inequality to apply
    template <typename Point>
    void enforce_ineq(Highs & highs,
                      MetabolicPolytope<Point> const& P,
                      Ineq const& ineq)
    {
        typedef typename MetabolicPolytope<Point>::VT VT;
        const VT& b_l = P.getLowerBounds();
        const VT& b_u = P.getUpperBounds();

        // Stores the old lp state.
        double lo = highs.getLp().col_lower_[ineq.k];
        double hi = highs.getLp().col_upper_[ineq.k];

        // Adds only a single side inequality.
        if (ineq.is_upper) {
            hi = std::isinf((double)b_u(ineq.k)) ? kHighsInf : (double)b_u(ineq.k);
        } else {
            lo = std::isinf((double)b_l(ineq.k)) ? -kHighsInf : (double)b_l(ineq.k);
        }
        highs.changeColBounds((HighsInt)ineq.k, lo, hi);
    }
    
    // Applies the solver options shared by every lp in this file.
    // @param highs the model to configure
    // @param config the simplification configuration
    inline void configure_highs(Highs & highs, Config const & config) {
        highs.setOptionValue("output_flag", false);
        highs.setOptionValue("solver", "simplex");
        highs.setOptionValue("presolve", "off");
        highs.setOptionValue("simplex_strategy", 4);
        highs.setOptionValue("primal_feasibility_tolerance", config.primal_feasibility_tol);
        highs.setOptionValue("dual_feasibility_tolerance", config.dual_feasibility_tol);
        highs.setOptionValue("simplex_iteration_limit", config.simplex_iter_limit);
        highs.setOptionValue("time_limit", config.time_limit);
    }

    // Builds the LP that describes the feasible region of the Polytope.
    // @tparam Point the point type of the polytope
    // @param P the MetabolicPolytope
    // @param highs the highs model
    template <typename Point>
    void build_lp_model(MetabolicPolytope<Point> const& P, 
                        Highs& highs,
                        Config const& config
                        ) 
    {
        typedef typename MetabolicPolytope<Point>::MT MT;
        typedef typename MetabolicPolytope<Point>::VT VT;

        const MT& A_eq = P.getEqualities();
        const VT& b_u = P.getUpperBounds();
        const VT& b_l = P.getLowerBounds();
        const VT& b_eq = P.getEqualityBounds();
        unsigned d = P.getDimension();

        // Sets solver options for HiGHS
        configure_highs(highs, config);

        for (unsigned j = 0; j < d; ++j) {
            double low = std::isinf((double)b_l(j)) ? -kHighsInf : (double)b_l(j);
            double high = std::isinf((double)b_u(j)) ? kHighsInf : (double)b_u(j);
            highs.addVar(low, high);
        }
            
        for (unsigned i = 0; i < (unsigned)A_eq.rows(); ++i) {
            std::vector<HighsInt> indices;
            std::vector<double> values;
            for (typename MT::InnerIterator it(A_eq, i); it; ++it) {
                indices.push_back((HighsInt)it.col());
                values.push_back((double)it.value());
            }
            highs.addRow((double)b_eq(i), (double)b_eq(i), indices.size(), indices.data(), values.data());
        }
    }

    // Builds the output MetabolicPolytope after simplification by reading the
    // column bounds and equality rows directly from the HiGHS model.
    // @tparam Point the point type of the polytope
    // @param highs the HiGHS model after simplification
    // @return the simplified MetabolicPolytope
    template <typename Point>
    void highs_to_polytope(Highs const& highs, 
                           MetabolicPolytope<Point> & P) 
    {
        typedef typename MetabolicPolytope<Point>::MT MT;
        typedef typename MetabolicPolytope<Point>::VT VT;
        typedef typename MetabolicPolytope<Point>::NT NT;
        typedef typename MetabolicPolytope<Point>::Triplet Triplet;

        const NT INF = std::numeric_limits<NT>::infinity();

        HighsLp lp = highs.getLp();
        
        unsigned d = (unsigned)highs.getNumCol();

        // Grabs the lower/upper bounds of the reaction variables.
        VT b_l_new(d), b_u_new(d);
        for (unsigned j = 0; j < d; ++j) {
            b_l_new(j) = lp.col_lower_[j] <=  -kHighsInf ? -INF : (NT)lp.col_lower_[j];
            b_u_new(j) = lp.col_upper_[j] >= kHighsInf ? INF : (NT)lp.col_upper_[j];
        }

        // Builds A_eq and b_eq.
        lp.a_matrix_.ensureRowwise();
        const auto& start = lp.a_matrix_.start_;
        const auto& indices = lp.a_matrix_.index_;
        const auto& values = lp.a_matrix_.value_;

        unsigned row_count = (unsigned)lp.num_row_;
        std::vector<Triplet> triplets;
        VT b_eq_new(row_count);

        for (unsigned i = 0; i < row_count; ++i) {
            for (int p = start[i]; p < start[i+1]; ++p) {
                triplets.push_back(Triplet(i, (unsigned)indices[p], (NT)values[p]));
            }
            b_eq_new(i) = (NT)lp.row_lower_[i];
        }

        MT A_eq_new(row_count, d);
        A_eq_new.setFromTriplets(triplets.begin(), triplets.end());
        A_eq_new.makeCompressed();

        P =  MetabolicPolytope<Point>(
            d,
            A_eq_new,
            b_l_new,
            b_u_new,
            b_eq_new
        );
    }

    // Converts every variable that the inequalities and bounds pin to
    // a single value, into an equality.
    //
    // @tparam Point the point type of the polytope
    // @param highs the model, already built from P
    // @param config the simplification configuration
    // @return the polytope with the degenerate dimensions moved into A_eq
    template<typename Point>
    MetabolicPolytope<Point> fix_dimensions(Highs & highs,
                                            MetabolicPolytope<Point> const& P,
                                            Config const& config) 
    {   
        typedef typename MetabolicPolytope<Point>::MT MT;
        typedef typename MetabolicPolytope<Point>::VT VT;
        typedef typename MetabolicPolytope<Point>::NT NT;

        MetabolicPolytope<Point> tP = P;
        MT const& A_eq    = P.getEqualities();
        unsigned const& d = P.getDimension();
        unsigned const& m = (unsigned)A_eq.rows();

        std::vector<std::vector<unsigned>> var_apps(d);

        for (unsigned i = 0; i < m; ++i) 
            for (typename MT::InnerIterator it(A_eq, i); it; ++it)
                var_apps[(unsigned)it.col()].push_back(i);

        std::queue<unsigned> pending;
        std::vector<bool> queued(d, true);
        for (unsigned k = 0; k < d; ++k) 
            pending.push(k);

        auto enqueue = [&](unsigned k) {
            if (queued[k]) 
                return;
            queued[k] = true;
            pending.push(k);

        };

        auto requeue_neighbours = [&](unsigned k) {
            for (unsigned i : var_apps[k])
                for (typename MT::InnerIterator it(A_eq, i); it; ++it)
                    if ((unsigned)it.col() != k) enqueue((unsigned)it.col());
        };

        auto fix_dimension = [&](unsigned k, NT val) {
            HighsInt id = k;
            double coeff = 1.0;
            highs.addRow((double)val, (double)val, 1, &id, &coeff);
            highs.changeColBounds(k, -kHighsInf, kHighsInf);
        };
        
        std::vector<double> u_observed(d, std::numeric_limits<double>::infinity());
        std::vector<double> l_observed(d, -std::numeric_limits<double>::infinity());

        auto observe = [&](HighsSolution const& sol) {
            const auto& sol_vals = sol.col_value;
            for (unsigned j = 0; j < d; ++j) {
                if (sol_vals[j] > l_observed[j]) l_observed[j] = sol_vals[j];
                if (sol_vals[j] < u_observed[j]) u_observed[j] = sol_vals[j];
            }
        };

        auto observe_variation = [&](unsigned k) {
            return std::abs(u_observed[k]-l_observed[k]) > config.dim_tolerance;
        };

        highs.run();
        if (highs.getModelStatus() != HighsModelStatus::kOptimal) {
            log_diagnostics(config.log_at(VerbosityLevel::Summary),
                            "clarkson: dimension fixing failed with status ",
                            (int)highs.getModelStatus());

            return P;
        }
        observe(highs.getSolution());

        while (!pending.empty()) {
            unsigned k = pending.front();
            pending.pop();
            queued[k] = 0;

            double l = highs.getLp().col_lower_[k];
            double u = highs.getLp().col_upper_[k];

            if (l <= -kHighsInf && u >= kHighsInf) continue;

            if (l > -kHighsInf && u < kHighsInf && std::abs(u-l) < config.dim_tolerance) {
                fix_dimension(k, (NT)((l+u)/2.0));
                requeue_neighbours(k);
                continue;
            }

            if (var_apps[k].empty()) continue;

            if (observe_variation(k)) continue;

            highs.changeColCost((HighsInt)k, 1.0);
            highs.changeObjectiveSense(ObjSense::kMaximize);
            highs.run();
            if (highs.getModelStatus() != HighsModelStatus::kOptimal) {
                highs.changeColCost((HighsInt)k, 0.0);
                continue;
            }
            double max_val = highs.getObjectiveValue();
            observe(highs.getSolution());

            if (observe_variation(k)) {
                highs.changeColCost((HighsInt)k, 0.0);
                continue;
            }

            highs.changeObjectiveSense(ObjSense::kMinimize);
            highs.run();
            if (highs.getModelStatus() != HighsModelStatus::kOptimal) {
                highs.changeColCost((HighsInt)k, 0.0);
                continue;
            }
            double min_val = highs.getObjectiveValue();
            observe(highs.getSolution());

            highs.changeColCost((HighsInt)k, 0.0);

            if (std::abs(max_val - min_val) < config.dim_tolerance) {
                NT mid = (max_val+min_val)/NT(2);
                fix_dimension(k, mid);
                requeue_neighbours(k);
            }
        }
        highs_to_polytope(highs, tP);
        return tP;
    }
    
    // Finds a point in the interior of P by maximizing a uniform slack variable against all bounds.
    //
    // The LP solved is the following:
    //
    // max y s.t. A_eq x = b_eq, b_l+y <= x <= b_u-y, 0 <= y <= 1
    // @tparam Point the point of the polytope
    // @tparam ZT the vector type of z
    // @param config the simplification configuration
    // @param z set to the interior point on success
    // @param success set to true if an interior point was found
    template<typename Point, typename ZT>
    void find_interior_point(MetabolicPolytope<Point> const& P,
                             Config const& config,
                             ZT& z,
                             bool & success)
    {
        typedef typename MetabolicPolytope<Point>::MT MT;
        typedef typename MetabolicPolytope<Point>::VT VT;

        const MT& A_eq = P.getEqualities();
        const VT& b_eq = P.getEqualityBounds();
        const VT& b_l = P.getLowerBounds();
        const VT& b_u = P.getUpperBounds();
        unsigned d = P.getDimension();

        Highs highs;
        configure_highs(highs, config);

        for (unsigned j = 0; j < d; ++j) {
            highs.addVar(-kHighsInf, kHighsInf);
        }
        highs.addVar(0.0, 1.0);

        for (unsigned i = 0; i < (unsigned)A_eq.rows(); ++i) {
            std::vector<HighsInt> indices;
            std::vector<double> values;
            for (typename MT::InnerIterator it(A_eq, i); it; ++it) {
                indices.push_back((HighsInt)it.col());
                values.push_back((double)it.value());
            }
            highs.addRow((double)b_eq(i), (double)b_eq(i), indices.size(), indices.data(), values.data());
        }

        for (unsigned j = 0; j < d; ++j) {
            if (!std::isinf((double)b_l(j))) {
                HighsInt idx[2] = {(HighsInt)j, (HighsInt)d};
                double val[2] = {1.0, -1.0};
                highs.addRow((double)b_l(j), kHighsInf, 2, idx, val);
            }
            if (!std::isinf((double)b_u(j))) {
                HighsInt idx[2] = {(HighsInt)j, (HighsInt)d};
                double val[2] = {1.0, 1.0};
                highs.addRow(-kHighsInf, (double)b_u(j), 2, idx, val);
            }
        }
        highs.changeColCost(d, 1.0);
        highs.changeObjectiveSense(ObjSense::kMaximize);
        highs.run();

        if (highs.getModelStatus() != HighsModelStatus::kOptimal) {
            log_diagnostics(config.log_at(VerbosityLevel::Summary),
                            "clarkson: interior LP status ",
                            (int)highs.getModelStatus());

            success = false;
            return;
        }

        if (highs.getObjectiveValue() < config.interior_tolerance) {
            log_diagnostics(config.log_at(VerbosityLevel::Summary),
                            "clarkson: slack ", highs.getObjectiveValue(),
                            " not significant after dimension fixing");
        
            success = false;
            return;
        }

        log_diagnostics(config.log_at(VerbosityLevel::Summary),
                        "clarkson: interior LP optimal, slack objective = ",
                        highs.getObjectiveValue());

        const auto& sol = highs.getSolution().col_value;
        z.resize(d);
        for (unsigned j = 0; j < d; ++j)
            z(j) = (typename ZT::Scalar)sol[j];
        
        success = true;
    }

    
    // Shoots the ray z+t*r, t >= 0, and returns the first box bound it crosses.
    // @tparam Point the point type of the polytope
    // @tparam ZT the vector type of z
    // @param P the polytope
    // @param z the ray origin, a interior point of P
    // @param r the ray direction
    // @param config the simplification configuration
    // @param success false if the ray escapes without hiting a facet
    // @return the facet hit first, meaningful only when success is true
    template <typename Point, typename ZT>
    Ineq ray_shoot(MetabolicPolytope<Point> const& P,
                   ZT const& z,
                   ZT const& r,
                   Config const& config,
                   bool & success)
    {
        typedef typename MetabolicPolytope<Point>::VT VT;
        const VT& b_l = P.getLowerBounds();
        const VT& b_u = P.getUpperBounds();
        unsigned d = P.getDimension();

        double best = std::numeric_limits<double>::infinity();

        Ineq hit;
        bool found = false;

        // Goes over all variables.
        for (unsigned k = 0; k < d; ++k) {
            double rk = (double)r(k);
            if (std::abs(rk) < config.ray_tolerance) continue;

            // Goes over both inequalities.
            for (unsigned side = 0; side < 2; ++side) {
                Ineq c{k, side==1};

                double tr = c.is_upper ? rk : -rk;
                if (tr <= config.ray_tolerance) continue;

                double rhs = get_row_rhs(c, c.is_upper ? b_u : b_l);
                if (std::isinf(rhs)) continue;

                double tz = get_row_value(c, z);
                double t = (rhs-tz)/tr;
                if (t < 0.0) continue;

                if (!found || t < best) {
                    best = t;
                    hit = c;
                    found = true;
                }
            }
        }
        
        success = found;
        return hit;
    }

    // Tests whether the side ineq is redundant given the essential set I. The model
    // already carries I, so only the tested constraint is temporarily applied.
    //
    // The tested bound is relaxed by `config.relaxation` rather than removed, and
    // if the derived solution x* is feasible for the original LP, then the constraint
    // is marked as redundant.
    // @tparam Point the point type of the metabolic polytope.
    // @param highs the model (with `I` applied) 
    // @param P the polytope
    // @param ineq the constraint to be tested
    // @param config the simplification configuration
    // @param success a variable tracking if the LP failed
    // @return whether ineq is redundant, and the LP optimum
    template <typename Point>
    std::pair<bool, typename MetabolicPolytope<Point>::VT> test_redundancy(Highs & highs,
                                                                           MetabolicPolytope<Point> const& P,
                                                                           Ineq const& ineq,
                                                                           Config const& config,
                                                                           bool & success)
    {   
        typedef typename MetabolicPolytope<Point>::VT VT;
        const VT& b_l = P.getLowerBounds();
        const VT& b_u = P.getUpperBounds();
        unsigned d = P.getDimension();

        double old_u = highs.getLp().col_upper_[ineq.k];
        double old_l = highs.getLp().col_lower_[ineq.k];
        double u = ineq.is_upper ? (double)b_u(ineq.k)+config.relaxation_gap : old_u;
        double l = !ineq.is_upper ? (double)b_l(ineq.k)-config.relaxation_gap : old_l;

        highs.changeColBounds((HighsInt)ineq.k, l, u);
        highs.changeColCost((HighsInt)ineq.k, 1.0);
        highs.changeObjectiveSense(ineq.is_upper ? ObjSense::kMaximize : ObjSense::kMinimize);
        highs.run();

        HighsModelStatus st = highs.getModelStatus();

        if (st != HighsModelStatus::kOptimal) {
            log_diagnostics(config.log_at(VerbosityLevel::Detailed),
                            "clarkson: LP status ", (int)st,
                            " on coordinate ", ineq.k,
                            (ineq.is_upper ? " upper" : " lower"));

            highs.changeColBounds((HighsInt)ineq.k, old_l, old_u);
            highs.changeColCost((HighsInt)ineq.k, 0.0);
            success = false;
            return {false, VT(0)};
        }

        success = true;

        const auto& sol = highs.getSolution().col_value;
        VT x_star(d);
        for (unsigned j = 0; j < d; ++j)
            x_star(j) = (typename VT::Scalar)sol[j];

        double rhs = get_row_rhs(ineq, ineq.is_upper ? b_u : b_l);
        double opt = get_row_value(ineq, x_star);

        highs.changeColBounds((HighsInt)ineq.k, old_l, old_u);
        highs.changeColCost((HighsInt)ineq.k, 0.0);
        return {opt <= rhs+config.facet_tolerance, x_star};

    }

    // Clarkson decides the fate of a single constraint, returning
    // @tparam Point the point type of the metabolic polytope
    // @tparam ZT the vector type of z
    // @param highs the model (with `I` applied) 
    // @param P the polytope
    // @param z a point in the interior of the polytope
    // @param k_ineq the candidate constraint
    // @param config the clarkson configuration
    // @param success false if the LP or the ray shot failed
    template <typename Point, typename ZT>
    std::pair<bool, Ineq> clarkson(Highs & highs,
                        MetabolicPolytope<Point> const& P,
                        ZT const& z,
                        Ineq const& k_ineq,
                        Config const& config,
                        bool & success) 
    {
        typedef typename MetabolicPolytope<Point>::VT VT;
        unsigned d = P.getDimension();

        auto [is_redundant, x_star] = test_redundancy(
            highs, P, k_ineq, config, success
        );

        // Handle the case were the Lp solver failed.
        if (!success) {
            return {false, Ineq{}};
        }

        if (!is_redundant) {
            VT r = x_star-(VT)z;
            Ineq hit = ray_shoot(P, z, r, config, success);
            if (!success) return {false, Ineq{}};
            return {false, hit};
        } else {
            return {true, k_ineq};
        }
    }

    // Removes redundant inequalities from the representation using Clarkson's algorithm.
    //
    // The model starts with every inequality relaxed and gains them back one at a time
    // as they are proved essential, so every LP is solved against the essential set I 
    // found so far rather than the full set of inequalities, keeping the LP sizes at a minimum.
    // @tparam Point the point tyoe if the polytope
    // @tparam ZT the vector type of z
    // @param highs the model, built from P
    // @param P the polytope
    // @param z an interior point of P
    // @param config the simplification configuration
    // @return the polytope with redundant bounds relaxed
    template<typename Point, typename ZT>
    MetabolicPolytope<Point> redundancy_removal_clarkson(Highs & highs,
                                                         MetabolicPolytope<Point> const& P,
                                                         ZT const& z,
                                                         Config const& config)
    {
        typedef typename MetabolicPolytope<Point>::MT MT;
        typedef typename MetabolicPolytope<Point>::VT VT;
        typedef typename MetabolicPolytope<Point>::NT NT;

        const NT INF = std::numeric_limits<NT>::infinity();
        unsigned d = P.getDimension();
        const MT& A_eq = P.getEqualities();
        const VT& b_eq = P.getEqualityBounds();
        const VT& b_l = P.getLowerBounds();
        const VT& b_u = P.getUpperBounds();

        // Starts with all inequalities relaxed.
        for (unsigned j = 0; j < d; ++j)
            highs.changeColBounds((HighsInt)j, -kHighsInf, kHighsInf);

        // Holds the inequalities with unknown redundancy status.
        std::set<Ineq> J;
        for (unsigned k = 0; k < d; ++k) {
            if (!std::isinf((double)b_l(k))) J.insert(Ineq{k, false});
            if (!std::isinf((double)b_u(k))) J.insert(Ineq{k, true});
        }

        std::vector<Ineq> I;

        std::mt19937 rng(config.clarkson_seed);
        std::vector<unsigned> fail_count(d, 0);

        while (!J.empty()) {
            // Picks constraints at random to make progress when LPs fail.
            std::uniform_int_distribution<std::size_t> pick(0, J.size()-1);
            auto it = J.begin();
            std::advance(it, pick(rng));

            Ineq k_ineq = *it;

            bool success = false;
            auto [is_redundant, ineq] = clarkson(highs, P, z, k_ineq, config, success);

            if (!success) {
                log_diagnostics(config.log_at(VerbosityLevel::Detailed),
                                "clarkson: LP Failed on coordinate ", k_ineq.k,
                                (k_ineq.is_upper ? " upper" : " lower"));

                if (++fail_count[k_ineq.k] > config.failed_iter_count) {
                    for (Ineq const& j : J) {
                        I.push_back(j);
                    }
                    break;
                }
                continue;
            }

            if (is_redundant) {
                J.erase(k_ineq);
            } else {
                if (!J.erase(ineq)) {
                    log_diagnostics(config.log_at(VerbosityLevel::Detailed),
                                    "clarkson: hit already essential inequality ",
                                    ineq.k, (ineq.is_upper ? " upper" : " lower"),
                                    ", marking ", k_ineq.k,
                                    (k_ineq.is_upper ? " upper" : " lower"),
                                    " as essential instead");
                
                    I.push_back(k_ineq);
                    enforce_ineq(highs, P, k_ineq);
                    J.erase(k_ineq);
                } else {
                    I.push_back(ineq);
                    enforce_ineq(highs, P, ineq);
                }
            }
        }

        std::vector<bool> keep_lo(d, 0), keep_hi(d, 0);
        for (Ineq const& in : I) {
            if (in.is_upper) keep_hi[in.k] = true;
            else keep_lo[in.k] = true;
        }

        VT b_l_new(d), b_u_new(d);
        for (unsigned j = 0; j < d; ++j) {
            b_l_new(j) = keep_lo[j] ? b_l(j) : -INF;
            b_u_new(j) = keep_hi[j] ? b_u(j) : INF;
        }

        return MetabolicPolytope<Point>(d, A_eq, b_l_new, b_u_new, b_eq);
    }
    
    // Simplifies the polytope by fixing the degenerate dimensions and running
    // Clarkson's algorithm, which removes redundant constraints.
    // @tparam Point the point type of the polytope
    // @param config the simplification configuration
    // @return the simplified polytope, the counts, and whether it succeeded
    template<typename Point>
    Result<Point> simplify(MetabolicPolytope<Point> const& P,
                           Config const& config = Config{}) 
    {
        typedef typename MetabolicPolytope<Point>::VT VT;

        Highs highs;
        Result<Point> res;
        build_lp_model(P, highs, config);

        // Verifies the LP is not empty
        highs.run();
        if (highs.getModelStatus() != HighsModelStatus::kOptimal) {
            log_diagnostics(config.log_at(VerbosityLevel::Summary),
                            "clarkson: initial LP failed with status ",
                            (int)highs.getModelStatus());

            res.success = false;
            res.P = P;
            return res;
        }

        log_diagnostics(config.log_at(VerbosityLevel::Summary),
                        "clarkson: starting simplification on polytope with ",
                        P.getDimension(), " reactions, ", 
                        P.getNumEqualities(), " metabolites, and ",
                        P.getNumFiniteBounds(), " finite bounds");

        // Removes degenerate facets.
        if (config.fix_dimensions) {
            res.P = fix_dimensions(highs, P, config);
            res.dims_fixed = res.P.getNumEqualities()-P.getNumEqualities();

            log_diagnostics(config.log_at(VerbosityLevel::Summary),
                            "clarkson: fixed ", res.dims_fixed, " dimensions");
        } else {
            res.P = P;
            res.dims_fixed = 0;

            log_diagnostics(config.log_at(VerbosityLevel::Summary),
                            "clarkson: dimension fixing disabled");
        }


        // Looks for an interior point of P, such a point is essential for clarkson.
        VT z;
        bool success;
        find_interior_point(res.P, config, z, success);
  
        if (success) {
            res.P = redundancy_removal_clarkson(highs, res.P, z, config);
        } else {
            log_diagnostics(config.log_at(VerbosityLevel::Summary),
                            "clarkson: interior point not found, falling back to the exhaustive method");

            exhaustive_simplification::Config exhaustive_config;
            exhaustive_config.verbosity = config.verbosity;
            exhaustive_config.log_stream = config.log_stream;
            exhaustive_config.facet_tolerance = config.facet_tolerance;
            exhaustive_config.dim_tolerance = config.dim_tolerance;
            exhaustive_config.fix_dimensions = false;                   // Either fixed or the caller doesn't want dimension fixing.

            auto exhaustive_res = exhaustive_simplification::simplify(
                res.P, exhaustive_config
            );

            if (!exhaustive_res.success) {
                log_diagnostics(config.log_at(VerbosityLevel::Summary),
                                "clarkson: exhaustive simplification also failed, returning the dimension fixed polytope");

                res.success = false;
                return res;

            }
            res.P = exhaustive_res.P;
        }
        
        // Collects the statistics.
        res.bounds_relaxed = P.getNumFiniteBounds()-res.P.getNumFiniteBounds();
        res.success = true;

        log_diagnostics(config.log_at(VerbosityLevel::Summary),
                        "clarkson: simplification finished, fixed ", res.dims_fixed,
                        " dimensions, relaxed ", res.bounds_relaxed,
                        " bounds");

        return res;
    }
}
#endif