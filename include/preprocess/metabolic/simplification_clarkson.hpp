// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2026 Vissarion Fisikopoulos
// Copyright (c) 2018-2026 Apostolos Chalkis
// Copyright (c) 2026      Dimitrios Pavlou

// Contributed and/or modified by Dimitrios Pavlou, as part of Google Summer of Code 2026 program

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef METABOLIC_SIMPLIFICATION_CLARKSON_HPP
#define METABOLIC_SIMPLIFICATION_CLARKSON_HPP

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
#include "preprocess/metabolic/simplification_exhaustive.hpp"
#include "common.hpp"
#include "Highs.h"

// Configuration parameters controlling the simplification process.
struct ClarksonConfig : ExhaustiveConfig {
    // The error tolerance for the interior point.
    double interior_tolerance = 1e-9;
    
    // The error tolerance for the ray shooting stage of clarkson.
    double ray_tolerance = 1e-9;

    // The gap by which a bound is relaxed in the redundancy LP.
    double relaxation_gap = 1.0;

    // The bound on the number of failed iteration's in clarkson.
    unsigned failed_iter_count = 50;

    // The seed based on which clarkson picks inequalities.
    unsigned clarkson_seed = 0;
};

// Simplifies a MetabolicPolytope using Clarkson's algorithm. The model starts with
// every bound relaxed and iteratively finds and applies essential constraints to the model.
//
// Clarkson needs a point in the interior of the polytope, which only exists once
// degenerate dimensions have been fixed. If no interior point is found the exhaustive
// method is used as fallback.
// @tparam Point the point type of the polytope
template <typename Point>
class ClarksonSimplifier {
    public:
        // Builds the LP model of the given polytope.
        // @param P_in the polytope to simplify
        // @param config_in the simplification configuration
        ClarksonSimplifier(MetabolicPolytope<Point> const& P_in,
                        ClarksonConfig const& config_in = ClarksonConfig{})
            : config(config_in), P(P_in)
        {
            build_lp_model(P, highs);
            configure_highs(highs, config);
        }

        // Runs the simplification.
        // @return the simplified polytope, and false if the polytope was empty
        std::pair<MetabolicPolytope<Point>, bool> simplify() {
            MetabolicPolytope<Point> Ps = P;

            if (config.fix_dimensions) {
                Ps = fix_degenerate_dimensions();
                P = Ps;
            }
            
            if (find_interior_point()) {
                Ps = redundancy_removal_clarkson();
            } else {
                // Without an interior point clarkson cannot run.
                ExhaustiveConfig exhaustive_config = config;
                exhaustive_config.fix_dimensions = false;

                ExhaustiveSimplifier<Point> fallback(P, exhaustive_config);
                auto [Pnew, ok] = fallback.simplify();
                if (!ok) return {P, false};
                Ps = Pnew;
            }

            return {Ps, true};
        }

    private:
        // Types.
        typedef typename MetabolicPolytope<Point>::MT MT;
        typedef typename MetabolicPolytope<Point>::VT VT;
        typedef typename MetabolicPolytope<Point>::NT NT;

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

        // The simplification configuration.
        ClarksonConfig config;

        // The input polytope.
        MetabolicPolytope<Point> P;

        // The LP model.
        Highs highs;

        // The interior point.
        VT z;
        // Runs the LP stored in highs.
        // @return true if solved to optimality
        bool run_lp() {
            highs.run();
            return highs.getModelStatus() == HighsModelStatus::kOptimal;
        }

        // Evaluates the left hand side of the constraint a x <= b. Trivially
        // returns x_k or -x_k depending on the side of the inequality.
        // @param a the inequality
        // @param the point to evaluate at
        // @return the inner dot product <a,x>
        static double row_value(Ineq const& a, VT const& x) {
            double xk = (double)x(a.k);
            return a.is_upper ? xk : -xk;
        }

        // Returns the right hand side of the box bound c written
        // as a row a x <= b. Notice that the lower bound b_l(k) <= x_k becomes -x_k
        // <= -b_l(k).
        // @param a the inequality
        // @param b the bound vector
        // @return the right hand side of the row
        inline double row_rhs(Ineq const& a) {
            return a.is_upper ? (double)P.getUpperBounds()(a.k) : 
                               -(double)P.getLowerBounds()(a.k);
        }

        // Collapses the bounds of a variable, by adding x_k = val to A_eq.
        // @param k the variable index
        // @param val the value to fix it at
        void fix_dimension(unsigned k, NT val) {
            HighsInt id = k;
            double coeff = 1.0;
            highs.addRow((double)val, (double)val, 1, &id, &coeff);
            highs.changeColBounds(k, -kHighsInf, kHighsInf);
        }

        // Applies one bound of P to the highs model.
        // @param ineq the inequality to apply
        void enforce_ineq(Ineq const& ineq)
        {
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

        // Shoots the ray z+t*r, t >= 0, and returns the first box bound it crosses.
        // @param r the ray direction
        // @param config the simplification configuration
        // @param success false if the ray escapes without hiting a facet
        // @return the facet hit first, meaningful only when success is true
        bool ray_shoot(VT const& r, Ineq & hit)
        {
            unsigned d = P.getDimension();
            double best = std::numeric_limits<double>::infinity();
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

                    double rhs = row_rhs(c);
                    if (std::isinf(rhs)) continue;

                    double tz = row_value(c, z);
                    double t = (rhs-tz)/tr;
                    if (t < 0.0) continue;

                    if (!found || t < best) {
                        best = t;
                        hit = c;
                        found = true;
                    }
                }
            }
            
            return found;
        }

        // Tests whether the side ineq is redundant given the essential set I. The model
        // already carries I, so only the tested constraint is temporarily applied.
        //
        // The tested bound is relaxed by `config.relaxation` rather than removed, and
        // if the derived solution x* is feasible for the original LP, then the constraint
        // is marked as redundant.
        // @param ineq the constraint to be tested
        // @param solved a variable tracking if the LP failed
        // @return whether ineq is redundant, and the LP optimum
        std::pair<bool, VT> test_redundancy(Ineq const& ineq, bool & solved)
        {   
            const VT& b_u = P.getUpperBounds();
            const VT& b_l = P.getLowerBounds();
            unsigned d = P.getDimension();

            double old_u = highs.getLp().col_upper_[ineq.k];
            double old_l = highs.getLp().col_lower_[ineq.k];
            double u = ineq.is_upper ? (double)b_u(ineq.k)+config.relaxation_gap : old_u;
            double l = !ineq.is_upper ? (double)b_l(ineq.k)-config.relaxation_gap : old_l;

            highs.changeColBounds((HighsInt)ineq.k, l, u);
            highs.changeColCost((HighsInt)ineq.k, 1.0);
            highs.changeObjectiveSense(ineq.is_upper ? ObjSense::kMaximize : ObjSense::kMinimize);
            
            solved = run_lp();

            VT x_star(d);
            bool redundant = false;
            if (solved) {
                const auto& sol = highs.getSolution().col_value;
                for (unsigned j = 0; j < d; ++j)
                    x_star(j) = (typename VT::Scalar)sol[j];

                redundant = row_value(ineq, x_star) <= row_rhs(ineq)+config.facet_tolerance;
            }

            highs.changeColBounds((HighsInt)ineq.k, old_l, old_u);
            highs.changeColCost((HighsInt)ineq.k, 0.0);

            return {redundant, x_star};
        }

        // Converts every variable that the inequalities and bounds pin to
        // a single value, into an equality.
        // @return the polytope with the degenerate dimensions moved into A_eq
        MetabolicPolytope<Point> fix_degenerate_dimensions() 
        {   
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


            std::vector<double> u_observed(d, std::numeric_limits<double>::infinity());
            std::vector<double> l_observed(d, -std::numeric_limits<double>::infinity());

            auto observe = [&]() {
                const auto& sol = highs.getSolution().col_value;
                for (unsigned j = 0; j < d; ++j) {
                    if (sol[j] > l_observed[j]) l_observed[j] = sol[j];
                    if (sol[j] < u_observed[j]) u_observed[j] = sol[j];
                }
            };

            auto observe_variation = [&](unsigned k) {
                return std::abs(u_observed[k]-l_observed[k]) > config.dim_tolerance;
            };

            if (!run_lp()) return tP;
            observe();

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

                if (!run_lp()) {
                    highs.changeColCost((HighsInt)k, 0.0);
                    continue;
                }
                double max_val = highs.getObjectiveValue();
                observe();

                if (observe_variation(k)) {
                    highs.changeColCost((HighsInt)k, 0.0);
                    continue;
                }

                highs.changeObjectiveSense(ObjSense::kMinimize);
                if (!run_lp()) {
                    highs.changeColCost((HighsInt)k, 0.0);
                    continue;
                }
                double min_val = highs.getObjectiveValue();
                observe();

                highs.changeColCost((HighsInt)k, 0.0);

                if (std::abs(max_val - min_val) < config.dim_tolerance) {
                    fix_dimension(k, (max_val+min_val)/NT(2));
                    requeue_neighbours(k);
                }
            }
            build_polytope_from_highs(highs, tP);
            return tP;
        }


    // Finds a point in the interior of P by maximizing a uniform slack variable against all bounds.
    //
    // The LP solved is the following:
    //
    // max y s.t. A_eq x = b_eq, b_l+y <= x <= b_u-y, 0 <= y <= 1
    // Note: the point is stored in z
    // @return true if a point with significant slack was found
    bool find_interior_point()
    {
        const MT& A_eq = P.getEqualities();
        const VT& b_eq = P.getEqualityBounds();
        const VT& b_l = P.getLowerBounds();
        const VT& b_u = P.getUpperBounds();
        unsigned d = P.getDimension();

        Highs slack_highs;
        configure_highs(slack_highs, config);

        for (unsigned j = 0; j < d; ++j) {
            slack_highs.addVar(-kHighsInf, kHighsInf);
        }
        slack_highs.addVar(0.0, 1.0);

        for (unsigned i = 0; i < (unsigned)A_eq.rows(); ++i) {
            std::vector<HighsInt> indices;
            std::vector<double> values;
            for (typename MT::InnerIterator it(A_eq, i); it; ++it) {
                indices.push_back((HighsInt)it.col());
                values.push_back((double)it.value());
            }
            slack_highs.addRow((double)b_eq(i), (double)b_eq(i), indices.size(), indices.data(), values.data());
        }

        for (unsigned j = 0; j < d; ++j) {
            if (!std::isinf((double)b_l(j))) {
                HighsInt idx[2] = {(HighsInt)j, (HighsInt)d};
                double val[2] = {1.0, -1.0};
                slack_highs.addRow((double)b_l(j), kHighsInf, 2, idx, val);
            }
            if (!std::isinf((double)b_u(j))) {
                HighsInt idx[2] = {(HighsInt)j, (HighsInt)d};
                double val[2] = {1.0, 1.0};
                slack_highs.addRow(-kHighsInf, (double)b_u(j), 2, idx, val);
            }
        }
        slack_highs.changeColCost(d, 1.0);
        slack_highs.changeObjectiveSense(ObjSense::kMaximize);
        slack_highs.run();

        if (slack_highs.getModelStatus() != HighsModelStatus::kOptimal ||
            slack_highs.getObjectiveValue() < config.interior_tolerance)
            return false;

        const auto& sol = slack_highs.getSolution().col_value;
        z.resize(d);
        for (unsigned j = 0; j < d; ++j)
            z(j) = (typename VT::Scalar)sol[j];
        
        return true;
    }

    // Removes redundant inequalities from the representation using Clarkson's algorithm.
    //
    // The model starts with every inequality relaxed and gains them back one at a time
    // as they are proved essential, so every LP is solved against the essential set I 
    // found so far rather than the full set of inequalities, keeping the LP sizes at a minimum.
    // @return the simplified polytope
    MetabolicPolytope<Point> redundancy_removal_clarkson()
    {
        const NT INF = std::numeric_limits<NT>::infinity();
        unsigned d = P.getDimension();
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

            bool solved = false;
            auto [is_redundant, x_star] = test_redundancy(k_ineq, solved);

            if (!solved) {
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
                continue;
            }

            Ineq hit;
            if (!ray_shoot(x_star-z, hit) || !J.erase(hit)) {
                I.push_back(k_ineq);
                enforce_ineq(k_ineq);
                J.erase(k_ineq);
            } else {
                I.push_back(hit);
                enforce_ineq(hit);
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

        return MetabolicPolytope<Point>(d, P.getEqualities(), b_l_new, b_u_new, 
                                        P.getEqualityBounds());
    }
};

#endif