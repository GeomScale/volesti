// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2026 Vissarion Fisikopoulos
// Copyright (c) 2018-2026 Apostolos Chalkis
// Copyright (c) 2026      Dimitrios Pavlou

// Contributed and/or modified by Dimitrios Pavlou, as part of Google Summer of Code 2026 program

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef METABOLIC_SIMPLIFICATION_EXHAUSTIVE_HPP
#define METABOLIC_SIMPLIFICATION_EXHAUSTIVE_HPP

#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include "convex_bodies/metabolic_polytope.hpp"
#include "preprocess/metabolic/common.hpp"
#include "Highs.h"

// Configuration parameters controlling the simplification process.
struct ExhaustiveConfig {
    // Tolerance for marking a bound as redundant. A bound is relaxed if
    // moving it does not change the optimum by more than this quantity.
    double facet_tolerance = 1e-7;

    // Tolerance for marking a dimension as degenerate. A dimension is
    // fixed when its max and min differ less than this quantity.
    double dim_tolerance = 1e-7;

    // If true, degenerate dimensions are fixed.
    bool fix_dimensions = false;

    // The primal feasibility tolerance of the LP solver.
    double primal_feasibility_tol = 1e-7;

    // The dual feasibility tolerance of the LP solver.
    double dual_feasibility_tol = 1e-7;

    // Iteration limit for the simplex solves.
    double simplex_iter_limit = 1000;

    // Solver time limit.
    double time_limit = 200;
};

// Applies the solver options shared by every lp in this file.
// @param highs the model to configure
// @param config the simplification configuration
inline void configure_highs(Highs & highs, ExhaustiveConfig const & config) {
    highs.setOptionValue("output_flag", false);
    highs.setOptionValue("log_to_console", false);
    highs.setOptionValue("solver", "simplex");
    highs.setOptionValue("simplex_strategy", 4);
    highs.setOptionValue("primal_feasibility_tolerance", config.primal_feasibility_tol);
    highs.setOptionValue("dual_feasibility_tolerance", config.dual_feasibility_tol);
    highs.setOptionValue("simplex_iteration_limit", config.simplex_iter_limit);
    highs.setOptionValue("time_limit", config.time_limit);
}

// Simplifies a MetabolicPolytope by testing every box bound for redundancy and
// every dimension for degeneracy. A redundant bound is relaxed to +-infinity
// making the variable essentially free on that side, and a degenerate dimension
// is moved out of the box bounds into A_eq as a new equality row x_k = mid.
//
// Variables are visited iterativly until no further simplifications happen in a pass.
// @tparam Point the point type of the polytope
template <typename Point>
class ExhaustiveSimplifier {
    public:
        // Builds the LP model of the given polytope.
        // @param P_in the polytope to simplify
        // @param config_in the simplification configuration
        ExhaustiveSimplifier(MetabolicPolytope<Point> const& P_in, ExhaustiveConfig const& config_in
                            = ExhaustiveConfig{})
            : config(config_in), P(P_in)
        {
            build_lp_model(P, highs);
            configure_highs(highs, config);
        }

        // Runs the simplification.
        // @return the simplified polytope, and false if the polytope was empty
        std::pair<MetabolicPolytope<Point>, bool> simplify() 
        {               
            // Verifies the LP is not empty
            if (!run_lp()) return {P, false};

            bool simplified = false;
            while (!simplified) {
                simplified = true;

                for (unsigned k = 0; k < P.getDimension(); ++k) {
                    if (simplify_variable(k))
                        simplified = false;
                }
            }

            MetabolicPolytope<Point> Pnew;
            build_polytope_from_highs(highs, Pnew);
            return {Pnew, true};
        }
        
    private:
        // Number type. 
        typedef typename MetabolicPolytope<Point>::NT NT;

        // The configuration.
        ExhaustiveConfig config;

        // The polytope.
        MetabolicPolytope<Point> P;

        // The highs model.
        Highs highs;

        // Returns whether both bounds of a variable are infinite.
        // @param bl the lower bound
        // @param bu the upper bound
        // @bool true if the variable is free
        static bool is_free(double bl, double bu) {
            return bu >= kHighsInf && bl <= -kHighsInf;
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

        // Runs the lp currently stored in highs.
        // @return true if it solved to optimality
        bool run_lp() {
            highs.run();
            if (highs.getModelStatus() == HighsModelStatus::kOptimal) {
                return true;
            } else {
                return false;
            }
        }

        // Tests whether a bound is redundant, by relaxing it by +-1.
        // @param k the variable index
        // @param is_upper whether it is an upper bound
        // @param opt the optimum of the variable when the bound is unrelaxed
        // @return true if the bound can be relaxed
        bool bound_is_redundant(unsigned k, bool is_upper, NT opt) {
            double bl = highs.getLp().col_lower_[k];
            double bu = highs.getLp().col_upper_[k];

            if (is_upper ? !(bu < kHighsInf) : !(bl > -kHighsInf))
                return false;

            if (is_upper) {
                highs.changeObjectiveSense(ObjSense::kMaximize);
                highs.changeColBounds((HighsInt)k, bl, bu+1.0);
            } else {
                highs.changeObjectiveSense(ObjSense::kMinimize);
                highs.changeColBounds((HighsInt)k, bl-1.0, bu);
            }

            bool solved = run_lp();
            NT opt_relaxed = (NT)highs.getObjectiveValue();
            highs.changeColBounds((HighsInt)k, bl, bu);
            return solved && std::abs(opt_relaxed-opt) < config.facet_tolerance;
        }

        // Simplifies a single variable, fixing it if it is degenerate
        // and otherwise relaxing it if possible.
        // @param k the variable
        // @return true if the model was simplified
        bool simplify_variable(unsigned k) {
            double bl = highs.getLp().col_lower_[k];
            double bu = highs.getLp().col_upper_[k];

            if (is_free(bl, bu)) return false;

            highs.changeColCost(k, 1.0);

            // LP1: max with current bounds
            highs.changeObjectiveSense(ObjSense::kMaximize);
            if (!run_lp()) {
                highs.changeColCost(k, 0.0);
                return false;
            }

            NT max_val = (NT)highs.getObjectiveValue();
            bool upper_redundant = bound_is_redundant(k, true, max_val);

            // LP3: min with current bounds
            highs.changeObjectiveSense(ObjSense::kMinimize);
            if (!run_lp()) {
                highs.changeColCost(k, 0.0);
                return false;
            }

            NT min_val = (NT)highs.getObjectiveValue();
            bool lower_redundant = bound_is_redundant(k, false, min_val);

            highs.changeColCost(k, 0.0);

            bool tight = std::abs(max_val - min_val) < config.dim_tolerance;

            // Case1: Fix dimensions
            if (config.fix_dimensions && tight) {
                NT mid = (max_val+min_val)/NT(2);
                fix_dimension(k, mid);
                return true;
            }

            bool has_redundant_bound = false;
            // Case2: Relax bounds
            if (upper_redundant) {
                highs.changeColBounds(k, bl, kHighsInf);
                bu = kHighsInf;
                has_redundant_bound = true;
            } 
            if (lower_redundant && bl > -kHighsInf) {
                highs.changeColBounds(k, -kHighsInf, bu);
                has_redundant_bound = true;
            }

            return has_redundant_bound;
        }
};
#endif