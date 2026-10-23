// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2026 Vissarion Fisikopoulos
// Copyright (c) 2018-2026 Apostolos Chalkis
// Copyright (c) 2026      Dimitrios Pavlou

// Contributed and/or modified by Dimitrios Pavlou, as part of Google Summer of Code 2026 program

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef METABOLIC_EXHAUSTIVE_SIMPLIFICATION_HPP
#define METABOLIC_EXHAUSTIVE_SIMPLIFICATION_HPP

#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include "convex_bodies/metabolic_polytope.hpp"
#include "Highs.h"

namespace exhaustive_simplification {
    // How much diagnostic information to print during simplification.
    enum class VerbosityLevel {
        Silent = 0, // Prints nothing
        Summary,    // Prints a summary at the end of simplification
        Detailed    // Prints a summary at the end of simplification and
                    // information about each LP solved.
    };

    // Writes diagnostic information to the given stream if it is not null.
    // @tparam Args the types of the arguments to print
    // @param log_stream the stream to write to, or nullptr
    // @param args the arguments to print
    template <typename... Args>
    inline void log_diagnostics(std::ostream* log_stream, Args const&... args) {
        if (log_stream) {
            ((*log_stream << args), ...);
            *log_stream << std::endl;
        }
    }

    // Configuration parameters controlling the simplification process.
    struct Config {
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

        // How much information is printed during simplification.
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

        // Number of simplification passes performed
        unsigned passes = 0;

        // Tracks if simplification was successful.
        bool success = true;
    };

    // Applies the solver options shared by every lp in this file.
    // @param highs the model to configure
    // @param config the simplification configuration
    inline void configure_highs(Highs & highs, Config const & config) {
        highs.setOptionValue("output_flag", false);
        highs.setOptionValue("solver", "simplex");
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
                        Highs & highs,
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
    void build_simplified_polytope(Highs const& highs, 
                                   MetabolicPolytope<Point> & P) 
    {
        typedef typename MetabolicPolytope<Point>::MT MT;
        typedef typename MetabolicPolytope<Point>::VT VT;
        typedef typename MetabolicPolytope<Point>::NT NT;
        typedef typename MetabolicPolytope<Point>::Triplet Triplet;

        const NT INF = std::numeric_limits<NT>::infinity();

        HighsLp lp = highs.getLp();
        
        unsigned d = (unsigned)highs.getNumCol();

        // Grabs the lower/upper bounds of the reaction variables
        
        VT b_l_new(d), b_u_new(d);
        for (unsigned j = 0; j < d; ++j) {
            b_l_new(j) = lp.col_lower_[j] <=  -kHighsInf ? -INF : (NT)lp.col_lower_[j];
            b_u_new(j) = lp.col_upper_[j] >= kHighsInf ? INF : (NT)lp.col_upper_[j];
        }

        // Builds A_eq and b_eq

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

    // Simplifies the MetabolicPolytope by detecting redundant bounds and fixing degenerate
    // dimensions. Redundant bounds are relaxed to +-inf (making the variables free), and
    // fixed dimensions are moved out of the box bounds and into A_eq as new equality rows
    // x_k = mid. 
    //
    // The returned polytope will generally have more free variables in b_l/b_u and more rows
    // in A_eq than the input while describing the same region (up to numerical tolerance).
    // @tparam Point the point type of the polytope
    // @param P the input MetabolicPolytope
    // @param config the simplification config
    // @return the simplified MetabolicPolytope
    template <typename Point>
    Result<Point> simplify(MetabolicPolytope<Point> const& P, 
                           Config const& config = Config{}) 
    {               
        typedef typename MetabolicPolytope<Point>::NT NT;

        auto is_free = [](double bl, double bu){
            return bu >= kHighsInf && bl <= -kHighsInf;
        };

        auto fix_dimension = [](Highs & highs, unsigned k, NT val) {
            HighsInt id = k;
            double coeff = 1.0;
            highs.addRow((double)val, (double)val, 1, &id, &coeff);
            highs.changeColBounds(k, -kHighsInf, kHighsInf);
        };

        unsigned d = P.getDimension();

        Highs highs;
        build_lp_model(P, highs, config);

        // Verifies the LP is not empty
        highs.run();
        if (highs.getModelStatus() != HighsModelStatus::kOptimal) {
            log_diagnostics(config.log_at(VerbosityLevel::Summary),
                            "exhaustive: initial LP failed with status ",
                            (int)highs.getModelStatus());

            Result<Point> result;
            result.success = false;
            result.P = P;
            return result;
        }
        
        log_diagnostics(config.log_at(VerbosityLevel::Summary),
                        "exhaustive: starting simplification on polytope with ",
                        d, " reactions, ", 
                        P.getNumEqualities(), " metabolites, and ",
                        P.getNumFiniteBounds(), " finite bounds");

        bool simplified = false;

        unsigned pass = 0;
        while (!simplified) {
            simplified = true;

            unsigned fixed_this_pass = 0;
            unsigned relaxed_this_pass = 0;

            for (unsigned k = 0; k < d; ++k) {
                HighsLp const& lp = highs.getLp();
                double bl = lp.col_lower_[k];
                double bu = lp.col_upper_[k];

                if (is_free(bl, bu)) continue;     // Skips removed constraints or equalities

                highs.changeColCost(k, 1.0);

                // Runs the LP and documents failures
                auto run_lp = [&](char const* lp_type) {
                    highs.run();
                    if (highs.getModelStatus() == HighsModelStatus::kOptimal) 
                        return true;

                    log_diagnostics(config.log_at(VerbosityLevel::Detailed),
                                "exhaustive: LP ", lp_type, 
                                " failed on variable ", k,
                                " with status ", (int)highs.getModelStatus());
                    
                    return false;
                };

                // LP1: max with current bounds
                highs.changeObjectiveSense(ObjSense::kMaximize);
                if (!run_lp("max")) {
                    highs.changeColCost(k, 0.0);
                    continue;
                }

                NT max_val = (NT)highs.getObjectiveValue();

                // LP2: max with relaxed bound
                // Skipped if there is no upper bound
                bool upper_redundant = false;
                if (bu < kHighsInf) {
                    highs.changeColBounds(k, bl, bu+1.0);
                    bool lp_solved = run_lp("relaxed max");
                    NT max_val_relaxed = (NT)highs.getObjectiveValue();
                    highs.changeColBounds(k, bl, bu);
                    upper_redundant = lp_solved && std::abs(max_val_relaxed-max_val) < config.facet_tolerance;
                }

                // LP3: min with current bounds
                highs.changeObjectiveSense(ObjSense::kMinimize);
                if (!run_lp("min")) {
                    highs.changeColCost(k, 0.0);
                    continue;
                }
                NT min_val = (NT)highs.getObjectiveValue();

                // LP4: min with relaxed bound
                // Skipped if there is no lower bound
                bool lower_redundant = false;
                if (bl > -kHighsInf) {
                    highs.changeColBounds(k, bl-1.0, bu);
                    bool lp_solved = run_lp("relaxed min");
                    NT min_val_relaxed = (NT)highs.getObjectiveValue();
                    highs.changeColBounds(k, bl, bu);
                    lower_redundant = lp_solved && std::abs(min_val_relaxed-min_val) < config.facet_tolerance;
                }

                highs.changeColCost(k, 0.0);

                bool tight = std::abs(max_val - min_val) < config.dim_tolerance;

                // Case1: Fix dimensions
                if (config.fix_dimensions && tight) {
                    NT mid = (max_val+min_val)/NT(2);
                    fix_dimension(highs, k, mid);
                    simplified = false;
                    ++fixed_this_pass;
                    relaxed_this_pass += 2;
                    continue;
                }

                // Case2: Relax bounds
                if (upper_redundant && bu < kHighsInf) {
                    highs.changeColBounds(k, bl, kHighsInf);
                    bu = highs.getLp().col_upper_[k];
                    simplified = false;
                    ++relaxed_this_pass;
                } 
                if (lower_redundant && bl > -kHighsInf) {
                    highs.changeColBounds(k, -kHighsInf, bu);
                    simplified = false;
                    ++relaxed_this_pass;
                }
            }
            ++pass;

            log_diagnostics(config.log_at(VerbosityLevel::Detailed),
                            "exhaustive: pass ", pass,
                            " relaxed ", relaxed_this_pass,
                            " bounds, fixed ", fixed_this_pass,
                            " dimensions");
        }


        Result<Point> result;
        MetabolicPolytope<Point> Pnew;
        build_simplified_polytope(highs, Pnew);
        result.dims_fixed = Pnew.getNumEqualities()-P.getNumEqualities();
        result.bounds_relaxed = P.getNumFiniteBounds()-Pnew.getNumFiniteBounds();
        result.passes = pass;
        result.P = Pnew;

        log_diagnostics(config.log_at(VerbosityLevel::Summary),
                        "exhaustive: simplification finished in ", pass,
                        " passes, fixed ", result.dims_fixed,
                        " dimensions, relaxed ", result.bounds_relaxed,
                        " bounds");

        return result;
    }
}
#endif