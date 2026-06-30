// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2026 Vissarion Fisikopoulos
// Copyright (c) 2018-2026 Apostolos Chalkis
// Copyright (c) 2026      Dimitrios Pavlou

// Contributed and/or modified by Dimitrios Pavlou, as part of Google Summer of Code 2026 program

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef SIMPLIFICATION_WARM_START_HPP
#define SIMPLIFICATION_WARM_START_HPP

#include "simplification/common.hpp"
#include "convex_bodies/metabolic_polytope.h"
#include "Highs.h"
#include <vector>
#include <cmath>
#include <limits>

namespace simplification {
    // Builds the LP that describes the feasible region of the Polytope.
    // @tparam Point the point type of the polytope
    // @param P the MetabolicPolytope
    // @param highs the highs model
    template <typename Point>
    void build_lp_model(MetabolicPolytope<Point> const& P, 
                        Highs& highs) 
    {
        typedef typename MetabolicPolytope<Point>::MT MT;
        typedef typename MetabolicPolytope<Point>::VT VT;

        const MT& A_eq = P.getEqualities();
        const VT& b_u = P.getUpperBounds();
        const VT& b_l = P.getLowerBounds();
        const VT& b_eq = P.getEqualityBounds();
        unsigned d = P.getDimension();

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

        auto is_free = [](double bl, double bu){return std::isinf(bu) && std::isinf(bl);};
        auto fix_dimension = [](Highs & highs, unsigned k, NT val) {
            HighsInt id = k;
            double coeff = 1.0;
            highs.addRow((double)val, (double)val, 1, &id, &coeff);
            highs.changeColBounds(k, -kHighsInf, kHighsInf);
        };

        unsigned d = P.getDimension();

        Highs highs;
        highs.setOptionValue("output_flag", false);
        highs.setOptionValue("solver", "simplex");
        highs.setOptionValue("simplex_strategy", 4);
        build_lp_model(P, highs);

        // Verifies the LP is not empty
        highs.run();
        if (highs.getModelStatus() != HighsModelStatus::kOptimal) {
            Result<Point> result;
            result.success = false;
            result.P = P;
            return result;
        }
        
        bool simplified = false;
        while (!simplified) {
            simplified = true;
            for (unsigned k = 0; k < d; ++k) {
                double bl = highs.getLp().col_lower_[k];
                double bu = highs.getLp().col_upper_[k];

                if (is_free(bl, bu)) continue;     // Skips removed constraints or equalities

                highs.changeColCost(k, 1.0);

                // LP1: max with current bounds
                highs.changeObjectiveSense(ObjSense::kMaximize);
                highs.run();
                NT max_val = (NT)highs.getObjectiveValue();

                // LP2: max with relaxed bound
                // Skipped if there is no upper bound
                bool upper_redundant = false;
                if (!std::isinf(bu)) {
                    highs.changeColBounds(k, bl, bu+1.0);
                    highs.run();
                    NT max_val_relaxed = (NT)highs.getObjectiveValue();
                    highs.changeColBounds(k, bl, bu);
                    upper_redundant = std::abs(max_val_relaxed-max_val) < config.facet_tolerance;
                }

                // LP3: min with current bounds
                highs.changeObjectiveSense(ObjSense::kMinimize);
                highs.run();
                NT min_val = (NT)highs.getObjectiveValue();

                // LP4: min with relaxed bound
                // Skipped if there is no lower bound
                bool lower_redundant = false;
                if (!std::isinf(bl)) {
                    highs.changeColBounds(k, bl-1.0, bu);
                    highs.run();
                    NT min_val_relaxed = (NT)highs.getObjectiveValue();
                    highs.changeColBounds(k, bl, bu);
                    lower_redundant = std::abs(min_val_relaxed-min_val) < config.facet_tolerance;
                }

                highs.changeColCost(k, 0.0);

                bool tight = std::abs(max_val - min_val) < config.dim_tolerance;

                // Case1: Fix dimensions
                if (config.fix_dimensions && tight) {
                    NT mid = (max_val+min_val)/NT(2);
                    fix_dimension(highs, k, mid);
                    simplified = false;
                    continue;
                }

                // Case2: Relax bounds
                if (upper_redundant && !std::isinf(bu)) {
                    highs.changeColBounds(k, bl, kHighsInf);
                    bu = highs.getLp().col_upper_[k];
                    simplified = false;
                } 
                if (lower_redundant && !std::isinf(bl)) {
                    highs.changeColBounds(k, -kHighsInf, bu);
                    simplified = false;
                }
            }

        }


        Result<Point> result;
        MetabolicPolytope<Point> Pnew;
        build_simplified_polytope(highs, Pnew);
        result.dims_fixed = Pnew.getNumEqualities()-P.getNumEqualities();
        result.bounds_relaxed = P.getNumFiniteBounds()-Pnew.getNumFiniteBounds();
        result.P = Pnew;

        return result;
    }
}
#endif