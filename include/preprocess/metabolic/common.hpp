
// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2026 Vissarion Fisikopoulos
// Copyright (c) 2018-2026 Apostolos Chalkis
// Copyright (c) 2026      Dimitrios Pavlou

// Contributed and/or modified by Dimitrios Pavlou, as part of Google Summer of Code 2026 program

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef METABOLIC_COMMON_HPP
#define METABOLIC_COMMON_HPP

#include "convex_bodies/metabolic_polytope.hpp"
#include "Highs.h"

// Builds the LP that describes the feasible region of the Polytope.
// @tparam Point the point type of the polytope
// @param P the MetabolicPolytope
// @param highs the highs model
template <typename Point>
void build_lp_model(MetabolicPolytope<Point> const& P, 
                    Highs & highs) 
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

// Builds the output MetabolicPolytope by reading the column bounds 
// and equality rows directly from the HiGHS model.
// @tparam Point the point type of the polytope
// @param highs the HiGHS model after simplification
// @return the simplified MetabolicPolytope
template <typename Point>
void build_polytope_from_highs(Highs const& highs, 
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
#endif