// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2026 Vissarion Fisikopoulos
// Copyright (c) 2018-2026 Apostolos Chalkis
// Copyright (c) 2026      Dimitrios Pavlou

// Contributed and/or modified by Dimitrios Pavlou, as part of Google Summer of Code 2026 program

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef METABOLIC_SIMPLIFY_AND_TRANSFORM_HPP
#define METABOLIC_SIMPLIFY_AND_TRANSFORM_HPP

#include <iostream>
#include <utility>
#include <tuple>
#include <Eigen/Eigen>
#include "convex_bodies/metabolic_polytope.hpp"
#include "convex_bodies/hpolytope.h"
#include "preprocess/metabolic/scaling.hpp"
#include "preprocess/metabolic/transformation.hpp"

// Simplifies a metabolic network and transforms what remains into a
// full dimensional H-polytope. The network is scaled, simplified with the
// given method, and the result is transformed into an H-polytope.
// @tparam Point the point type of the polytope
// @tparam Simplifier the simplifier type, e.g. ClarksonSimplifier<Point>
// @tparam SimplifierConfig the configuration type of the simplifier
// @tparam ScalingPolicy the scaling policy e.g. MaxBoundScaling
// @param P the metabolic network to simplify and transform
// @param config the configuration of the simplifier
// @param scaling_policy the scaling policy to be used
// @return a tuple (HP, shift, N, simplified):
//         - HPolytope   : the full-dimensional polytope
//         - shift       : a particular solution of the A_eq x = b_eq
//         - N           : a basis of the nullspace of A_eq
//         - simplified  : false if the simplification failed
template <typename Point,
          typename Simplifier,
          typename SimplifierConfig,
          typename ScalingPolicy = MaxBoundScaling>
std::tuple<HPolytope<Point>,
           typename MetabolicPolytope<Point>::VT,
           Eigen::Matrix<typename MetabolicPolytope<Point>::NT, Eigen::Dynamic, Eigen::Dynamic>,
           bool> 
simplify_and_transform(
          MetabolicPolytope<Point> const& P,
          SimplifierConfig const& config = SimplifierConfig{},
          ScalingPolicy scaling_policy = ScalingPolicy{})
{
    Scaling<Point> s;
    MetabolicPolytope<Point> Ps = scale(P, s, scaling_policy);

    Simplifier simplifier(Ps, config);
    auto [Ps_simplified, simplified] = simplifier.simplify();

    if (!simplified) {
        #ifdef VOLESTI_DEBUG
        std::cout << "simplification failed in simplify_and_transform" << std::endl;
        #endif
        return {HPolytope<Point>{}, {}, {}, false};
    }
    
    pad_rows(s, Ps_simplified.getNumEqualities());

    MetabolicPolytope<Point> P_final = rescale(Ps_simplified, s, true);

    auto [HP, shift, N] = transform(P_final);
    return {HP, shift, N, simplified};
}
#endif