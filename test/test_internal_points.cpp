// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2024 Vissarion Fisikopoulos
// Copyright (c) 2018-2024 Apostolos Chalkis

// Contributed by Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

#include "doctest.h"
#include <fstream>
#include <iostream>

#include <boost/random.hpp>

#include "misc/misc.h"
#include "cartesian_geom/cartesian_kernel.h"
#include "convex_bodies/hpolytope.h"

#include "preprocess/max_inscribed_ball.hpp"

#include "generators/known_polytope_generators.h"
#include "generators/h_polytopes_generator.h"

template <typename NT>
void call_test_max_ball() {
    typedef Cartesian <NT> Kernel;
    typedef typename Kernel::Point Point;
    typedef HPolytope <Point> Hpolytope;
    typedef boost::mt19937 PolyRNGType;
    Hpolytope P;

    std::cout << "\n--- Testing Chebychev ball for skinny H-polytope" << std::endl;
    bool pre_rounding = true; // round random polytope before applying the skinny transformation
    P = skinny_random_hpoly<Hpolytope, NT, PolyRNGType>(7, 200, pre_rounding, 2000.0, 14027);
    P.normalize();
    std::pair<Point, NT> InnerBall = P.ComputeInnerBall();

    auto [center, radius, converged] =  max_inscribed_ball(P.get_mat(), P.get_vec(), 500, 1e-08);
    
    CHECK(P.is_in(Point(center)) == -1);
    CHECK(std::abs(radius - InnerBall.second) <= 1e-06);
    CHECK(converged);
}

template <typename NT>
void call_test_max_ball_feasibility() {
    typedef Cartesian <NT> Kernel;
    typedef typename Kernel::Point Point;
    typedef HPolytope <Point> Hpolytope;
    typedef boost::mt19937 PolyRNGType;
    Hpolytope P;

    std::cout << "\n--- Testing feasibility point for skinny H-polytope" << std::endl;
    bool pre_rounding = true; // round random polytope before applying the skinny transformation 
    P = skinny_random_hpoly<Hpolytope, NT, PolyRNGType>(50, 500, pre_rounding, 2000.0, 127);
    P.normalize();

    bool feasibility_only = true; // compute only a feasible point
    NT tol = 1e-08;
    unsigned int maxiter = 500;
    auto [center, radius, converged] =  max_inscribed_ball(P.get_mat(), P.get_vec(), maxiter, tol, feasibility_only);

    CHECK(P.is_in(Point(center)) == -1);
    CHECK(converged);
}

TEST_CASE("test_max_ball") {
    call_test_max_ball<double>();
    call_test_max_ball_feasibility<double>();
}