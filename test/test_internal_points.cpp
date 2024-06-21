// VolEsti (volume computation and sampling library)

// Copyright (c) 2024 Vissarion Fisikopoulos
// Copyright (c) 2024 Apostolos Chalkis
// Copyright (c) 2024 Elias Tsigaridas

// Licensed under GNU LGPL.3, see LICENCE file

#include "doctest.h"
#include <fstream>
#include <iostream>

#include <boost/random.hpp>

#include "misc/misc.h"
#include "cartesian_geom/cartesian_kernel.h"
#include "convex_bodies/hpolytope.h"

#include "preprocess/max_inscribed_ball.hpp"
#include "preprocess/analytic_center_linear_ineq.h"

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
    NT max_min_eig_ratio = NT(2000);
    P = skinny_random_hpoly<Hpolytope, NT, PolyRNGType>(4, 180, pre_rounding, max_min_eig_ratio, 127);
    P.normalize();
    std::pair<Point, NT> InnerBall = P.ComputeInnerBall();

    NT tol = 1e-08;
    unsigned int maxiter = 500;
    auto [center, radius, converged] =  max_inscribed_ball(P.get_mat(), P.get_vec(), maxiter, tol);
    
    CHECK(P.is_in(Point(center)) == -1);
    CHECK(std::abs(radius - InnerBall.second) <= 1e-03);
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
    NT max_min_eig_ratio = NT(2000);
    P = skinny_random_hpoly<Hpolytope, NT, PolyRNGType>(50, 500, pre_rounding, max_min_eig_ratio, 127);
    P.normalize();

    bool feasibility_only = true; // compute only a feasible point
    NT tol = 1e-08;
    unsigned int maxiter = 500;
    auto [center, radius, converged] =  max_inscribed_ball(P.get_mat(), P.get_vec(), maxiter, tol, feasibility_only);

    CHECK(P.is_in(Point(center)) == -1);
    CHECK(converged);
}

template <typename NT>
void call_test_analytic_center() {
    typedef Cartesian <NT> Kernel;
    typedef typename Kernel::Point Point;
    typedef HPolytope <Point> Hpolytope;
    typedef typename Hpolytope::MT MT;
    typedef typename Hpolytope::VT VT;
    typedef boost::mt19937 PolyRNGType;
    Hpolytope P;

    std::cout << "\n--- Testing analytic center for skinny H-polytope" << std::endl;
    bool pre_rounding = true; // round random polytope before applying the skinny transformation 
    NT max_min_eig_ratio = NT(100);
    P = skinny_random_hpoly<Hpolytope, NT, PolyRNGType>(3, 15, pre_rounding, max_min_eig_ratio, 127);
    P.normalize();
    
    auto [Hessian, analytic_center, converged] = analytic_center_linear_ineq<MT, VT, NT>(P.get_mat(), P.get_vec());
    
    CHECK(P.is_in(Point(analytic_center)) == -1);
    CHECK(converged);
    CHECK(std::abs(analytic_center(0) + 4.75912) < 1e-04);
    CHECK(std::abs(analytic_center(1) + 4.28762) < 1e-04);
    CHECK(std::abs(analytic_center(2) - 7.54156) < 1e-04);
}

TEST_CASE("test_max_ball") {
    call_test_max_ball<double>();
}

TEST_CASE("test_feasibility_point") {
    call_test_max_ball_feasibility<double>();
}

TEST_CASE("test_analytic_center") {
    call_test_analytic_center<double>();
}
