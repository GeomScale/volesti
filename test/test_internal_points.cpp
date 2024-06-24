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

#include "convex_bodies/orderpolytope.h"
#include "misc/poset.h"

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
void call_test_max_ball_sparse() {
    typedef Cartesian <NT> Kernel;
    typedef typename Kernel::Point Point;
    typedef HPolytope <Point> Hpolytope;
    typedef boost::mt19937 PolyRNGType;
    typedef typename OrderPolytope<Point>::VT VT;
    typedef typename OrderPolytope<Point>::MT MT;
    typedef typename Poset::RT RT;
    typedef typename Poset::RV RV;
    typedef Eigen::SparseMatrix<NT> SpMT;

    // Create Poset, 4 elements, a0 <= a1, a0 <= a2, a1 <= a3
    RV poset_data{{0, 1}, {0, 2}, {1, 3}};
    Poset poset(4, poset_data);
    
    // Initialize order polytope from the poset
    OrderPolytope<Point> OP(poset);
    OP.normalize();
    SpMT Asp = OP.get_mat();
    MT A = MT(OP.get_mat());
    VT b = OP.get_vec();

    std::cout << "\n--- Testing Chebychev ball for sparse order Polytope" << std::endl;

    NT tol = 1e-08;
    unsigned int maxiter = 500;
    auto [center, radius, converged] =  max_inscribed_ball(Asp, b, maxiter, tol);
    auto [center2, radius2, converged2] =  max_inscribed_ball(A, b, maxiter, tol);

    VT center_(4);
    center_ << 0.207107, 0.5, 0.593398, 0.792893;

    CHECK(OP.is_in(Point(center)) == -1);
    auto [E, x0, round_val] = inscribed_ellipsoid_rounding<MT, VT, NT>(OP, Point(center));
    
    CHECK((center - center_).norm() <= 1e-06);
    CHECK(std::abs(radius - 0.207107) <= 1e-06);
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
    typedef Eigen::SparseMatrix<NT> SpMT;
    Hpolytope P;

    std::cout << "\n--- Testing analytic center for skinny H-polytope" << std::endl;
    bool pre_rounding = true; // round random polytope before applying the skinny transformation 
    NT max_min_eig_ratio = NT(100);
    P = skinny_random_hpoly<Hpolytope, NT, PolyRNGType>(3, 15, pre_rounding, max_min_eig_ratio, 127);
    P.normalize();
    
    auto [Hessian, analytic_center, converged] = analytic_center_linear_ineq<MT, MT, VT, NT>(P.get_mat(), P.get_vec());
    SpMT Asp = P.get_mat().sparseView();
    
    auto [Hessian_sp, analytic_center2, converged2] = analytic_center_linear_ineq<MT, SpMT, VT, NT>(Asp, P.get_vec());

    CHECK(P.is_in(Point(analytic_center)) == -1);
    CHECK(converged);
    CHECK(std::abs(analytic_center(0) + 4.75912) < 1e-04);
    CHECK(std::abs(analytic_center(1) + 4.28762) < 1e-04);
    CHECK(std::abs(analytic_center(2) - 7.54156) < 1e-04);

    CHECK(P.is_in(Point(analytic_center2)) == -1);
    CHECK(converged2);
    CHECK(std::abs(analytic_center(0) - analytic_center2(0)) < 1e-12);
    CHECK(std::abs(analytic_center(1) - analytic_center2(1)) < 1e-12);
    CHECK(std::abs(analytic_center(2) - analytic_center2(2)) < 1e-12);

    CHECK((Hessian - Hessian_sp).norm() < 1e-12);
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

TEST_CASE("test_max_ball_sparse") {
    call_test_max_ball_sparse<double>();
}
