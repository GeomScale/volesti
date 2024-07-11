// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

#include "doctest.h"
#include <fstream>
#include <iostream>

#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>

#include "misc/misc.h"
#include "random_walks/random_walks.hpp"

#include "volume/volume_sequence_of_balls.hpp"
#include "volume/volume_cooling_gaussians.hpp"
#include "volume/volume_cooling_balls.hpp"

#include "preprocess/min_sampling_covering_ellipsoid_rounding.hpp"
#include "preprocess/inscribed_ellipsoid_rounding.hpp"
#include "preprocess/svd_rounding.hpp"

#include "generators/known_polytope_generators.h"
#include "generators/h_polytopes_generator.h"

template <typename NT>
NT factorial(NT n)
{
    return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

template <typename NT>
void test_values(NT volume, NT expected, NT exact)
{
    std::cout << "Computed volume " << volume << std::endl;
    std::cout << "Expected volume = " << expected << std::endl;
    std::cout << "Relative error (expected) = "
              << std::abs((volume-expected)/expected) << std::endl;
    std::cout << "Relative error (exact) = "
              << std::abs((volume-exact)/exact) << std::endl;
    CHECK((std::abs((volume - exact)/exact) < 0.2 || 
            std::abs((volume - expected)/expected) < 0.00001));
}

template <class Polytope>
void rounding_min_ellipsoid_test(Polytope &HP,
                                 double const& expectedBall,
                                 double const& expectedCDHR,
                                 double const& expectedRDHR,
                                 double const& expectedBilliard,
                                 double const& exact)
{
    typedef typename Polytope::PointType Point;
    typedef typename Point::FT NT;
    typedef typename Polytope::MT MT;
    typedef typename Polytope::VT VT;

    int d = HP.dimension();

    typedef BoostRandomNumberGenerator<boost::mt19937, NT, 5> RNGType;
    RNGType rng(d);

    std::pair<Point, NT> InnerBall = HP.ComputeInnerBall();
    std::tuple<MT, VT, NT> res = min_sampling_covering_ellipsoid_rounding<CDHRWalk, MT, VT>(HP, InnerBall,
                                                                                            10 + 10 * d, rng);

    // Setup the parameters
    int walk_len = 1;
    NT e = 0.1;

    // Estimate the volume
    std::cout << "Number type: " << typeid(NT).name() << std::endl;


    //TODO: low accuracy in high dimensions
    //NT volume = res.second * volume_cooling_balls<BallWalk, RNGType>(HP, e, walk_len);
    //test_values(volume, expectedBall, exact);

    NT volume = std::get<2>(res) * volume_cooling_balls<CDHRWalk, RNGType>(HP, e, walk_len).second;
    test_values(volume, expectedCDHR, exact);

    volume = std::get<2>(res) * volume_cooling_balls<RDHRWalk, RNGType>(HP, e, 2*walk_len).second;
    test_values(volume, expectedRDHR, exact);

    volume = std::get<2>(res) * volume_cooling_balls<BilliardWalk, RNGType>(HP, e, walk_len).second;
    test_values(volume, expectedBilliard, exact);
}


template <class Polytope>
void rounding_max_ellipsoid_test(Polytope &HP,
                                 double const& expectedBall,
                                 double const& expectedCDHR,
                                 double const& expectedRDHR,
                                 double const& expectedBilliard,
                                 double const& exact)
{
    typedef typename Polytope::PointType Point;
    typedef typename Point::FT NT;
    typedef typename Polytope::MT MT;
    typedef typename Polytope::VT VT;

    int d = HP.dimension();

    typedef BoostRandomNumberGenerator<boost::mt19937, NT, 5> RNGType;
    RNGType rng(d);
    std::pair<Point, NT> InnerBall = HP.ComputeInnerBall();
    std::tuple<MT, VT, NT> res = inscribed_ellipsoid_rounding<MT, VT, NT>(HP, InnerBall.first);
    // Setup the parameters
    int walk_len = 1;
    NT e = 0.1;

    // Estimate the volume
    std::cout << "Number type: " << typeid(NT).name() << std::endl;

    NT volume = std::get<2>(res) * volume_cooling_balls<BilliardWalk, RNGType>(HP, e, walk_len).second;
    test_values(volume, expectedBilliard, exact);
}

template <class Polytope>
void rounding_max_ellipsoid_sparse_test(double const& expectedBilliard,
                                        double const& expected)
{
    typedef typename Polytope::PointType Point;
    typedef typename Point::FT NT;
    typedef typename Polytope::MT MT;
    typedef typename Polytope::VT VT;
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
    

    NT tol = 1e-08;
    unsigned int maxiter = 500;
    auto [center, radius, converged] =  max_inscribed_ball(Asp, OP.get_vec(), maxiter, tol);
    CHECK(OP.is_in(Point(center)) == -1);
    auto [E, x0, round_val] = inscribed_ellipsoid_rounding<MT, VT, NT>(OP, Point(center));

    MT A = MT(OP.get_mat());
    VT b = OP.get_vec();
    int d = OP.dimension();

    Polytope HP(d, A, b);

    typedef BoostRandomNumberGenerator<boost::mt19937, NT, 5> RNGType;
    // Setup the parameters
    int walk_len = 1;
    NT e = 0.1;

    // Estimate the volume
    std::cout << "Number type: " << typeid(NT).name() << std::endl;

    NT volume = round_val * volume_cooling_balls<BilliardWalk, RNGType>(HP, e, walk_len).second;
    test_values(volume, expectedBilliard, expected);
}

template <class Polytope>
void rounding_log_barrier_test(Polytope &HP,
                                 double const& expectedBall,
                                 double const& expectedCDHR,
                                 double const& expectedRDHR,
                                 double const& expectedBilliard,
                                 double const& exact)
{
    typedef typename Polytope::PointType Point;
    typedef typename Point::FT NT;
    typedef typename Polytope::MT MT;
    typedef typename Polytope::VT VT;

    int d = HP.dimension();

    typedef BoostRandomNumberGenerator<boost::mt19937, NT, 5> RNGType;
    RNGType rng(d);
    std::pair<Point, NT> InnerBall = HP.ComputeInnerBall();
    std::tuple<MT, VT, NT> res = inscribed_ellipsoid_rounding<MT, VT, NT, Polytope, Point, EllipsoidType::LOG_BARRIER>(HP, InnerBall.first);
    // Setup the parameters
    int walk_len = 1;
    NT e = 0.1;

    // Estimate the volume
    std::cout << "Number type: " << typeid(NT).name() << std::endl;

    NT volume = std::get<2>(res) * volume_cooling_balls<BilliardWalk, RNGType>(HP, e, walk_len).second;
    test_values(volume, expectedBilliard, exact);
}

template <class Polytope>
void rounding_volumetric_barrier_test(Polytope &HP,
                                 double const& expectedBall,
                                 double const& expectedCDHR,
                                 double const& expectedRDHR,
                                 double const& expectedBilliard,
                                 double const& exact)
{
    typedef typename Polytope::PointType Point;
    typedef typename Point::FT NT;
    typedef typename Polytope::MT MT;
    typedef typename Polytope::VT VT;

    int d = HP.dimension();

    typedef BoostRandomNumberGenerator<boost::mt19937, NT, 5> RNGType;
    RNGType rng(d);
    std::pair<Point, NT> InnerBall = HP.ComputeInnerBall();
    std::tuple<MT, VT, NT> res = inscribed_ellipsoid_rounding<MT, VT, NT, Polytope, Point, EllipsoidType::VOLUMETRIC_BARRIER>(HP, InnerBall.first);
    // Setup the parameters
    int walk_len = 1;
    NT e = 0.1;

    // Estimate the volume
    std::cout << "Number type: " << typeid(NT).name() << std::endl;

    NT volume = std::get<2>(res) * volume_cooling_balls<BilliardWalk, RNGType>(HP, e, walk_len).second;
    test_values(volume, expectedBilliard, exact);
}


template <class Polytope>
void rounding_svd_test(Polytope &HP,
                       double const& expectedBall,
                       double const& expectedCDHR,
                       double const& expectedRDHR,
                       double const& expectedBilliard,
                       double const& exact)
{
    typedef typename Polytope::PointType Point;
    typedef typename Point::FT NT;
    typedef typename Polytope::MT MT;
    typedef typename Polytope::VT VT;

    int d = HP.dimension();

    typedef BoostRandomNumberGenerator<boost::mt19937, NT, 5> RNGType;
    RNGType rng(d);

    std::pair<Point, NT> InnerBall = HP.ComputeInnerBall();
    std::tuple<MT, VT, NT> res = svd_rounding<CDHRWalk, MT, VT>(HP, InnerBall, 10 + 10 * d, rng);

    // Setup the parameters
    int walk_len = 1;
    NT e = 0.1;

    // Estimate the volume
    std::cout << "Number type: " << typeid(NT).name() << std::endl;

    NT volume = std::get<2>(res) * volume_cooling_balls<BilliardWalk, RNGType>(HP, e, walk_len).second;
    test_values(volume, expectedBilliard, exact);
}


template <typename NT>
void call_test_min_ellipsoid() {
    typedef Cartesian <NT> Kernel;
    typedef typename Kernel::Point Point;
    typedef HPolytope <Point> Hpolytope;
    Hpolytope P;

    std::cout << "\n--- Testing min ellipsoid rounding of H-skinny_cube5" << std::endl;
    P = generate_skinny_cube<Hpolytope>(5);
    rounding_min_ellipsoid_test(P, 0, 3070.64, 3188.25, 3140.6, 3200.0);

    std::cout << "\n--- Testing min ellipsoid rounding of H-skinny_cube10" << std::endl;

    P = generate_skinny_cube<Hpolytope>(10);
    rounding_min_ellipsoid_test(P, 0, 122550, 108426, 105003.0, 102400.0);
}


template <typename NT>
void call_test_max_ellipsoid() {
    typedef Cartesian <NT> Kernel;
    typedef typename Kernel::Point Point;
    typedef HPolytope <Point> Hpolytope;
    Hpolytope P;

    std::cout << "\n--- Testing max ellipsoid rounding of H-skinny_cube5" << std::endl;
    P = generate_skinny_cube<Hpolytope>(5);
    rounding_max_ellipsoid_test(P, 0, 3070.64, 3188.25, 3262.61, 3200.0);
}

template <typename NT>
void call_test_max_ellipsoid_sparse() {
    typedef Cartesian <NT> Kernel;
    typedef typename Kernel::Point Point;
    typedef HPolytope <Point> Hpolytope;

    std::cout << "\n--- Testing max ellipsoid rounding of sparse Order Polytope" << std::endl;
    rounding_max_ellipsoid_sparse_test<Hpolytope>(0.13979, 3070.64);
}

template <typename NT>
void call_test_log_barrier() {
    typedef Cartesian <NT> Kernel;
    typedef typename Kernel::Point Point;
    typedef HPolytope <Point> Hpolytope;
    Hpolytope P;

    std::cout << "\n--- Testing log-barrier rounding of H-skinny_cube5" << std::endl;
    P = generate_skinny_cube<Hpolytope>(5);
    rounding_log_barrier_test(P, 0, 3070.64, 3188.25, 3262.77, 3200.0);
}

template <typename NT>
void call_test_volumetric_barrier() {
    typedef Cartesian <NT> Kernel;
    typedef typename Kernel::Point Point;
    typedef HPolytope <Point> Hpolytope;
    Hpolytope P;

    std::cout << "\n--- Testing volumetric barrier rounding of H-skinny_cube5" << std::endl;
    P = generate_skinny_cube<Hpolytope>(5);
    rounding_volumetric_barrier_test(P, 0, 3070.64, 3188.25, 3262.77, 3200.0);
}


template <typename NT>
void call_test_svd() {
    typedef Cartesian <NT> Kernel;
    typedef typename Kernel::Point Point;
    typedef HPolytope <Point> Hpolytope;
    Hpolytope P;

    std::cout << "\n--- Testing SVD rounding of H-skinny_cube5" << std::endl;
    P = generate_skinny_cube<Hpolytope>(5);
    rounding_svd_test(P, 0, 3070.64, 3188.25, 3140.6, 3200.0);
}


TEST_CASE("round_min_ellipsoid") {
    call_test_min_ellipsoid<double>();
}

TEST_CASE("round_max_ellipsoid") {
    call_test_max_ellipsoid<double>();
}

TEST_CASE("round_max_ellipsoid_sparse") {
    call_test_max_ellipsoid_sparse<double>();
}

TEST_CASE("round_log_barrier_test") {
    call_test_log_barrier<double>();
}

TEST_CASE("round_volumetric_barrier_test") {
    call_test_volumetric_barrier<double>();
}

TEST_CASE("round_svd") {
    call_test_svd<double>();
}
