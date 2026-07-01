// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

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
#include "generators/known_polytope_generators.h"

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
void test_volume(Polytope &HP,
                 double const& expectedBall,
                 double const& expectedCDHR,
                 double const& expectedRDHR,
                 double const& exact,
                 bool birk = false)
{
    typedef typename Polytope::PointType Point;
    typedef typename Point::FT NT;

    // Setup the parameters
    int walk_len = 10 + HP.dimension()/10;
    NT e=0.1;
    NT volume;

    // Estimate the volume
    std::cout << "Number type: " << typeid(NT).name() << std::endl;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT, 3> RNGType;

    // TODO: low accuracy in high-dimensions
    if (!birk) {
        volume = volume_cooling_gaussians<GaussianBallWalk, RNGType>(HP, e, walk_len);
        test_values(volume, expectedBall, exact);
    }

    volume = volume_cooling_gaussians<GaussianCDHRWalk, RNGType>(HP, e, walk_len);
    test_values(volume, expectedCDHR, exact);

    volume = volume_cooling_gaussians<GaussianRDHRWalk, RNGType>(HP, e, walk_len);
    test_values(volume, expectedRDHR, exact);
}

template <typename NT>
void call_test_cube(){
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef HPolytope<Point> Hpolytope;
    Hpolytope P;

    std::cout << "--- Testing volume of H-cube10" << std::endl;
    P = generate_cube<Hpolytope>(10, false);
    test_volume(P, 1079.56, 1110.92, 1113.93, 1024);

    std::cout << "--- Testing volume of H-cube20" << std::endl;
    P = generate_cube<Hpolytope>(20, false);
    test_volume(P, 1.1025e+06, 1.05174e+06, 995224, 1048576);
}

template <typename NT>
void call_test_cube_float(){
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef HPolytope<Point> Hpolytope;
    Hpolytope P;
/*
    std::cout << "--- Testing volume of H-cube10" << std::endl;
    P = generate_cube<Hpolytope>(10, false);
    test_volume(P, 1000.55, 1024);

    std::cout << "--- Testing volume of H-cube20" << std::endl;
    P = generate_cube<Hpolytope>(20, false);
    test_volume(P, 1114192.7854272256, 1048576);
    */
}

template <typename NT>
void call_test_cross(){
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef HPolytope<Point> Hpolytope;

    typedef BoostRandomNumberGenerator<boost::mt19937, NT, 123> RNGType;

    std::cout << "--- Testing volume of H-cross10" << std::endl;
    Hpolytope P = generate_cross<Hpolytope>(10, false);
    test_volume(P,
                0.000292199,
                0.000274014,
                0.000294463,
                0.0002821869);
}

template <typename NT>
void call_test_birk() {
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef HPolytope<Point> Hpolytope;
    Hpolytope P;

    typedef BoostRandomNumberGenerator<boost::mt19937, NT, 123> RNGType;

    std::cout << "--- Testing volume of H-birk3" << std::endl;
    P = generate_birkhoff<Hpolytope>(3);
    test_volume(P, 0.116678, 0.122104, 0.11326, 0.125, true);

    std::cout << "--- Testing volume of H-birk4" << std::endl;
    P = generate_birkhoff<Hpolytope>(4);
    test_volume(P,
                0.000450761,
                0.00108943,
                0.00110742,
                0.000970018,
                true);

    std::cout << "--- Testing volume of H-birk5" << std::endl;
    P = generate_birkhoff<Hpolytope>(5);
    test_volume(P,
                2.97522e-08,
                2.00743e-07,
                2.05779e-07,
                2.25  * std::pow(10,-7),
                true);

    std::cout << "--- Testing volume of H-birk6" << std::endl;
    P = generate_birkhoff<Hpolytope>(6);
    test_volume(P,
                3.66375e-19,
                7.51051 * std::pow(10,-13),
                8.20587e-13,
                9.455459196 * std::pow(10,-13),
                true);
}

template <typename NT>
void call_test_prod_simplex() {
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;

    typedef HPolytope<Point> Hpolytope;
    Hpolytope P;

    std::cout << "--- Testing volume of H-prod_simplex5" << std::endl;
    P = generate_prod_simplex<Hpolytope>(5);
    test_volume(P,
                6.3448 * std::pow(10,-5),
                6.94695 * std::pow(10,-5),
                6.13242e-05,
                std::pow(1.0 / factorial(5.0), 2));

    std::cout << "--- Testing volume of H-prod_simplex10" << std::endl;
    P = generate_prod_simplex<Hpolytope>(10);
    test_volume(P,
                1.66017 * std::pow(10,-14),
                8.48116 * std::pow(10,-14),
                6.90898e-14,
                std::pow(1.0 / factorial(10.0), 2));

    std::cout << "--- Testing volume of H-prod_simplex15" << std::endl;
    P = generate_prod_simplex<Hpolytope>(15);
    test_volume(P,
                2.0232 * std::pow(10,-29),
                5.4624 * std::pow(10,-25),
                6.95082e-25,
                std::pow(1.0 / factorial(15.0), 2));
}

template <typename NT>
void call_test_simplex() {
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;

    typedef HPolytope<Point> Hpolytope;
    Hpolytope P;

    std::cout << "--- Testing volume of H-simplex10" << std::endl;
    P = generate_simplex<Hpolytope>(10, false);
    test_volume(P,
                2.14048 * std::pow(10,-7),
                2.70598 * std::pow(10,-7),
                2.53893e-07,
                1.0 / factorial(10.0));

    std::cout << "--- Testing volume of H-simplex20" << std::endl;
    P = generate_simplex<Hpolytope>(20, false);
    test_volume(P,
                2.00646 * std::pow(10,-21),
                4.16845 * std::pow(10,-19),
                3.79918e-19,
                1.0 / factorial(20.0));

    std::cout << "--- Testing volume of H-simplex30" << std::endl;
    P = generate_simplex<Hpolytope>(30, false);
    test_volume(P,
                2.31348 * std::pow(10,-35),
                4.02288 * std::pow(10,-33),
                3.47743e-33,
                1.0 / factorial(30.0));
}

template <typename NT>
void call_test_skinny_cube() {
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT, 123> RNGType;

    typedef HPolytope<Point> Hpolytope;
    Hpolytope P;

    std::cout << "--- Testing volume of H-skinny_cube10" << std::endl;
    //TODO: needs rounding
    //P = gen_skinny_cube<Hpolytope>(10);
    //test_volume(P, 15591.1, 102400.0);

    //std::cout << "--- Testing volume of H-skinny_cube20" << std::endl;
    //P = gen_skinny_cube<Hpolytope>(20);
    //test_volume(P, 104857600, 104857600.0);
}

template <typename NT>
void call_test_sparse_simplex() {
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;

    typedef HPolytope<Point, Eigen::SparseMatrix<NT>> Hpolytope;
    Hpolytope SP;

    std::cout << "--- Testing volume of sparse H-simplex10" << std::endl;
    SP = generate_simplex<Hpolytope>(10, false);
    test_volume(SP,
                2.14048 * std::pow(10,-7),
                2.70598 * std::pow(10,-7),
                2.53893e-07,
                1.0 / factorial(10.0));

    std::cout << "--- Testing volume of sparse H-simplex20" << std::endl;
    SP = generate_simplex<Hpolytope>(20, false);
    test_volume(SP,
                2.00646 * std::pow(10,-21),
                4.16845 * std::pow(10,-19),
                5.0348e-19,
                1.0 / factorial(20.0));
}



TEST_CASE("cube") {
    call_test_cube<double>();
    call_test_cube_float<float>();
}

TEST_CASE("cross") {
    call_test_cross<double>();
}

TEST_CASE("birk") {
    call_test_birk<double>();
}

TEST_CASE("prod_simplex") {
    call_test_prod_simplex<double>();
}

TEST_CASE("simplex") {
    call_test_simplex<double>();
}

TEST_CASE("skinny_cube") {
    call_test_skinny_cube<double>();
}

TEST_CASE("sparse_simplex") {
    call_test_sparse_simplex<double>();
}
