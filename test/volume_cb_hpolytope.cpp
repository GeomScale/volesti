// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file


// Edited by HZ on 11.06.2020 - mute doctest.h
#include "doctest.h"
#include <fstream>
#include <iostream>
#include "misc.h"
#include "random.hpp"
#include "random/uniform_int.hpp"
#include "random/normal_distribution.hpp"
#include "random/uniform_real_distribution.hpp"

#include "random_walks/random_walks.hpp"

#include "volume/volume_sequence_of_balls.hpp"
#include "volume/volume_cooling_gaussians.hpp"
#include "volume/volume_cooling_balls.hpp"
#include "known_polytope_generators.h"

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
    CHECK((std::abs((volume - exact)/exact) < 0.3 || 
           std::abs((volume - expected)/expected) < 0.00001));
}

template <class Polytope>
void test_volume(Polytope &HP,
                 double const& expectedBall,
                 double const& expectedCDHR,
                 double const& expectedRDHR,
                 double const& expectedBilliard,
                 double const& exact,
                 bool birk = false)
{
    typedef typename Polytope::PointType Point;
    typedef typename Point::FT NT;

    // Setup the parameters
    int walk_len = 10 + HP.dimension()/10;
    NT e=0.1, volume;

    // Estimate the volume
    std::cout << "Number type: " << typeid(NT).name() << std::endl;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT, 3> RNGType;

    //TODO: low accuracy in high dimensions
    if (!birk) {
        volume = volume_cooling_balls<BallWalk, RNGType>(HP, e, walk_len);
        test_values(volume, expectedBall, exact);
    }

    volume = volume_cooling_balls<CDHRWalk, RNGType>(HP, e, walk_len);
    test_values(volume, expectedCDHR, exact);

    volume = volume_cooling_balls<RDHRWalk, RNGType>(HP, e, walk_len);
    test_values(volume, expectedRDHR, exact);

    volume = volume_cooling_balls<BilliardWalk, RNGType>(HP, e, walk_len);
    test_values(volume, expectedBilliard, exact);
}

template <typename NT>
void call_test_cube(){
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef HPolytope<Point> Hpolytope;
    Hpolytope P;

    std::cout << "--- Testing volume of H-cube10" << std::endl;
    P = gen_cube<Hpolytope>(10, false);
    test_volume(P, 1118.63, 1163.36, 1119.15, 1094.33, 1024);

    std::cout << "--- Testing volume of H-cube20" << std::endl;
    P = gen_cube<Hpolytope>(20, false);
    test_volume(P, 965744, 1051230, 1006470, 987076, 1048576);
}

template <typename NT>
void call_test_cube_float(){
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef HPolytope<Point> Hpolytope;
    Hpolytope P;
/*
    std::cout << "--- Testing volume of H-cube10" << std::endl;
    P = gen_cube<Hpolytope>(10, false);
    test_volume(P, 1000.55, 1024);

    std::cout << "--- Testing volume of H-cube20" << std::endl;
    P = gen_cube<Hpolytope>(20, false);
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
    Hpolytope P = gen_cross<Hpolytope>(10, false);
    test_volume(P,
                0.000291034,
                0.000281135,
                0.000294805,
                0.000286491,
                0.0002821869);
}

template <typename NT>
void call_test_birk()
{
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef HPolytope<Point> Hpolytope;
    Hpolytope P;

    typedef BoostRandomNumberGenerator<boost::mt19937, NT, 123> RNGType;

    std::cout << "--- Testing volume of H-birk3" << std::endl;
    P = gen_birk<Hpolytope>(3);
    test_volume(P, 0.114343, 0.125548, 0.113241, 0.112446, 0.125, true);

    std::cout << "--- Testing volume of H-birk4" << std::endl;
    P = gen_birk<Hpolytope>(4);
    test_volume(P, 0.00112956, 0.00109593, 0.00108152, 0.000845192,
                0.000970018, true);

    std::cout << "--- Testing volume of H-birk5" << std::endl;
    P = gen_birk<Hpolytope>(5);
    test_volume(P,
                1.97968e-07,
                2.55884e-07,
                1.39042e-07,
                2.11924e-07,
                0.000000225, 
                true);

    std::cout << "--- Testing volume of H-birk6" << std::endl;
    P = gen_birk<Hpolytope>(6);
    test_volume(P,
                7.84351e-13,
                6.10783e-13,
                4.39192e-13,
                6.62349e-13,
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
    P = gen_prod_simplex<Hpolytope>(5);
    test_volume(P,
                6.40072 * std::pow(10,-5),
                6.69062 * std::pow(10,-5),
                6.20744e-05,
                6.31986 * std::pow(10,-5),
                std::pow(1.0 / factorial(5.0), 2));

    std::cout << "--- Testing volume of H-prod_simplex10" << std::endl;
    P = gen_prod_simplex<Hpolytope>(10);
    test_volume(P,
                6.83631 * std::pow(10,-14),
                8.19581 * std::pow(10,-14),
                9.35005e-14,
                6.57309e-14,
                std::pow(1.0 / factorial(10.0), 2));

    std::cout << "--- Testing volume of H-prod_simplex15" << std::endl;
    P = gen_prod_simplex<Hpolytope>(15);
    test_volume(P,
                6.25978 * std::pow(10,-25),
                9.33162 * std::pow(10,-25),
                5.04617e-25,
                5.72542e-25,
                std::pow(1.0 / factorial(15.0), 2));
}

template <typename NT>
void call_test_simplex() {
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;

    typedef HPolytope<Point> Hpolytope;
    Hpolytope P;

    std::cout << "--- Testing volume of H-simplex10" << std::endl;
    P = gen_simplex<Hpolytope>(10, false);
    test_volume(P,
                3.22432 * std::pow(10,-7),
                2.90617 * std::pow(10,-7),
                2.93392 * std::pow(10,-7),
                3.00286e-07,
                1.0 / factorial(10.0));

    std::cout << "--- Testing volume of H-simplex20" << std::endl;
    P = gen_simplex<Hpolytope>(20, false);
    test_volume(P,
                4.03788 * std::pow(10,-19),
                4.14182 * std::pow(10,-19),
                4.5877e-19,
                4.54245e-19,
                1.0 / factorial(20.0));

    std::cout << "--- Testing volume of H-simplex30" << std::endl;
    P = gen_simplex<Hpolytope>(30, false);
    test_volume(P,
                2.5776 * std::pow(10,-33),
                3.5157 * std::pow(10,-33),
                2.74483e-33,
                3.08769e-33,
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

