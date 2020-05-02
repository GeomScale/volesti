// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

#include "doctest.h"
#include <fstream>
#include <iostream>
#include "misc.h"
#include "new_volume.hpp"
#include "new_gaussian_volume.hpp"
#include "new_cooling_balls.hpp"
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
    CHECK(std::abs((volume - expected)/expected) < 0.00001);
}

template <class Polytope>
void test_volume(Polytope &HP,
                 double const& expectedBall,
                 double const& expectedCDHR,
                 double const& expectedRDHR,
                 double const& expectedBilliard,
                 double const& exact)
{
    typedef typename Polytope::PointType Point;
    typedef typename Point::FT NT;

    // Setup the parameters
    int walk_len = 10 + HP.dimension()/10;
    NT e=0.1;

    // Estimate the volume
    std::cout << "Number type: " << typeid(NT).name() << std::endl;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT, 3> RNGType;

    //TODO: low accuracy in high dimensions
    NT volume = volume_cooling_balls<BallWalk, RNGType>(HP, e, walk_len);
    test_values(volume, expectedBall, exact);

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
    test_volume(P, 1107.81, 1199.29, 1097.47, 1097.26, 1024);

    std::cout << "--- Testing volume of H-cube20" << std::endl;
    P = gen_cube<Hpolytope>(20, false);
    test_volume(P, 1051340, 1017430, 1083870, 1094580, 1048576);
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
                0.000290427,
                0.000283703,
                0.000293307,
                0.000280408,
                0.0002821869);
}

template <typename NT>
void call_test_birk()
{
    //TODO: All but Ball walk compute wrong results

    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef HPolytope<Point> Hpolytope;
    Hpolytope P;

    typedef BoostRandomNumberGenerator<boost::mt19937, NT, 123> RNGType;

    std::cout << "--- Testing volume of H-birk3" << std::endl;
    std::ifstream inp;
    std::vector<std::vector<NT> > Pin;
    inp.open("../R-proj/inst/extdata/birk3.ine",std::ifstream::in);
    read_pointset(inp,Pin);
    P.init(Pin);
    test_volume(P, 0.118464, 0.115994, 0.119691, 0.12394, 0.125);

    std::cout << "--- Testing volume of H-birk4" << std::endl;
    std::ifstream inp2;
    std::vector<std::vector<NT> > Pin2;
    inp2.open("../R-proj/inst/extdata/birk4.ine",std::ifstream::in);
    read_pointset(inp2,Pin2);
    P.init(Pin2);
    test_volume(P, 0.000856744, 0.000840153, 0.00122184, 0.000957595,
                0.000970018);

    std::cout << "--- Testing volume of H-birk5" << std::endl;
    std::ifstream inp3;
    std::vector<std::vector<NT> > Pin3;
    inp3.open("../R-proj/inst/extdata/birk5.ine",std::ifstream::in);
    read_pointset(inp3,Pin3);
    P.init(Pin3);
    test_volume(P,
                6.12191 * std::pow(10,-8),
                1.97637 * std::pow(10,-7),
                2.5704 * std::pow(10,-7),
                2.08666 * std::pow(10,-7),
                0.000000225);

    std::cout << "--- Testing volume of H-birk6" << std::endl;
    std::ifstream inp4;
    std::vector<std::vector<NT> > Pin4;
    inp4.open("../R-proj/inst/extdata/birk6.ine",std::ifstream::in);
    read_pointset(inp4,Pin4);
    P.init(Pin4);
    test_volume(P,
                1.99583 * std::pow(10,-18),
                7.24084 * std::pow(10,-13),
                5.73665 * std::pow(10,-13),
                1.0156 * std::pow(10,-12),
                0.0000000000009455459196);
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
                7.19668 * std::pow(10,-5),
                6.42808 * std::pow(10,-5),
                7.3549 * std::pow(10,-5),
                6.65353 * std::pow(10,-5),
                std::pow(1.0 / factorial(5.0), 2));

    std::cout << "--- Testing volume of H-prod_simplex10" << std::endl;
    P = gen_prod_simplex<Hpolytope>(10);
    test_volume(P,
                3.48279 * std::pow(10,-14),
                8.06992 * std::pow(10,-14),
                9.55515 * std::pow(10,-14),
                7.11823 * std::pow(10,-14),
                std::pow(1.0 / factorial(10.0), 2));

    std::cout << "--- Testing volume of H-prod_simplex15" << std::endl;
    P = gen_prod_simplex<Hpolytope>(15);
    test_volume(P,
                1.12426 * std::pow(10,-25),
                7.62343 * std::pow(10,-25),
                4.59091 * std::pow(10,-25),
                5.3557 * std::pow(10,-25),
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
                2.27161 * std::pow(10,-7),
                3.18695 * std::pow(10,-7),
                2.82644 * std::pow(10,-7),
                2.68846 * std::pow(10,-7),
                1.0 / factorial(10.0));

    std::cout << "--- Testing volume of H-simplex20" << std::endl;
    P = gen_simplex<Hpolytope>(20, false);
    test_volume(P,
                1.62349 * std::pow(10,-19),
                4.63344 * std::pow(10,-19),
                4.09729 * std::pow(10,-19),
                4.09423 * std::pow(10,-19),
                1.0 / factorial(20.0));

    std::cout << "--- Testing volume of H-simplex30" << std::endl;
    P = gen_simplex<Hpolytope>(30, false);
    test_volume(P,
                2.0167 * std::pow(10,-33),
                3.88403 * std::pow(10,-33),
                3.79774 * std::pow(10,-33),
                3.57908 * std::pow(10,-33),
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

