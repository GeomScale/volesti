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

    // TODO: low accuracy in high-dimensions
    NT volume = volume_gaussian_annealing<GaussianBallWalk, RNGType>(HP, e, walk_len);
    test_values(volume, expectedBall, exact);

    volume = volume_gaussian_annealing<GaussianCDHRWalk, RNGType>(HP, e, walk_len);
    test_values(volume, expectedCDHR, exact);

    volume = volume_gaussian_annealing<GaussianRDHRWalk, RNGType>(HP, e, walk_len);
    test_values(volume, expectedRDHR, exact);
}

template <typename NT>
void call_test_cube(){
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef HPolytope<Point> Hpolytope;
    Hpolytope P;

    std::cout << "--- Testing volume of H-cube10" << std::endl;
    P = gen_cube<Hpolytope>(10, false);
    test_volume(P, 1102.47, 1119.38, 1037.34, 1024);

    std::cout << "--- Testing volume of H-cube20" << std::endl;
    P = gen_cube<Hpolytope>(20, false);
    test_volume(P, 1104980, 1051120, 989794, 1048576);
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
                0.000280513,
                0.000274237,
                0.000296688,
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
    std::ifstream inp;
    std::vector<std::vector<NT> > Pin;
    inp.open("../R-proj/inst/extdata/birk3.ine",std::ifstream::in);
    read_pointset(inp,Pin);
    P.init(Pin);
    test_volume(P, 0.101546, 0.121466, 0.117802, 0.125);

    std::cout << "--- Testing volume of H-birk4" << std::endl;
    std::ifstream inp2;
    std::vector<std::vector<NT> > Pin2;
    inp2.open("../R-proj/inst/extdata/birk4.ine",std::ifstream::in);
    read_pointset(inp2,Pin2);
    P.init(Pin2);
    test_volume(P,
                0.000942906,
                0.00114121,
                0.00104026,
                0.000970018);

    std::cout << "--- Testing volume of H-birk5" << std::endl;
    std::ifstream inp3;
    std::vector<std::vector<NT> > Pin3;
    inp3.open("../R-proj/inst/extdata/birk5.ine",std::ifstream::in);
    read_pointset(inp3,Pin3);
    P.init(Pin3);
    test_volume(P,
                1.08184 * std::pow(10,-7),
                2.31411 * std::pow(10,-7),
                2.39421 * std::pow(10,-7),
                2.25  * std::pow(10,-7));

    std::cout << "--- Testing volume of H-birk6" << std::endl;
    std::ifstream inp4;
    std::vector<std::vector<NT> > Pin4;
    inp4.open("../R-proj/inst/extdata/birk6.ine",std::ifstream::in);
    read_pointset(inp4,Pin4);
    P.init(Pin4);
    test_volume(P,
                5.35067 * std::pow(10,-17),
                8.44699 * std::pow(10,-13),
                9.58445 * std::pow(10,-13),
                9.455459196 * std::pow(10,-13));
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
                8.61492 * std::pow(10,-5),
                7.36265 * std::pow(10,-5),
                7.07551 * std::pow(10,-5),
                std::pow(1.0 / factorial(5.0), 2));

    std::cout << "--- Testing volume of H-prod_simplex10" << std::endl;
    P = gen_prod_simplex<Hpolytope>(10);
    test_volume(P,
                8.9345 * std::pow(10,-14),
                8.49473 * std::pow(10,-14),
                7.36811 * std::pow(10,-14),
                std::pow(1.0 / factorial(10.0), 2));

    std::cout << "--- Testing volume of H-prod_simplex15" << std::endl;
    P = gen_prod_simplex<Hpolytope>(15);
    test_volume(P,
                6.37869 * std::pow(10,-25),
                5.55858 * std::pow(10,-25),
                5.93187 * std::pow(10,-25),
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
                3.0842 * std::pow(10,-7),
                2.59229 * std::pow(10,-7),
                2.68156 * std::pow(10,-7),
                1.0 / factorial(10.0));

    std::cout << "--- Testing volume of H-simplex20" << std::endl;
    P = gen_simplex<Hpolytope>(20, false);
    test_volume(P,
                3.41016 * std::pow(10,-19),
                4.16827 * std::pow(10,-19),
                4.43788 * std::pow(10,-19),
                1.0 / factorial(20.0));

    std::cout << "--- Testing volume of H-simplex30" << std::endl;
    P = gen_simplex<Hpolytope>(30, false);
    test_volume(P,
                3.30097 * std::pow(10,-33),
                4.15532 * std::pow(10,-33),
                3.54527 * std::pow(10,-33),
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
