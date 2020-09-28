// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

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
    CHECK(((std::abs((volume - expected)/expected) < 0.00001) || (std::abs((volume - exact)/exact) < 0.2)));
}

template <class Polytope>
void test_volume(Polytope &P,
                 double const& expectedBall,
                 double const& expectedCDHR,
                 double const& expectedRDHR,
                 double const& exact)
{
    typedef typename Polytope::PointType Point;
    typedef typename Point::FT NT;

    // Setup the parameters
    int walk_len = 1;
    NT e=0.1;

    // Estimate the volume
    std::cout << "Number type: " << typeid(NT).name() << std::endl;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT, 3> RNGType;

    //TODO: low accuracy in high dimensions
    Polytope P1(P.dimension(), P.get_mat(), P.get_vec());
    NT volume = volume_cooling_gaussians<GaussianBallWalk, RNGType>(P1, e, walk_len);
    test_values(volume, expectedBall, exact);

    Polytope P2(P.dimension(), P.get_mat(), P.get_vec());
    volume = volume_cooling_gaussians<GaussianCDHRWalk, RNGType>(P2, e, walk_len);
    test_values(volume, expectedCDHR, exact);

    Polytope P3(P.dimension(), P.get_mat(), P.get_vec());
    volume = volume_cooling_gaussians<GaussianRDHRWalk, RNGType>(P3, e, walk_len);
    test_values(volume, expectedRDHR, exact);
}

template <typename NT>
void call_test_cube(){
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef boost::mt19937 RNGType;
    typedef VPolytope<Point> Vpolytope;
    Vpolytope P;

    std::cout << "--- Testing volume of V-cube2" << std::endl;
    P = generate_cube<Vpolytope>(2, true);
    test_volume(P, 3.88715, 4.04324, 4.33111, 4.0);

    std::cout << "--- Testing volume of V-cube4" << std::endl;
    P = generate_cube<Vpolytope>(4, true);
    test_volume(P,
                15.5448,
                15.4739,
                15.1616,
                16.0);
}

/*template <typename NT>
void call_test_cube_float(){
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef boost::mt19937 RNGType;
    typedef VPolytope<Point> Vpolytope;
    Vpolytope P;

    std::cout << "--- Testing volume of V-cube10 (float)" << std::endl;
    P = generate_cube<Vpolytope>(10, false);
    test_volume(P, 1000.55, 1024, 1024, 1024);

    std::cout << "--- Testing volume of V-cube20 (float)" << std::endl;
    P = generate_cube<Vpolytope>(20, false);
    test_volume(P, 1114192.7854272256,
                1048576,
                1048576,
                1048576);
}*/

template <typename NT>
void call_test_cross(){
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef boost::mt19937 RNGType;
    typedef VPolytope<Point> Vpolytope;
    Vpolytope P;

    std::cout << "--- Testing volume of V-cross5" << std::endl;
    P = generate_cross<Vpolytope>(5, true);
    test_volume(P,
                0.281565,
                0.23355,
                0.254434,
                0.266666667);

    /*std::cout << "--- Testing volume of V-cross10" << std::endl;
    P = generate_cross<Vpolytope>(10, true);
    test_volume(P,
                0.000309838,
                0.000311191,
                0.000299492,
                0.0002821869);*/
// both slow and inaccurate for CG
//    std::cout << "--- Testing volume of V-cross20" << std::endl;
//    P = gen_cross<Vpolytope>(20, true);
//    test_volume(P,
//                9.27132 * std::pow(10,-13),
//                4.42692 * std::pow(10,-13),
//                4.19453 * std::pow(10,-13),
//                std::pow(2.0,20.0) / factorial(20.0));
}

template <typename NT>
void call_test_simplex() {
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;

    typedef boost::mt19937 RNGType;
    typedef VPolytope<Point> Vpolytope;
    Vpolytope P;

    std::cout << "--- Testing volume of V-simplex5" << std::endl;
    P = generate_simplex<Vpolytope>(5, true);
    test_volume(P,
                0.00694196,
                0.00885402,
                0.00842199,
                1.0 / factorial(5.0));
// too slow for CG
/*
    std::cout << "--- Testing volume of V-simplex10" << std::endl;
    P = gen_simplex<Vpolytope>(10, true);
    test_volume(P,
                2.20484 * std::pow(10,-7),
                2.87992 * std::pow(10,-7),
                3.09976 * std::pow(10,-7),
                1.0 / factorial(10.0));

    std::cout << "--- Testing volume of V-simplex20" << std::endl;
    P = gen_simplex<Vpolytope>(20, true);
    test_volume(P,
                8.30917 * std::pow(10,-20),
                3.54876 * std::pow(10,-19),
                3.72986 * std::pow(10,-19),
                1.0 / factorial(20.0));
*/
}


TEST_CASE("cube") {
    call_test_cube<double>();
    //call_test_cube_float<float>();
}

TEST_CASE("cross") {
    call_test_cross<double>();
}

//TEST_CASE("simplex") {
//    call_test_simplex<double>();
//}
