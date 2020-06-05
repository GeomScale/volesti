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

template <class Polytope>
void test_volume(Polytope &P, double const& expected, double const& exact)
{
    typedef typename Polytope::PointType Point;
    typedef typename Point::FT NT;

    // Setup the parameters
    int walk_len = 10 + P.dimension()/10;
    NT e=1;

    // Estimate the volume
    std::cout << "Number type: " << typeid(NT).name() << std::endl;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT, 3> RNGType;

    Polytope P1(P.dimension(), P.get_mat(), P.get_vec());
    NT volume = volume_sequence_of_balls<RDHRWalk, RNGType>(P1, e, walk_len);

    //TODO: test other walks

    std::cout << "Computed volume " << volume << std::endl;
    std::cout << "Expected volume = " << expected << std::endl;
    std::cout << "Relative error (expected) = "
              << std::abs((volume-expected)/expected) << std::endl;
    std::cout << "Relative error (exact) = "
              << std::abs((volume-exact)/exact) << std::endl;
    CHECK(std::abs((volume - expected)/expected) < 0.00001);
}


template <typename NT>
void call_test_cube(){
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef boost::mt19937 RNGType;
    typedef VPolytope<Point> Vpolytope;
    Vpolytope P;

    std::cout << "--- Testing volume of V-cube10" << std::endl;
    P = generate_cube<Vpolytope>(10, true);
    test_volume(P, 1096.5089688155, 1024);

    std::cout << "--- Testing volume of V-cube20" << std::endl;
    P = generate_cube<Vpolytope>(20, true);
    test_volume(P, 967352.7854272256, 1048576);
}

template <typename NT>
void call_test_cube_float(){
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef boost::mt19937 RNGType;
    typedef VPolytope<Point> Vpolytope;

    std::cout << "--- Testing volume of V-cube10 (float)" << std::endl;
    Vpolytope P1 = generate_cube<Vpolytope>(10, false);
    test_volume(P1, 1000.55, 1024);

    std::cout << "--- Testing volume of V-cube20 (float)" << std::endl;
    Vpolytope P2 = generate_cube<Vpolytope>(20, false);
    test_volume(P2, 1114192.7854272256, 1048576);
}

template <typename NT>
void call_test_cross(){
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef boost::mt19937 RNGType;
    typedef VPolytope<Point> Vpolytope;

    std::cout << "--- Testing volume of V-cross5" << std::endl;
    Vpolytope P1 = gen_cross<Vpolytope>(5, true);
    test_volume(P1, 0.276845, 0.266666667);

    std::cout << "--- Testing volume of V-cross10" << std::endl;
    Vpolytope P2 = gen_cross<Vpolytope>(10, true);
    test_volume(P2, 0.000291003, 0.0002821869);
}

template <typename NT>
void call_test_simplex() {
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;

    typedef boost::mt19937 RNGType;
    typedef VPolytope<Point> Vpolytope;
    Vpolytope P;

    std::cout << "--- Testing volume of V-simplex5" << std::endl;
    P = gen_simplex<Vpolytope>(5, true);
    test_volume(P, 0.00810133, 1.0 / factorial(5.0));


    // too slow for SoB
//    std::cout << "--- Testing volume of V-simplex10" << std::endl;
//    P = gen_simplex<Vpolytope>(10, true);
//    test_volume(P, 2.99056 * std::pow(10,-7), 1.0 / factorial(10.0));

}


TEST_CASE("cube") {
    //TODO: Runtime error, check ComputeInnerBall()
//    call_test_cube<double>();
    //call_test_cube_float<float>();
}

TEST_CASE("cross") {
    call_test_cross<double>();
}

TEST_CASE("simplex") {
    call_test_simplex<double>();
}
