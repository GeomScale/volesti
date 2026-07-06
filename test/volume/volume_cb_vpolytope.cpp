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
    CHECK(std::abs((volume - exact)/exact) < 0.2);
}

template <class Polytope>
void test_volume(Polytope &P,
                 double const& expectedBall,
                 double const& expectedCDHR,
                 double const& expectedRDHR,
                 double const& expectedBilliard,
                 double const& exact)
{
    typedef typename Polytope::PointType Point;
    typedef typename Point::FT NT;

    // Setup the parameters
    int walk_len = 10 + P.dimension()/10;
    NT e=0.1;

    // Estimate the volume
    std::cout << "Number type: " << typeid(NT).name() << std::endl;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT, 3> RNGType;

    //TODO: low accuracy in high dimensions
    Polytope P1(P.dimension(), P.get_mat(), P.get_vec());
    NT volume = volume_cooling_balls<BallWalk, RNGType>(P1, e, walk_len).second;
    test_values(volume, expectedBall, exact);

    Polytope P2(P.dimension(), P.get_mat(), P.get_vec());
    volume = volume_cooling_balls<CDHRWalk, RNGType>(P2, e, walk_len).second;
    test_values(volume, expectedCDHR, exact);

    Polytope P3(P.dimension(), P.get_mat(), P.get_vec());
    volume = volume_cooling_balls<RDHRWalk, RNGType>(P3, e, walk_len).second;
    test_values(volume, expectedRDHR, exact);

    Polytope P4(P.dimension(), P.get_mat(), P.get_vec());
    volume = volume_cooling_balls<BilliardWalk, RNGType>(P4, e, walk_len).second;
    test_values(volume, expectedBilliard, exact);
}

template <typename NT>
void call_test_cube(){
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef boost::mt19937 RNGType;
    typedef VPolytope<Point> Vpolytope;


    std::cout << "--- Testing volume of V-cube2" << std::endl;
    Vpolytope P1 = generate_cube<Vpolytope>(2, true);
    test_volume(P1, 4.43443, 4.129, 4.43443, 4.40191, 4);

    std::cout << "--- Testing volume of V-cube5" << std::endl;
    Vpolytope P2 = generate_cube<Vpolytope>(5, true);
    test_volume(P2, 32, 32, 32, 32, 32);


}

template <typename NT>
void call_test_cube_float(){
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef boost::mt19937 RNGType;
    typedef VPolytope<Point> Vpolytope;
    Vpolytope P;

    std::cout << "--- Testing volume of V-cube10 (float)" << std::endl;
    P = generate_cube<Vpolytope>(10, true);
    test_volume(P, 1000.55, 1024, 1024, 1024, 1024);

    std::cout << "--- Testing volume of V-cube20 (float)" << std::endl;
    P = generate_cube<Vpolytope>(20, true);
    test_volume(P, 1114192.7854272256,
                1048576,
                1048576,
                1048576,
                1048576);
}

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
                0.28425,
                0.273255,
                0.28413,
                0.286071,
                0.266666667);

    std::cout << "--- Testing volume of V-cross10" << std::endl;
    P = generate_cross<Vpolytope>(10, true);
    test_volume(P,
                0.000283841,
                0.00031188,
                0.000284841,
                0.00027759,
                0.0002821869);

    std::cout << "--- Testing volume of V-cross20" << std::endl;
    P = generate_cross<Vpolytope>(20, true);
    test_volume(P,
                4.16807 * std::pow(10,-13),
                4.42692 * std::pow(10,-13),
                4.19453 * std::pow(10,-13),
                4.63423 * std::pow(10,-13),
                std::pow(2.0,20.0) / factorial(20.0));
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
                0.00846587,
                0.0096107,
                0.00842591,
                0.00855401,
                1.0 / factorial(5.0));

    std::cout << "--- Testing volume of V-simplex10" << std::endl;
    P = generate_simplex<Vpolytope>(10, true);
    test_volume(P,
                2.35669 * std::pow(10,-7),
                3.00778 * std::pow(10,-7),
                3.0366 * std::pow(10,-7),
                2.72952 * std::pow(10,-7),
                1.0 / factorial(10.0));
/* too slow
    std::cout << "--- Testing volume of V-simplex20" << std::endl;
    P = gen_simplex<Vpolytope>(20, true);
    test_volume(P,
                1.13981 * std::pow(10,-19),
                3.63355 * std::pow(10,-19),
                4.46469 * std::pow(10,-19),
                4.22932 * std::pow(10,-19),
                1.0 / factorial(20.0));
*/
}


TEST_CASE("cube") {
    //TODO: Runtime error, check ComputeInnerBall()
    call_test_cube<double>();
    //call_test_cube_float<float>();
}

TEST_CASE("cross") {
    call_test_cross<double>();
}

TEST_CASE("simplex") {
    call_test_simplex<double>();
}
