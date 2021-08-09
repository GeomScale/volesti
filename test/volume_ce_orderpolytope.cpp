// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis
// Copyright (c) 2021 Vaibhav Thakkar

// Licensed under GNU LGPL.3, see LICENCE file

#include "doctest.h"
#include <fstream>
#include <iostream>
#include "misc.h"
#include "random.hpp"
#include "random/uniform_int.hpp"
#include "random/normal_distribution.hpp"
#include "random/uniform_real_distribution.hpp"

#include "convex_bodies/orderpolytope.h"
#include "misc/poset.h"
#include "random_walks/gaussian_accelerated_billiard_walk.hpp"
#include "volume/volume_cooling_ellipsoids.hpp"
#include "generators/order_polytope/order_polytope_generator.h"

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
    CHECK((std::abs((volume - exact)/exact) < 0.35 ||
           std::abs((volume - expected)/expected) < 0.00001));
}

template <class Polytope>
void test_volume(Polytope &OP,
                 double const& expectedGaussianAcceleratedBilliard,
                 double const& exact)
{
    typedef typename Polytope::PointType Point;
    typedef typename Point::FT NT;

    // Setup the parameters
    int walk_len = 10 + OP.dimension()/10;
    NT e=0.1, volume;

    // Estimate the volume
    std::cout << "Number type: " << typeid(NT).name() << std::endl;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT, 5> RNGType;

    volume = volume_cooling_ellipsoids<GaussianAcceleratedBilliardWalk, RNGType>(OP, e, walk_len).second;
    test_values(volume, expectedGaussianAcceleratedBilliard, exact);
}

template <typename NT>
void call_test_cube(){
    typedef Cartesian<NT> Kernel;
    typedef typename Kernel::Point Point;
    typedef OrderPolytope<Point> Orderpolytope;

    std::cout << "--- Testing volume of OrderPolytope-cube10" << std::endl;
    Orderpolytope OP = generate_cube_orderpoly<Orderpolytope>(10);
    OP.normalize();
    test_volume(OP, 1.07956, 1.0);

    std::cout << "--- Testing volume of OrderPolytope-cube20" << std::endl;
    OP = generate_cube_orderpoly<Orderpolytope>(20);
    OP.normalize();
    test_volume(OP, 0.958769, 1.0);
}

template <typename NT>
void call_test_bipartite(){
    typedef Cartesian<NT> Kernel;
    typedef typename Kernel::Point Point;
    typedef OrderPolytope<Point> Orderpolytope;

    std::cout << "--- Testing volume of OrderPolytope-biparitite_0.5_008_0" << std::endl;
    std::string filepath = "../../include/generators/order_polytope/instances/bipartite_0.5_008_0.txt";
    Orderpolytope OP = generate_orderpoly<Orderpolytope>(filepath);
    OP.normalize();
    test_volume(OP, 0.0382686, 1504.0/factorial(8));    // 1504 is the number of linear extensions of this poset
}

template <typename NT>
void call_test_andes(){
    typedef Cartesian<NT> Kernel;
    typedef typename Kernel::Point Point;
    typedef OrderPolytope<Point> Orderpolytope;

    std::cout << "--- Testing volume of OrderPolytope-bayesiannetwork_andes_008_0" << std::endl;
    std::string filepath = "../../include/generators/order_polytope/instances/bayesiannetwork_andes_008_0.txt";
    Orderpolytope OP = generate_orderpoly<Orderpolytope>(filepath);
    OP.normalize();
    test_volume(OP, 0.000741068, 28.0/factorial(8));    // 28 is the number of linear extensions of this poset
}


TEST_CASE("cube") {
    call_test_cube<double>();
}

TEST_CASE("bipartite") {
    call_test_bipartite<double>();
}

TEST_CASE("andes") {
    call_test_andes<double>();
}

