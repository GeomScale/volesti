#include "Eigen/Eigen"
#include <iostream>
#include "cartesian_geom/cartesian_kernel.h"
#include "convex_bodies/hpolytope.h"
#include "generators/known_polytope_generators.h"
#include "random_walks/random_walks.hpp"
#include "misc/misc.h"
#include "doctest.h"  

typedef double NT;
typedef Cartesian<NT> Kernel;
typedef typename Kernel::Point Point;
typedef BoostRandomNumberGenerator<boost::mt19937, NT, 3> RNGType;
typedef HPolytope<Point> Polytope;

void test_gaussian_hmc_reflections(const std::string& polytope_type, RNGType& rng, unsigned int walk_length) {
    Polytope P;
    if (polytope_type == "cube") {
        P = generate_cube<Polytope>(3, false, 10.0);
    } else if (polytope_type == "cross") {
        P = generate_cross<Polytope>(3, false);
    } else if (polytope_type == "simplex") {
        P = generate_simplex<Polytope>(3, false);
    } else {
        std::cout << "Unknown polytope type: " << polytope_type << std::endl;
        return;
    }

    Point interior_point = P.ComputeInnerBall().first;
    double L = 1.0; // Step size
    unsigned int rho = 100; // Max reflections
    typename GaussianHamiltonianMonteCarloExactWalk::parameters params(L, true, rho, true);
    typename GaussianHamiltonianMonteCarloExactWalk::Walk<Polytope, RNGType> walk(P, interior_point, 1.0, rng, params);
    
    walk.apply(P, interior_point, 1.0, walk_length, rng);
    
    CHECK(P.is_in(interior_point) == -1);  // Use CHECK from doctest
    std::cout << "Test passed for " << polytope_type << ": Point remains inside the polytope after reflections.\n";
}

TEST_CASE("GHMC Reflections Test - Cube") {
    RNGType rng(123456);
    test_gaussian_hmc_reflections("cube", rng, 10);
}

TEST_CASE("GHMC Reflections Test - Cross Polytope") {
    RNGType rng(123456);
    test_gaussian_hmc_reflections("cross", rng, 10);
}

TEST_CASE("GHMC Reflections Test - Simplex") {
    RNGType rng(123456);
    test_gaussian_hmc_reflections("simplex", rng, 10);
}
