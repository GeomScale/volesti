#include "Eigen/Eigen"
#include <iostream>
#include "random.hpp"
#include "random/uniform_int.hpp"
#include "random/normal_distribution.hpp"
#include "cartesian_geom/cartesian_kernel.h"
#include "convex_bodies/hpolytope.h"
#include "generators/known_polytope_generators.h"
#include "generators/h_polytopes_generator.h" // Include your polytope generator
#include "random_walks/random_walks.hpp"
#include "misc/misc.h"
#include "doctest.h"

typedef double NT;
typedef Cartesian<NT> Kernel;
typedef typename Kernel::Point Point;
typedef BoostRandomNumberGenerator<boost::mt19937, NT, 3> RNGType;
typedef HPolytope<Point> Polytope;

void test_gaussian_hmc_reflections(const std::string& polytope_type, RNGType& rng, unsigned int walk_length, unsigned int dim = 3, unsigned int m = 0, double seed = std::numeric_limits<double>::signaling_NaN()) {
    Polytope P;
    if (polytope_type == "cube") {
        P = generate_cube<Polytope>(dim, false, 10.0);
    } else if (polytope_type == "cross") {
        P = generate_cross<Polytope>(dim, false);
    } else if (polytope_type == "simplex") {
        P = generate_simplex<Polytope>(dim, false);
    } else if (polytope_type == "random") {
       P = random_hpoly<Polytope, boost::mt19937>(dim, m, 123456);
    }  else {
        std::cout << "Unknown polytope type: " << polytope_type << std::endl;
        return;
    }

    Point interior_point = P.ComputeInnerBall().first;
    double L = 1.0; // Step size
    unsigned int rho = 100; // Max reflections
    typename GaussianHamiltonianMonteCarloExactWalk::parameters params(L, true, rho, true);
    typename GaussianHamiltonianMonteCarloExactWalk::Walk<Polytope, RNGType> walk(P, interior_point, 1.0, rng, params);
    
    walk.apply(P, interior_point, 1.0, walk_length, rng);
    
    CHECK(P.is_in(interior_point) == -1);
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

TEST_CASE("GHMC Reflections Test - Random Polytope") {
    RNGType rng(123456);
    test_gaussian_hmc_reflections("random", rng, 10, 3, 30); // 3D polytope with 30 hyperplanes
}