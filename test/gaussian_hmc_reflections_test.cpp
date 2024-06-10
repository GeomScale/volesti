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

void test_gaussian_hmc_reflections(const std::string& polytope_type, unsigned int walk_length) {
    Polytope P;
    unsigned int dim = 3; //dimension
    if (polytope_type == "cube") 
    {
        std::cout<<std::endl<<"CUBE"<<std::endl<<std::endl;
        P = generate_cube<Polytope>(dim, false);
    } 
    else if (polytope_type == "cross") 
    {
        P = generate_cross<Polytope>(dim, false);
    }
    else if (polytope_type == "simplex") 
    {
        P = generate_simplex<Polytope>(dim, false);
    } 
    else if (polytope_type == "randomHPoly") 
    {
        std::cout<<std::endl<<"RANDOM HPPOLY"<<std::endl<<std::endl;
        unsigned int m = 100; //number of hyperplanes
        P = random_hpoly<Polytope, boost::mt19937>(dim, m, 123456);
    }
    else {
        std::cout << "Unknown polytope type: " << polytope_type << std::endl;
        return;
    }
    
    Point interior_point = P.ComputeInnerBall().first;
    double L = 100.0; // Step size
    unsigned int rho = 100; // Max reflections
    RNGType rng(123456);
    //typename GaussianHamiltonianMonteCarloExactWalk::parameters params(L, true, rho, true);
    //typename GaussianHamiltonianMonteCarloExactWalk::Walk<Polytope, RNGType> walk(P, interior_point, 1.0, rng, params);
    typename GaussianHamiltonianMonteCarloExactWalk::Walk<Polytope, RNGType> walk(P, interior_point, 1.0, rng);
    walk.apply(P, interior_point, 1.0, walk_length, rng);
}









TEST_CASE("GHMC Reflections Test - Cube") {
    unsigned int walk_length = 100;
    test_gaussian_hmc_reflections("cube", walk_length);
}



TEST_CASE("GHMC Reflections Test - Random HPolytope") {
    RNGType rng(123456);
    test_gaussian_hmc_reflections("randomHPoly", 10); // 3D polytope with 30 hyperplanes
}


TEST_CASE("GHMC Reflections Test - Simplex") {
    RNGType rng(123456);
    test_gaussian_hmc_reflections("simplex", 10);
}


TEST_CASE("GHMC Reflections Test - Cross Polytope") {
    RNGType rng(1234567);
    test_gaussian_hmc_reflections("cross", 1000);
}
