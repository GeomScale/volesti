#include "doctest.h"
#include "simple_MC_integration.hpp"
#include "Eigen/Eigen"
#include <vector>
#include "cartesian_geom/cartesian_kernel.h"
#include "hpolytope.h"
#include "known_polytope_generators.h"
#include "ode_solvers/oracle_functors.hpp"
#include "random_walks/random_walks.hpp"
#include <iostream>
#include <fstream>
#include "misc.h"

typedef double NT;
typedef Cartesian<NT> Kernel;
typedef typename Kernel::Point Point;
typedef std::vector<Point> Points;
typedef HPolytope<Point> HPOLYTOPE;
typedef VPolytope<Point> VPOLYTOPE;
typedef boost::mt19937 RNGType;
typedef BoostRandomNumberGenerator<RNGType, NT> RandomNumberGenerator;

NT exp_normsq(Point X){
	return exp(-X.squared_length()) ;
}

NT simple_polynomial_1D(Point X){
	return (X[0]-1) * (X[0]-2) * (X[0]-3);
}

NT logx_natural_1D(Point X){
	return log(X[0]);
}

NT rooted_squaresum_2D(Point X){
	return sqrt(X[0]*X[0]+X[1]*X[1]);
}

template <typename NT>
void test_values(NT computed, NT expected, NT exact)
{
    std::cout << "Computed integration value = " << computed << std::endl;
    std::cout << "Expected integration value = " << expected << std::endl;
	std::cout << "Exact integration value = " << exact << std::endl;
    std::cout << "Relative error (expected) = "
              << std::abs((computed-expected)/expected) << std::endl;
    std::cout << "Relative error (exact) = "
              << std::abs((computed-exact)/exact) << std::endl ;
    CHECK(((std::abs((computed - expected)/expected) < 0.00001) || (std::abs((computed - exact)/exact) < 0.2)));
}

void call_test_simple_mc_over_hyperrectangle(){
	
	NT integration_value;
	std::cout << "\nTESTS FOR SIMPLE MC INTEGRATION OVER DEFINED INTEGRATION LIMITS USING UNIFORM RANDOM WALKS\n";

	integration_value = simple_mc_integrate<BallWalk>(exp_normsq, 10, 100000);
	test_values(integration_value, 54.8, 55.25);

	integration_value = simple_mc_integrate<BilliardWalk>(exp_normsq, 8, 100000);
	test_values(integration_value, 24.8, 24.76);
	
	integration_value = simple_mc_integrate<BilliardWalk>(exp_normsq, 5, 100000);
	test_values(integration_value, 7.49, 7.46);

	Limit LL{-1};
	Limit UL{6};
	integration_value = simple_mc_integrate(simple_polynomial_1D, 1, 100000, LL, UL);
	test_values(integration_value, 39.7, 40.25);

	Limit LL1{0.5};
	Limit UL1{10};
	integration_value = simple_mc_integrate(logx_natural_1D, 1, 100000, LL1, UL1);
	test_values(integration_value, 13.65, 13.872);

	Limit LL2{-1, -1};
	Limit UL2{1, 1};
	integration_value = simple_mc_integrate(rooted_squaresum_2D, 2, 100000, LL2, UL2);
	test_values(integration_value, 2.99, 3.0607);

}

void call_test_simple_mc_over_polytope(){

	std::cout << "\nTESTS FOR SIMPLE MC INTEGRATION OVER CONVEX-BODIES USING UNIFORM RANDOM WALKS\n";

	NT integration_value;
	// H-Polytope Integration Test:1 for 2D Polytope around the origin
	HPOLYTOPE HP = generate_cube<HPOLYTOPE>(2, false);
	integration_value = simple_mc_polytope_integrate<BilliardWalk, HPOLYTOPE>(exp_normsq, 2, HP, 100000, SOB, 10, 0.01);
	test_values(integration_value, 2.20, 2.230);

	// H-Polytope Integration Test:2 for 2D Polytope shifted to (1,1) from origin
	std::vector<NT> origin{1,1};
	Point newOrigin(2,origin);
	integration_value = simple_mc_polytope_integrate<BilliardWalk, HPOLYTOPE>(exp_normsq, 2, HP, 100000, SOB, 1, 0.01, newOrigin);
	test_values(integration_value, 0.78, 0.777);

	// H-Polytope Integration Test:3
	HPOLYTOPE HP1 = generate_cube<HPOLYTOPE>(10, false);
	integration_value = simple_mc_polytope_integrate<BilliardWalk, HPOLYTOPE>(exp_normsq, 10, HP1, 100000, SOB);
	test_values(integration_value, 54.7, 55.25);

	// H-Polytope Integration Test:4
	HPOLYTOPE HP2 = generate_cube<HPOLYTOPE>(15, false);
	integration_value = simple_mc_polytope_integrate<BilliardWalk, HPOLYTOPE>(exp_normsq, 15, HP2, 100000, SOB);
	test_values(integration_value, 405.9, 410.690);

	// H-Polytope Integration Test:5
	HPOLYTOPE HP3 = generate_cube<HPOLYTOPE>(20, false);
	integration_value = simple_mc_polytope_integrate<BilliardWalk, HPOLYTOPE>(exp_normsq, 20, HP3, 100000, SOB);
	test_values(integration_value, 3050.0, 3052.71);

	// V-Polytope Integration Test:6
	// VPOLYTOPE VP = generate_cross<VPOLYTOPE>(2, true);
	// integration_value = simple_mc_polytope_integrate<BilliardWalk,VPOLYTOPE>(exp_normsq, VP, 100000, CB);
	// std::cout << "Integration value: " << integration_value << std::endl;
	// test_values(integration_value,expected,exact);

	// Polytope Integration Test:3 Reading a HPolytope from ine file for 20 Dimensions
	// std::string fileName("cube10.ine");
	// std::ifstream inp;
	// std::vector<std::vector<NT> > Pin;
	// inp.open(fileName, std::ifstream::in);
	// read_pointset(inp,Pin);
	// HPOLYTOPE HP2(Pin);
	// SimpleMCPolytopeIntegrate(exp_N_dim, HP2, 15000, SOB);
	// inp.close();
}                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 

TEST_CASE("hyperrectangles") {
    call_test_simple_mc_over_hyperrectangle();
}

TEST_CASE("polytopes") {
    call_test_simple_mc_over_polytope();
}
