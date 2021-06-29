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
	// Experiment to check integration values for exp_normsq over [-1,1]^n for n=1,2,...,10
	for(int i=1 ; i<=10 ; i++){
		Point LL1(i), UL1(i);
		for( int j=0; j<i ; j++) LL1.set_coord(j,-1);
		for( int j=0; j<i ; j++) UL1.set_coord(j,1);
		std::cout << "For dimension " << i << std::endl;
		integration_value = simple_mc_integrate<BilliardWalk>(exp_normsq,100000, LL1, UL1, 10, 0.01);
		std::cout << "Integration value by Billiard Walks: " << integration_value << std::endl;
		NT sum=0; Uint N=100000;
		for(int i=0; i < N ; i++){
			sum=sum+exp_normsq(samplerBWLimits(LL1,UL1));
		}
		std::cout << "Integration value of C++ <random> sampling: " << sum/N * hyper_rect_volume(LL1,UL1) << std::endl;
	}

	std::cout << "\nTESTS FOR SIMPLE MC INTEGRATION OVER DEFINED INTEGRATION LIMITS USING UNIFORM RANDOM WALKS\n";

	std::vector<NT> ll{-1,-1,-1};
	std::vector<NT> ul{1,1,1};
	Point LL(3,ll), UL(3,ul);
	integration_value = simple_mc_integrate(exp_normsq, 100000, LL, UL);
	test_values(integration_value, 3.31, 3.332);

	std::vector<NT> ll1{-1,-1,-1,-1,-1};
	std::vector<NT> ul1{1,1,1,1,1};
	Point LL1(5,ll1), UL1(5,ul1);
	integration_value = simple_mc_integrate(exp_normsq, 100000, LL1, UL1);
	test_values(integration_value, 7.49, 7.46);

	std::vector<NT> ll2{-1};
	std::vector<NT> ul2{6};
	Point LL2(1, ll2), UL2(1, ul2);
	integration_value = simple_mc_integrate(simple_polynomial_1D, 100000, LL2, UL2);
	test_values(integration_value, 39.7, 40.25);

	std::vector<NT> ll3{0.5};
	std::vector<NT> ul3{10};
	Point LL3(1, ll3), UL3(1, ul3);
	integration_value = simple_mc_integrate(logx_natural_1D, 100000, LL3, UL3);
	test_values(integration_value, 13.65, 13.872);

	std::vector<NT> ll4{-1,-1};
	std::vector<NT> ul4{1,1};
	Point LL4(2, ll4), UL4(2, ul4);
	integration_value = simple_mc_integrate(rooted_squaresum_2D, 100000, LL4, UL4);
	test_values(integration_value, 2.99, 3.0607);

}

void call_test_simple_mc_over_polytope(){

	std::cout << "\nTESTS FOR SIMPLE MC INTEGRATION OVER CONVEX-BODIES USING UNIFORM RANDOM WALKS\n";
	NT integration_value;
	// H-Polytope Integration Test:1 for 2D Polytope around the origin
	HPOLYTOPE HP = generate_cube<HPOLYTOPE>(2, false);
	integration_value = simple_mc_polytope_integrate<BilliardWalk,HPOLYTOPE>(exp_normsq, HP, 100000, SOB, 10, 0.01);
	test_values(integration_value, 2.20, 2.230);

	// H-Polytope Integration Test:2 for 2D Polytope shifted to (1,1) from origin
	std::vector<NT> origin{1,1};
	Point newOrigin(2,origin);
	integration_value = simple_mc_polytope_integrate<BilliardWalk,HPOLYTOPE>(exp_normsq, HP, 100000, SOB, 1, 0.01, newOrigin);
	test_values(integration_value, 0.78, 0.777);

	// H-Polytope Integration Test:3
	HPOLYTOPE HP1 = generate_cube<HPOLYTOPE>(6, false);
	integration_value = simple_mc_polytope_integrate<BilliardWalk,HPOLYTOPE>(exp_normsq, HP1, 100000, SOB);
	test_values(integration_value, 10.9, 11.102);

	// H-Polytope Integration Test:4
	HPOLYTOPE HP2 = generate_cube<HPOLYTOPE>(8, false);
	integration_value = simple_mc_polytope_integrate<BilliardWalk,HPOLYTOPE>(exp_normsq, HP2, 100000, SOB);
	test_values(integration_value, 24.9, 24.76);

	// H-Polytope Integration Test:5
	HPOLYTOPE HP3 = generate_cube<HPOLYTOPE>(10, false);
	integration_value = simple_mc_polytope_integrate<BilliardWalk,HPOLYTOPE>(exp_normsq, HP3, 100000, SOB);
	test_values(integration_value, 54.7, 55.25);

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

TEST_CASE("simple_mc_integration_over_hyperrectangles") {
    call_test_simple_mc_over_hyperrectangle();
}

TEST_CASE("simple_mc_integration_over_polytopes") {
    call_test_simple_mc_over_polytope();
}
