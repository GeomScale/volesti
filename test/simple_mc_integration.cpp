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
typedef HPolytope<Point> HPOLYTOPE;
typedef boost::mt19937 RNGType;
typedef BoostRandomNumberGenerator<RNGType, NT> RandomNumberGenerator;
typedef typename HPolytope<Point>::MT MT;
typedef typename HPolytope<Point>::VT VT;

NT exp_N_Dim(Point X){
	return exp(-X.squared_length()) ;
}

void call_test_simple_mc_over_hyperrectangle(){
	srand(time(0));
	std::cout << "\nTESTS FOR SIMPLE MC INTEGRATION OVER HYPER-RECTANGLES USING INBUILT RANDOM SAMPLING\n";

	std::vector<NT> ll{-1,-1};
	std::vector<NT> ul{1,1};
	Point LL(2,ll), UL(2,ul);
	simple_mc_integrate(exp_N_Dim,10000,LL,UL);

	std::vector<NT> ll1{0,0};
	std::vector<NT> ul1{2,2};
	Point LL1(2,ll1), UL1(2,ul1);
	simple_mc_integrate(exp_N_Dim,10000,LL1,UL1);


}

void call_test_simple_mc_over_polytope(){
	std::cout << "\nTESTS FOR SIMPLE MC INTEGRATION OVER H-POLYTOPES USING UNIFORM SAMPLING\n";

	// Polytope Integration Test:1 for 2D Polytope around the origin
	HPOLYTOPE HP = generate_cube<HPOLYTOPE>(2, false);
	simple_mc_polytope_integrate<BallWalk>(exp_N_Dim, HP, 100, SOB);

	// Polytope Integration Test:2 for 2D Polytope shifted to (1,1) from origin
	// std::vector<NT> origin{1,1};
	// Point newOrigin(2,origin);
	// simple_mc_polytope_integrate<BallWalk>(exp_N_Dim, HP, 100, SOB, newOrigin);

	// Polytope Integration Test:2 Shifting Polytope relative to origin
	// HPOLYTOPE HP1 = generate_cube<HPOLYTOPE>(2, false);
	// simple_mc_polytope_integrate(exp_N_Dim, HP1, 100000);
	// HPOLYTOPE HP2 = generate_cube<HPOLYTOPE>(3, false);
	// simple_mc_polytope_integrate(exp_N_Dim, HP2, 100000);
	//simple_mc_polytope_integrate(exp_N_Dim, HP1, 15000, CG);

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
