#include "doctest.h"
#include "simple_MC_integration.hpp"
#include "Eigen/Eigen"
#include <vector>
#include "cartesian_geom/cartesian_kernel.h"
#include "hpolytope.h"
#include "known_polytope_generators.h"
#include <iostream>
#include <fstream>
#include "misc.h"

typedef double NT;
typedef Cartesian<NT> Kernel;
typedef typename Kernel::Point Point;
typedef std::vector<Point> pts;
typedef HPolytope<Point> HPOLYTOPE;
typedef BoostRandomNumberGenerator<boost::mt19937, NT> RNGType;
typedef typename HPolytope<Point>::MT MT;
typedef typename HPolytope<Point>::VT VT;

NT exp_N_Dim(VT X){
	return exp( -X.squaredNorm() );
}

void call_test_simple_mc_over_hyperrectangle(){
	srand(time(0));
	std::cout << "\nTESTS FOR SIMPLE MC INTEGRATION OVER HYPER-RECTANGLES USING INBUILT RANDOM SAMPLING\n";
	VT LL(1),UL(1) ;
	LL << -1;
	UL << 1 ;
	simple_mc_integrate(&exp_N_Dim,10000,LL,UL);
}

void call_test_simple_mc_over_polytope(){
	std::cout << "\nTESTS FOR SIMPLE MC INTEGRATION OVER H-POLYTOPES USING ReHMC SAMPLING\n";

	// Polytope Integration Test:1 for 2D Polytope around the origin
	HPOLYTOPE HP = generate_cube<HPOLYTOPE>(1, false);
	HP.print();
	simple_mc_polytope_integrate(exp_N_Dim, HP, 10000);

	// Polytope Integration Test:2 Shifting Polytope relative to origin
	VT newOrigin(4);
	newOrigin << 0.2, -0.7, 0.96, -0.79;
	HPOLYTOPE HP1 = generate_cube<HPOLYTOPE>(4, false);
	simple_mc_polytope_integrate(exp_N_Dim, HP1, 15000, SOB , newOrigin);
	simple_mc_polytope_integrate(exp_N_Dim, HP1, 15000, CG);

	// Polytope Integration Test:3 Reading a HPolytope from ine file for 20 Dimensions
	// std::string fileName("cube10.ine");
	// std::ifstream inp;
	// std::vector<std::vector<NT> > Pin;
	// inp.open(fileName, std::ifstream::in);
	// read_pointset(inp,Pin);
	// HPOLYTOPE HP2(Pin);
	// SimpleMCPolytopeIntegrate(expXY20D, HP2, 15000, SOB);
	// inp.close();
}                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 

TEST_CASE("simple_mc_integration_over_hyperrectangles") {
    call_test_simple_mc_over_hyperrectangle();
}

TEST_CASE("simple_mc_integration_over_polytopes") {
    call_test_simple_mc_over_polytope();
}
