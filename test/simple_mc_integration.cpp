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

NT expXY2D(VT X){
    return exp( -pow( X(0),2 )- pow( X(1),2) );
}

NT expXY4D(VT X){
    return exp( -normSquared(X) );
}

NT expXY20D(VT X){
    return exp( -normSquared(X));
}

void call_test_norm_squared(){
    std::cout << "\nTESTS FOR RETURNING VALUE OF normSquared() FUNCTION\n";
    VT normSq(7);
    normSq << 2,24,36,4,23,43,53;
    std::cout << normSquared(normSq) << "\n";
    VT normSq1(4);
    normSq1 << 21,2,32,4;
    std::cout << normSquared(normSq1) << "\n";
}

void call_test_legitimate_limits(){
    std::cout << "\nTESTS FOR LEGITIMATE LIMITS\n"; 
    VT LL0(6),UL0(6);
    LL0 << 1.2, 3.4, -5.8, 1.9, -2.01, 2.8;
    UL0 << 1.3, 3.9, 5.89, 3.1, 2.5, 2.91;
    std::cout << legitLimits(LL0,UL0) <<"\n";

    VT LL1(6),UL1(6);
    LL1 << 1.3, -3.9, 5.89, -3.1, 2.5, 2.91;
    UL1 << 1.2, 3.4, 5.8, 1.9, 2.01, 2.8;
    std::cout << legitLimits(LL1,UL1) <<"\n";

    VT LL2(4),UL2(6);
    LL2 << 1.2, 3.4, 5.8, 1.9;
    UL2 << 1.3, 3.9, 5.89, 3.1, 2.5, 2.91;
    std::cout << legitLimits(LL2,UL2) <<"\n";
}

void call_test_hyperrectangle_volume(){
    std::cout << "\nTESTS FOR HYPER-RECTANGULAR VOLUME BETWEEN TWO N-DIMENSIONAL POINTS\n";
    VT LL0(6),UL0(6);
    LL0 << 1.2, 3.4, -5.8, 1.9, -2.01, 2.8;
    UL0 << 1.3, 3.9, 5.89, 3.1, 2.5, 2.91;
    std::cout << hyperRectVolume(LL0,UL0)<<"\n";

    VT LL1(6),UL1(6);
    LL1 << 1.3, -3.9, 5.89, -3.1, 2.5, 2.91;
    UL1 << 1.2, 3.4, 5.8, 1.9, 2.01, 2.8;
    std::cout << hyperRectVolume(LL1,UL1)<<"\n";
}

void call_test_sampler_bw_two_points(){
    std::cout << "\nTESTS FOR SAMPLER BETWEEN N-DIMENSIONAL POINTS\n";
    VT LL(6),UL(6);
    LL << 1.2, 3.4, -5.8, 1.9, -2.01, 2.8;
    UL << 1.3, 3.9, 5.89, 3.1, 2.5, 2.91;
    for(int i =0 ; i<20 ; i++) samplerBWLimits(LL,UL);
}

void call_test_simple_mc_over_hyperrectangle(){
    srand(time(0));
    std::cout << "\nTESTS FOR SIMPLE MC INTEGRATION OVER HYPER-RECTANGLES USING INBUILT RANDOM SAMPLING\n";
    VT LL(2),UL(2) ;
    LL << -1.5, 0.5;
    UL << 2.5 , 3.6;
    SimpleMCIntegrate(expXY2D,10000,LL,UL);
}

void call_test_simple_mc_over_polytope(){
    std::cout << "\nTESTS FOR SIMPLE MC INTEGRATION OVER H-POLYTOPES USING ReHMC SAMPLING\n";

    // Polytope Integration Test:1 for 2D Polytope around the origin
    HPOLYTOPE HP = generate_cube<HPOLYTOPE>(2, false);
    SimpleMCPolytopeIntegrate(expXY2D, HP, 15000);

    // Polytope Integration Test:2 Shifting Polytope relative to origin
    VT newOrigin(4);
    newOrigin << 0.2, -0.7, 0.96, -0.79;
    HPOLYTOPE HP1 = generate_cube<HPOLYTOPE>(4, false);
    SimpleMCPolytopeIntegrate(expXY4D, HP1, 15000, SOB , newOrigin);
    SimpleMCPolytopeIntegrate(expXY4D, HP1, 15000, CG);

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

TEST_CASE("norm_squared") {
    call_test_norm_squared();
}

TEST_CASE("legitimate_integration_limits") {
    call_test_legitimate_limits();
}

TEST_CASE("volume_of_hyperrectangle") {
    call_test_hyperrectangle_volume();
}

TEST_CASE("sampler_between_two_points") {
    call_test_sampler_bw_two_points();
}

TEST_CASE("simple_mc_integration_over_hyperrectangles") {
    call_test_simple_mc_over_hyperrectangle();
}

TEST_CASE("simple_mc_integration_over_polytopes") {
    call_test_simple_mc_over_polytope();
}
