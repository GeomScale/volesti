#include "simple_MC_integration.hpp"
#include <iostream>
#include <cmath>
#include <vector>
#include <random>

typedef double NT;
typedef Cartesian<NT> Kernel;
typedef typename Kernel::Point Point;
typedef std::vector<Point> pts;
typedef HPolytope<Point> HPOLYTOPE;
typedef BoostRandomNumberGenerator<boost::mt19937, NT> RNGType;
typedef typename HPolytope<Point>::MT MT;
typedef typename HPolytope<Point>::VT VT;

double natural_exp2(VT X){
    return exp( -(pow( X(0),2 )+ pow( X(1),2) ));
}

int main(){

    srand(time(0)); // Mandatory code DO NOT REMOVE

    int N;
    std::cin >> N;

    // Limits of integration in VT format
    VT LL(2),UL(2);
    LL << -1, -1;  // Lower limits of integration
    UL << 1, 1;  // Upper limits of integration

    // Integration function over hyperrectangle
    SimpleMCIntegrate(natural_exp2,1000,LL,UL);
 
    // Creating a HPolytope and calculating the volume
    HPOLYTOPE HP_cube = generate_cube<HPOLYTOPE>(10, false);
    HPOLYTOPE HP_simplex = generate_simplex<HPOLYTOPE>(5, false);

    // Integration over Polytope and shifting to a new origin
    VT newOrigin(2);
    newOrigin << 0,0;
    SimpleMCPolytopeIntegrate(natural_exp2,HP_cube,1000,newOrigin,CG);

    return 0;
}
