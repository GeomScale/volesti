#include "simple_MC_integration.hpp"
#include <iostream>
#include <cmath>
#include <vector>
#include <random>

typedef double NT;
typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;

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
    //SimpleMCIntegrate(&natural_exp2,1000,LL,UL);

    // Integration over Polytope and shifting to a new origin
    VT newOrigin(2);
    newOrigin << 0,0;
    SimpleMCPolytopeIntegrate(&natural_exp2,10,N);

    return 0;
}
