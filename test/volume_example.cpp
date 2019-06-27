// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

// VolEsti example

#include <unistd.h>
#include "Eigen/Eigen"
#include <fstream>
#include "volume.h"
#include "polytope_generators.h"
#include <typeinfo>

int main()
{
    typedef double NT;
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef boost::mt19937    RNGType;
    typedef HPolytope<Point> Hpolytope;
    Hpolytope HP;

    std::cout << "Testing volume of example H-polytope" << std::endl;
    
    //Create the polytope in form Ax<=b
    typedef typename Hpolytope::MT    MT;
    typedef typename Hpolytope::VT    VT;
    MT A;
    VT b;
    unsigned int m;
    unsigned int dim=10;

    A.resize(2 * dim, dim);
    b.resize(2 * dim);
    for (unsigned int i = 0; i < dim; ++i) {
        b(i) = 1;
        for (unsigned int j = 0; j < dim; ++j) {
            if (i == j) {
                A(i, j) = 1;
            } else {
                A(i, j) = 0;
            }
        }
    }
    for (unsigned int i = 0; i < dim; ++i) {
        b(i + dim) = 1;
        for (unsigned int j = 0; j < dim; ++j) {
            if (i == j) {
                A(i + dim, j) = -1;
            } else {
                A(i + dim, j) = 0;
            }
        }
    }
    HP.init(dim,A,b);

    //Compute chebychev ball
    std::pair<Point,NT> CheBall;
    CheBall = HP.ComputeInnerBall();

    //Parameter setup
    int n = HP.dimension();
    int walk_len=10 + n/10;
    int e=1;
    int rnum = std::pow(e,-2) * 400 * n * std::log(n);
    NT C=2;
    NT ratio = 1.0-1.0/(NT(n));
    int N = 500 * ((int) C) + ((int) (n * n / 2));
    int W = 4*n*n+500;
    
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(seed);
    boost::normal_distribution<> rdist(0,1);
    boost::random::uniform_real_distribution<>(urdist);
    boost::random::uniform_real_distribution<> urdist1(-1,1);
    
    vars<NT, RNGType> var1(rnum,n,walk_len,1,0,1,0,0,0,0,rng,
             urdist,urdist1,-1.0,false,false,false,false,false,false,true,false,false);
    vars_g<NT, RNGType> var2(n,walk_len,N,W,1,0.2,CheBall.second,rng,C,0.1,ratio,-1,false,
                    false,false,false,false,false,false,true,false);

    // Estimate the volume
    NT vol1 = volume(HP, var1, CheBall);
    std::cout << "Computed volume (alg.1) = " << vol1 << std::endl;
    
    NT vol2 = volume_gaussian_annealing(HP, var2, var1, CheBall);
    std::cout << "Computed volume (alg.2) = " << vol2 << std::endl;
       
    return 0;
}
