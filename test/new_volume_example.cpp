// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

// VolEsti example

#include "Eigen/Eigen"
//#define VOLESTI_DEBUG
#include <fstream>
#include "random.hpp"
#include "random/uniform_int.hpp"
#include "random/normal_distribution.hpp"
#include "random/uniform_real_distribution.hpp"
#include "volume.h"
#include "new_volume.hpp"
#include "rotating.h"
#include "misc.h"
#include "linear_extensions.h"
#include "cooling_balls.h"
#include "cooling_hpoly.h"
#include "sample_only.h"
#include "exact_vols.h"
#include "generators/known_polytope_generators.h"

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

    // Setup the parameters
    int n = HP.dimension();
    int walk_len=10 + n/10;
    int nexp=1, n_threads=1;
    NT e=1, err=0.0000000001;
    int rnum = std::pow(e,-2) * 400 * n * std::log(n);
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(seed);
    boost::normal_distribution<> rdist(0,1);
    boost::random::uniform_real_distribution<>(urdist);
    boost::random::uniform_real_distribution<> urdist1(-1,1);



    Hpolytope P = gen_skinny_cube<Hpolytope>(5);
    P.ComputeInnerBall();

    // Estimate the volume
    NT vol1;
    double tstart, tstop;
    //std::cout << "Default (cube) = " << volume(HP) << std::endl;
    //std::cout << "Default (cube) = " << volume(HP, 0.5) << std::endl;
    //std::cout << "Default (cube) = " << volume(HP, 0.5, 2) << std::endl;

    typedef BoostRandomNumberGenerator<boost::mt19937, NT> RNG;
    RNG boost_rng(HP.dimension()-1, seed);

    tstart = (double)clock()/(double)CLOCKS_PER_SEC;
    std::cout << "BallWalk (cube) = "
              << volume<BallWalk<Hpolytope,RNG>>(HP, boost_rng, e, walk_len) << " , ";
    std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;

    tstart = (double)clock()/(double)CLOCKS_PER_SEC;
    std::cout << "CDHRWalk (cube) = "
              << volume<CDHRWalk<Hpolytope,RNG>>(HP, boost_rng, e, walk_len) << " , ";
    std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;

    tstart = (double)clock()/(double)CLOCKS_PER_SEC;
    std::cout << "RDHRWalk (cube) = "
              << volume<RDHRWalk<Hpolytope,RNG>>(HP, boost_rng, e, walk_len) << " , ";
    std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;

    tstart = (double)clock()/(double)CLOCKS_PER_SEC;
    std::cout << "BilliardWalk (cube) = "
              << volume<BilliardWalk<Hpolytope,RNG>>(HP, boost_rng, e, walk_len) << " , ";
    std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;

    //std::cout << "Default (cube) = " << volume_old(P, var, 1.0, walk_len, BilliardWalkOld<Point>()) << std::endl;
/*    vol1 = volume(HP, 1, walk_len, BallWalk<Point>());
    std::cout << "Ball (cube) = " << vol1 << std::endl;
    vol1 = volume(P, 1, walk_len, BallWalk<Point>());
    std::cout << "     (skinny) = " << vol1 << std::endl;

    vol1 = volume(HP, 1, walk_len, RDHRWalk<Point>());
    std::cout << "RDHNR (cube) = " << vol1 << std::endl;
    vol1 = volume(P, 1, walk_len, RDHRWalk<Point>());
    std::cout << "      (skinny) = " << vol1 << std::endl;

    vol1 = volume(HP, 1, walk_len, CDHRWalk<Point>());
    std::cout << "CDHNR (cube) = " << vol1 << std::endl;
    vol1 = volume(P, 1, walk_len, CDHRWalk<Point>());
    std::cout << "      (skinny) = " << vol1 << std::endl;

    vol1 = volume(HP, 1, walk_len, BilliardWalk<Point>());
    std::cout << "Billiard (cube) = " << vol1 << std::endl;
    //vol1 = volume(P, 1, walk_len, BilliardWalk<Point>());
    std::cout << "         (skinny) = " << vol1 << std::endl;
*/
    {
        tstart = (double)clock()/(double)CLOCKS_PER_SEC;
        vars<NT, RNGType> var(rnum,n,walk_len,n_threads,err,e,0,0,0,0,0.0,rng,
                 urdist,urdist1,-1.0,false,false,false,false,false,true,false,false,false);
        std::cout << "OLD (cube) = " << volume(HP, var, CheBall) << " , ";
        std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;
    }
    {
        tstart = (double)clock()/(double)CLOCKS_PER_SEC;
        vars<NT, RNGType> var(rnum,n,walk_len,n_threads,err,e,0,0,0,0,0.0,rng,
                 urdist,urdist1,-1.0,false,false,false,false,false,false,true,false,false);
        std::cout << "OLD (cube) = " << volume(HP, var, CheBall) << " , ";
        std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;
    }
    {
        tstart = (double)clock()/(double)CLOCKS_PER_SEC;
        vars<NT, RNGType> var(rnum,n,walk_len,n_threads,err,e,0,0,0,0,0.0,rng,
                 urdist,urdist1,-1.0,false,false,false,false,false,false,false,true,false);
        std::cout << "OLD (cube) = " << volume(HP, var, CheBall) << " , ";
        std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;
    }
    {
        tstart = (double)clock()/(double)CLOCKS_PER_SEC;
        vars<NT, RNGType> var(rnum,n,walk_len,n_threads,err,e,0,0,0,0,0.0,rng,
                 urdist,urdist1,-1.0,false,false,false,false,false,false,false,false,true);
        //std::cout << "OLD (cube) = " << volume(HP, var, CheBall) << " , ";
        std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;
    }

    return 0;
}
