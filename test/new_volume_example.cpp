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
#include "new_gaussian_volume.hpp"
#include "new_cooling_balls.hpp"
#include "rotating.h"
#include "misc.h"
#include "linear_extensions.h"
#include "cooling_balls.h"
//#include "cooling_hpoly.h"
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

    HP = gen_cube<Hpolytope>(40, false);

    //Compute chebychev ball
    std::pair<Point,NT> CheBall;
    CheBall = HP.ComputeInnerBall();

    // Setup the parameters
    int n = HP.dimension();
    int walk_len=10 + n/10;
    int n_threads=1;
    NT e=1, err=0.1;
    NT C=2.0,ratio,frac=0.1,delta=-1.0;
    int N = 500 * ((int) C) + ((int) (n * n / 2));
    int W = 6*n*n+800;
    ratio = 1.0-1.0/(NT(n));

    int rnum = std::pow(e,-2) * 400 * n * std::log(n);
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(seed);
    boost::normal_distribution<> rdist(0,1);
    boost::random::uniform_real_distribution<>(urdist);
    boost::random::uniform_real_distribution<> urdist1(-1,1);



    Hpolytope P = gen_skinny_cube<Hpolytope>(5);
    P.ComputeInnerBall();

    //Hpolytope P40 = gen_cube<Hpolytope>(40, false);
    //CheBall = P40.ComputeInnerBall();

    // Estimate the volume
    double tstart;
    //std::cout << "Default (cube) = " << volume(HP) << std::endl;
    //std::cout << "Default (cube) = " << volume(HP, 0.5) << std::endl;
    //std::cout << "Default (cube) = " << volume(HP, 0.5, 2) << std::endl;

    //typedef BoostRandomNumberGenerator<boost::mt19937, NT> RNG;
    typedef BoostRandomNumberGenerator<boost::mt11213b, NT> RNG;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT, 3> RNGFixed;
    //RNG boost_rng(HP.dimension());
    //RNG3 boost_rng3(HP.dimension());
/*
    tstart = (double)clock()/(double)CLOCKS_PER_SEC;
    std::cout << "BallWalk (cube) = "
              << volume_sequence_of_balls<BallWalk, RNG>(HP, e, walk_len) << " , ";
    std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;tstart = (double)clock()/(double)CLOCKS_PER_SEC;
    tstart = (double)clock()/(double)CLOCKS_PER_SEC;
    std::cout << "BallWalk (cube) = "
              << volume_sequence_of_balls<BallWalk, RNGFixed>(HP, e, walk_len) << " , ";
    std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;tstart = (double)clock()/(double)CLOCKS_PER_SEC;
    tstart = (double)clock()/(double)CLOCKS_PER_SEC;
    std::cout << "BallWalk (cube) = "
              << volume_sequence_of_balls<BallWalk>(HP, e, walk_len) << " , ";
    std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;tstart = (double)clock()/(double)CLOCKS_PER_SEC;

    tstart = (double)clock()/(double)CLOCKS_PER_SEC;
    std::cout << "CDHRWalk (cube) = "
              << volume_sequence_of_balls<CDHRWalk>(HP, e, walk_len) << " , ";
    std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;
    tstart = (double)clock()/(double)CLOCKS_PER_SEC;
    std::cout << "default (cube) = "
              << volume_sequence_of_balls<>(HP, e, walk_len) << " , ";
    std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;

    tstart = (double)clock()/(double)CLOCKS_PER_SEC;
    std::cout << "RDHRWalk (cube) = "
              << volume_sequence_of_balls<RDHRWalk>(HP, e, walk_len) << " , ";
    std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;

    tstart = (double)clock()/(double)CLOCKS_PER_SEC;
    std::cout << "BilliardWalk (cube) = "
              << volume_sequence_of_balls<BilliardWalk>(HP, e, walk_len) << " , ";
    std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;
    std::cout << std::endl;
*/

/*
 * Gaussian annealing
 *
    tstart = (double)clock()/(double)CLOCKS_PER_SEC;
    std::cout << "new GC Ball (cube) = "
              << volume_gaussian_annealing<GaussianBallWalk>(HP, e, walk_len) << " , ";
    std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;

    tstart = (double)clock()/(double)CLOCKS_PER_SEC;
    std::cout << "new GC RDHR (cube) = "
              << volume_gaussian_annealing<GaussianRDHRWalk>(HP, e, walk_len) << " , ";
    std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;

    tstart = (double)clock()/(double)CLOCKS_PER_SEC;
    std::cout << "new GC CDHR (cube) = "
              << volume_gaussian_annealing<GaussianCDHRWalk>(HP, e, walk_len) << " , ";
    std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;
*/

/*
    {
        tstart = (double)clock()/(double)CLOCKS_PER_SEC;
        vars<NT, RNGType> var(rnum,n,walk_len,n_threads,err,e,0,0,0,0,0.0,rng,
                 urdist,urdist1,-1.0,false,false,false,false,false,true,false,false,false);
        std::cout << "OLD Ball (cube) = " << volume(HP, var, CheBall) << " , ";
        std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;
    }
    {
        tstart = (double)clock()/(double)CLOCKS_PER_SEC;
        vars<NT, RNGType> var(rnum,n,walk_len,n_threads,err,e,0,0,0,0,0.0,rng,
                 urdist,urdist1,-1.0,false,false,false,false,false,false,true,false,false);
        std::cout << "OLD CDHR (cube) = " << volume(HP, var, CheBall) << " , ";
        std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;
    }
    {
        tstart = (double)clock()/(double)CLOCKS_PER_SEC;
        vars<NT, RNGType> var(rnum,n,walk_len,n_threads,err,e,0,0,0,0,0.0,rng,
                 urdist,urdist1,-1.0,false,false,false,false,false,false,false,true,false);
        std::cout << "OLD RDHR (cube) = " << volume(HP, var, CheBall) << " , ";
        std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;
    }
    {
        tstart = (double)clock()/(double)CLOCKS_PER_SEC;
        vars<NT, RNGType> var(rnum,n,walk_len,n_threads,err,e,0,0,0,0,0.0,rng,
                 urdist,urdist1,-1.0,false,false,false,false,false,false,false,false,true);
        //std::cout << "OLD (cube) = " << volume(HP, var, CheBall) << " , ";
        std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;
    }*/
    /*
    {
        tstart = (double)clock()/(double)CLOCKS_PER_SEC;
        vars<NT, RNGType> var2(rnum,n,10 + n/10,n_threads,err,e,0,0,0,0,0.0,rng,
                 urdist,urdist1,-1.0,false,false,false,false,false,true,false,false,false);
        vars_g<NT, RNGType> var1(n,walk_len,N,W,1,e,CheBall.second,rng,C,frac,ratio,delta,
                    false,false,false,false,true,false,false,false);
        std::cout << "old GC Ball (cube) = "
                  << volume_gaussian_annealing(HP, var1, var2, CheBall)
                  << " , ";
        std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;
    }
    {
        tstart = (double)clock()/(double)CLOCKS_PER_SEC;
        vars<NT, RNGType> var2(rnum,n,10 + n/10,n_threads,err,e,0,0,0,0,0.0,rng,
                 urdist,urdist1,-1.0,false,false,false,false,false,false,true,false,false);
        vars_g<NT, RNGType> var1(n,walk_len,N,W,1,e,CheBall.second,rng,C,frac,ratio,delta,
                    false,false,false,false,false,true,false,false);
        std::cout << "old GC CDHR (cube) = "
                  << volume_gaussian_annealing(HP, var1, var2, CheBall)
                  << " , ";
        std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;
    }
    {
        tstart = (double)clock()/(double)CLOCKS_PER_SEC;
        vars<NT, RNGType> var2(rnum,n,10 + n/10,n_threads,err,e,0,0,0,0,0.0,rng,
                 urdist,urdist1,-1.0,false,false,false,false,false,false,false,true,false);
        vars_g<NT, RNGType> var1(n,walk_len,N,W,1,e,CheBall.second,rng,C,frac,ratio,delta,
                    false,false,false,false,false,false,true,false);
        std::cout << "old GC RDHR (cube) = "
                  << volume_gaussian_annealing(HP, var1, var2, CheBall)
                  << " , ";
        std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;
    }
*/
    NT diameter;
    HP.comp_diam(diameter, CheBall.second);


    NT lb = 0.1, ub = 0.15, p = 0.75, rmax = 0.0, alpha = 0.2;
    int W2 = 500, NNu = 150, nu =10;
    bool win2 = false;
    vars_ban <NT> var_ban(lb, ub, p, rmax, alpha, W2, NNu, nu, win2);

    vars<NT, RNGType> var(rnum,n,walk_len,n_threads,err,e,0,0,0,
                          CheBall.second,diameter,rng,
                          urdist,urdist1,-1.0,false,false,false,
                          false,false,false,false,false,true);


    {
        vars<NT, RNGType> var(rnum,n,walk_len,n_threads,err,e,0,0,0,
                              CheBall.second,diameter,rng,
                              urdist,urdist1,-1.0,false,false,false,
                              false,false,true,false,false,false);

        tstart = (double)clock()/(double)CLOCKS_PER_SEC;
        std::cout << "OLD Ball = " << vol_cooling_balls(HP, var, var_ban, CheBall)
                  << " , ";
        std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;
    }
    {
        vars<NT, RNGType> var(rnum,n,walk_len,n_threads,err,e,0,0,0,
                              CheBall.second,diameter,rng,
                              urdist,urdist1,-1.0,false,false,false,
                              false,false,false,false,true,false);

        tstart = (double)clock()/(double)CLOCKS_PER_SEC;
        std::cout << "OLD RDHR = " << vol_cooling_balls(HP, var, var_ban, CheBall)
                  << " , ";
        std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;
    }
    {
        vars<NT, RNGType> var(rnum,n,walk_len,n_threads,err,e,0,0,0,
                              CheBall.second,diameter,rng,
                              urdist,urdist1,-1.0,false,false,false,
                              false,false,false,true,false,false);

        tstart = (double)clock()/(double)CLOCKS_PER_SEC;
        std::cout << "OLD CDHR = " << vol_cooling_balls(HP, var, var_ban, CheBall)
                  << " , ";
        std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;
    }
    {
        vars<NT, RNGType> var(rnum,n,walk_len,n_threads,err,e,0,0,0,
                              CheBall.second,diameter,rng,
                              urdist,urdist1,-1.0,false,false,false,
                              false,false,false,false,false,true);

        tstart = (double)clock()/(double)CLOCKS_PER_SEC;
        std::cout << "OLD Blrd = " << vol_cooling_balls(HP, var, var_ban, CheBall)
                  << " , ";
        std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;
    }

    tstart = (double)clock()/(double)CLOCKS_PER_SEC;
    std::cout << "NEW Ball = "
              << volume_cooling_balls<BallWalk>(HP, e, walk_len)
              << " , ";
    std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;

    tstart = (double)clock()/(double)CLOCKS_PER_SEC;
    std::cout << "NEW RDHR = "
              << volume_cooling_balls<RDHRWalk>(HP, e, walk_len)
              << " , ";
    std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;

    tstart = (double)clock()/(double)CLOCKS_PER_SEC;
    std::cout << "NEW CDHR = "
              << volume_cooling_balls<CDHRWalk>(HP, e, walk_len)
              << " , ";
    std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;

    tstart = (double)clock()/(double)CLOCKS_PER_SEC;
    std::cout << "NEW Blrd = "
              << volume_cooling_balls<BilliardWalk>(HP, e, walk_len)
              << " , ";
    std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;

    return 0;
}
