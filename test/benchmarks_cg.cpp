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
    typedef VPolytope<Point> Vpolytope;

    std::cout << "Volume algorithm: Gaussian Annealing" << std::endl << std::endl;

    Hpolytope HP = gen_cube<Hpolytope>(10, false);

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

    double tstart;

    ////////////////////////////////////////////////////////////////
    /// H-Polytopes
    ///
    ///
    ///

    std::cout << "Volume estimation on H-polytopes (cube-10)" << std::endl;

    // Estimate the volume

    typedef BoostRandomNumberGenerator<boost::mt11213b, NT> RNG;

    tstart = (double)clock()/(double)CLOCKS_PER_SEC;
    std::cout << "BallWalk (cube) = "
              << volume_gaussian_annealing<GaussianBallWalk>(HP, e, walk_len) << " , ";
    std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;
    tstart = (double)clock()/(double)CLOCKS_PER_SEC;

    tstart = (double)clock()/(double)CLOCKS_PER_SEC;
    std::cout << "CDHRWalk (cube) = "
              << volume_gaussian_annealing<GaussianCDHRWalk>(HP, e, walk_len) << " , ";
    std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;

    tstart = (double)clock()/(double)CLOCKS_PER_SEC;
    std::cout << "RDHRWalk (cube) = "
              << volume_gaussian_annealing<GaussianRDHRWalk>(HP, e, walk_len) << " , ";
    std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;

    // OLD Implementation
    {
        tstart = (double)clock()/(double)CLOCKS_PER_SEC;
        vars<NT, RNGType> var2(rnum,n,10 + n/10,n_threads,err,e,0,0,0,0,0.0,rng,
                 urdist,urdist1,-1.0,false,false,false,false,false,false,false,false,false);
        vars_g<NT, RNGType> var1(n,walk_len,N,W,1,e,CheBall.second,rng,C,frac,ratio,delta,
                    false,false,false,false,false,true,false,false);
        std::cout << "old GC Ball (cube) = "
                  << volume_gaussian_annealing(HP, var1, var2, CheBall)
                  << " , ";
        std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;
    }
    {
        tstart = (double)clock()/(double)CLOCKS_PER_SEC;
        vars<NT, RNGType> var2(rnum,n,10 + n/10,n_threads,err,e,0,0,0,0,0.0,rng,
                 urdist,urdist1,-1.0,false,false,false,false,false,false,false,false,false);
        vars_g<NT, RNGType> var1(n,walk_len,N,W,1,e,CheBall.second,rng,C,frac,ratio,delta,
                    false,false,false,false,false,false,true,false);
        std::cout << "old GC CDHR (cube) = "
                  << volume_gaussian_annealing(HP, var1, var2, CheBall)
                  << " , ";
        std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;
    }
    {
        tstart = (double)clock()/(double)CLOCKS_PER_SEC;
        vars<NT, RNGType> var2(rnum,n,10 + n/10,n_threads,err,e,0,0,0,0,0.0,rng,
                 urdist,urdist1,-1.0,false,false,false,false,false,false,false,false,false);
        vars_g<NT, RNGType> var1(n,walk_len,N,W,1,e,CheBall.second,rng,C,frac,ratio,delta,
                    false,false,false,false,false,false,false,true);
        std::cout << "old GC RDHR (cube) = "
                  << volume_gaussian_annealing(HP, var1, var2, CheBall)
                  << " , ";
        std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;
    }
/*
    ////////////////////////////////////////////////////////////////
    /// V-Polytopes
    ///
    ///
    ///

    std::cout << "Volume estimation on V-polytopes (cross-10)" << std::endl;


    Vpolytope VP;
    VP = gen_cross<Vpolytope>(10, true);

    // NEW IMPLEMENTATIOM

    // Estimate the volume

    VP.init(VP.dimension(), VP.get_mat(), VP.get_vec());
    tstart = (double)clock()/(double)CLOCKS_PER_SEC;
    std::cout << "Ball (cross) = "
              << volume_gaussian_annealing<GaussianBallWalk>(VP) << " , ";
    std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;

    VP.init(VP.dimension(), VP.get_mat(), VP.get_vec());
    tstart = (double)clock()/(double)CLOCKS_PER_SEC;
    std::cout << "RDHR (cross) = "
              << volume_gaussian_annealing<GaussianRDHRWalk>(VP) << " , ";
    std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;

    VP.init(VP.dimension(), VP.get_mat(), VP.get_vec());
    tstart = (double)clock()/(double)CLOCKS_PER_SEC;
    std::cout << "CDHR (cross) = "
              << volume_gaussian_annealing<GaussianCDHRWalk>(VP) << " , ";
    std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;

    // OLD Implementation

    VP.init(VP.dimension(), VP.get_mat(), VP.get_vec());
    auto VPCheBall = VP.ComputeInnerBall();
    {
        tstart = (double)clock()/(double)CLOCKS_PER_SEC;
        vars<NT, RNGType> var2(rnum,n,10 + n/10,n_threads,err,e,0,0,0,0,0.0,rng,
                 urdist,urdist1,-1.0,false,false,false,false,false,false,false,false,false);
        vars_g<NT, RNGType> var1(n,walk_len,N,W,1,e,VPCheBall.second,rng,C,frac,ratio,delta,
                    false,false,false,false,false,true,false,false);
        VP.init(VP.dimension(), VP.get_mat(), VP.get_vec());
        std::cout << "old GC Ball (cross) = "
                  << volume_gaussian_annealing(HP, var1, var2, CheBall)
                  << " , ";
        std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;
    }
    VP.init(VP.dimension(), VP.get_mat(), VP.get_vec());
    VPCheBall = VP.ComputeInnerBall();
    {
        tstart = (double)clock()/(double)CLOCKS_PER_SEC;
        vars<NT, RNGType> var2(rnum,n,10 + n/10,n_threads,err,e,0,0,0,0,0.0,rng,
                 urdist,urdist1,-1.0,false,false,false,false,false,false,false,false,false);
        vars_g<NT, RNGType> var1(n,walk_len,N,W,1,e,VPCheBall.second,rng,C,frac,ratio,delta,
                    false,false,false,false,false,false,true,false);
        VP.init(VP.dimension(), VP.get_mat(), VP.get_vec());
        std::cout << "old GC CDHR (cross) = "
                  << volume_gaussian_annealing(HP, var1, var2, CheBall)
                  << " , ";
        std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;
    }
    VP.init(VP.dimension(), VP.get_mat(), VP.get_vec());
    VPCheBall = VP.ComputeInnerBall();
    {
        tstart = (double)clock()/(double)CLOCKS_PER_SEC;
        vars<NT, RNGType> var2(rnum,n,10 + n/10,n_threads,err,e,0,0,0,0,0.0,rng,
                 urdist,urdist1,-1.0,false,false,false,false,false,false,false,false,false);
        vars_g<NT, RNGType> var1(n,walk_len,N,W,1,e,VPCheBall.second,rng,C,frac,ratio,delta,
                    false,false,false,false,false,false,false,true);
        VP.init(VP.dimension(), VP.get_mat(), VP.get_vec());
        std::cout << "old GC RDHR (cross) = "
                  << volume_gaussian_annealing(HP, var1, var2, CheBall)
                  << " , ";
        std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;
    }
*/
    return 0;
}
