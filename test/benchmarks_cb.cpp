// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

// VolEsti example

#include "Eigen/Eigen"
//#define VOLESTI_DEBUG
#include <fstream>
#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>

#include "random_walks/random_walks.hpp"

#include "volume/volume_sequence_of_balls.hpp"
#include "volume/volume_cooling_gaussians.hpp"
#include "volume/volume_cooling_balls.hpp"

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

    std::cout << "Volume algorithm: Cooling Balls" << std::endl << std::endl;

    Hpolytope HP = generate_cube<Hpolytope>(10, false);

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
              << volume_cooling_balls<BallWalk, RNG>(HP, e, walk_len) << " , ";
    std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;
    tstart = (double)clock()/(double)CLOCKS_PER_SEC;

    tstart = (double)clock()/(double)CLOCKS_PER_SEC;
    std::cout << "RDHRWalk (cube) = "
              << volume_cooling_balls<RDHRWalk, RNG>(HP, e, walk_len) << " , ";
    std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;

    tstart = (double)clock()/(double)CLOCKS_PER_SEC;
    std::cout << "CDHRWalk (cube) = "
              << volume_cooling_balls<CDHRWalk, RNG>(HP, e, walk_len) << " , ";
    std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;

    tstart = (double)clock()/(double)CLOCKS_PER_SEC;
    std::cout << "BirdWalk (cube) = "
              << volume_cooling_balls<BilliardWalk, RNG>(HP, e, walk_len) << " , ";
    std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;

#ifdef VOLESTI_OLD_IMPLEMENTATION

    // OLD Implementation

    NT diameter;
    HP.comp_diam(diameter, CheBall.second);

    NT lb = 0.1, ub = 0.15, p = 0.75, rmax = 0.0, alpha = 0.2;
    int W2 = 500, NNu = 150, nu =10;
    bool win2 = false;
    vars_ban <NT> var_ban(lb, ub, p, rmax, alpha, W2, NNu, nu, win2);

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
                              false,false,false,true,false,false);

        tstart = (double)clock()/(double)CLOCKS_PER_SEC;
        std::cout << "OLD RDHR = " << vol_cooling_balls(HP, var, var_ban, CheBall)
                  << " , ";
        std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;
    }
    {
        vars<NT, RNGType> var(rnum,n,walk_len,n_threads,err,e,0,0,0,
                              CheBall.second,diameter,rng,
                              urdist,urdist1,-1.0,false,false,false,
                              false,false,false,false,true,false);

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
#endif

    ////////////////////////////////////////////////////////////////
    /// V-Polytopes
    ///
    ///
    ///

    std::cout << "Volume estimation on V-polytopes (cross-10)" << std::endl;


    Vpolytope VP;
    VP = std::move(generate_cross<Vpolytope>(10, true));

    // NEW IMPLEMENTATIOM

    // Estimate the volume

    //VP.init(VP.dimension(), VP.get_mat(), VP.get_vec());
    tstart = (double)clock()/(double)CLOCKS_PER_SEC;
    std::cout << "Ball (cross) = "
              << volume_cooling_balls<BallWalk, RNG>(VP) << " , ";
    std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;

    //VP.init(VP.dimension(), VP.get_mat(), VP.get_vec());
    tstart = (double)clock()/(double)CLOCKS_PER_SEC;
    std::cout << "RDHR (cross) = "
              << volume_cooling_balls<RDHRWalk, RNG>(VP) << " , ";
    std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;

    //VP.init(VP.dimension(), VP.get_mat(), VP.get_vec());
    tstart = (double)clock()/(double)CLOCKS_PER_SEC;
    std::cout << "CDHR (cross) = "
              << volume_cooling_balls<CDHRWalk, RNG>(VP) << " , ";
    std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;

    // VP.init(VP.dimension(), VP.get_mat(), VP.get_vec());
    tstart = (double)clock()/(double)CLOCKS_PER_SEC;
    std::cout << "Blrd (cross) = "
              << volume_cooling_balls<BilliardWalk, RNG>(VP) << " , ";
    std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;

#ifdef VOLESTI_OLD_IMPLEMENTATION

    // OLD Implementation

    Vpolytope VP(VP.dimension(), VP.get_mat(), VP.get_vec());
    auto VPCheBall = VP.ComputeInnerBall();
    {
        vars<NT, RNGType> var(rnum,n,walk_len,n_threads,err,e,0,0,0,
                              VPCheBall.second,diameter,rng,
                              urdist,urdist1,-1.0,false,false,false,
                              false,false,true,false,false,false);

        VP.init(VP.dimension(), VP.get_mat(), VP.get_vec());
        tstart = (double)clock()/(double)CLOCKS_PER_SEC;
        std::cout << "OLD Ball = " << vol_cooling_balls(VP, var, var_ban, VPCheBall)
                  << " , ";
        std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;
    }
    Vpolytope VP(VP.dimension(), VP.get_mat(), VP.get_vec());
    VPCheBall = VP.ComputeInnerBall();
    {
        vars<NT, RNGType> var(rnum,n,walk_len,n_threads,err,e,0,0,0,
                              VPCheBall.second,diameter,rng,
                              urdist,urdist1,-1.0,false,false,false,
                              false,false,false,true,false,false);

        Vpolytope VP(VP.dimension(), VP.get_mat(), VP.get_vec());
        tstart = (double)clock()/(double)CLOCKS_PER_SEC;
        std::cout << "OLD CDHR = " << vol_cooling_balls(VP, var, var_ban, VPCheBall)
                  << " , ";
        std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;
    }
    Vpolytope VP(VP.dimension(), VP.get_mat(), VP.get_vec());
    VPCheBall = VP.ComputeInnerBall();
    {
        vars<NT, RNGType> var(rnum,n,walk_len,n_threads,err,e,0,0,0,
                              VPCheBall.second,diameter,rng,
                              urdist,urdist1,-1.0,false,false,false,
                              false,false,false,false,true,false);

        Vpolytope VP(VP.dimension(), VP.get_mat(), VP.get_vec());
        tstart = (double)clock()/(double)CLOCKS_PER_SEC;
        std::cout << "OLD RDHR = " << vol_cooling_balls(VP, var, var_ban, VPCheBall)
                  << " , ";
        std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;
    }
    Vpolytope VP(VP.dimension(), VP.get_mat(), VP.get_vec());
    VPCheBall = VP.ComputeInnerBall();
    {
        vars<NT, RNGType> var(rnum,n,walk_len,n_threads,err,e,0,0,0,
                              VPCheBall.second,diameter,rng,
                              urdist,urdist1,-1.0,false,false,false,
                              false,false,false,false,false,true);

        Vpolytope VP(VP.dimension(), VP.get_mat(), VP.get_vec());
        tstart = (double)clock()/(double)CLOCKS_PER_SEC;
        std::cout << "OLD Blrd = " << vol_cooling_balls(VP, var, var_ban, VPCheBall)
                  << " , ";
        std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;
    }

#endif

    return 0;
}
