// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

// VolEsti example

#include "Eigen/Eigen"
//#define VOLESTI_DEBUG
#include <fstream>
#include <boost/random.hpp>
#include "random/uniform_int.hpp"
#include "random/normal_distribution.hpp"
#include "random/uniform_real_distribution.hpp"
#include "volume/volume_sequence_of_balls.hpp"
#include "volume/volume_cooling_gaussians.hpp"
#include "volume/volume_cooling_balls.hpp"
#include "volume/rotating.hpp"
#include "misc/misc.h"
#include "misc/linear_extensions.h"
#include "sampling/sampling.hpp"
#include "exact_vols.h"
#include "generators/known_polytope_generators.h"
#include "generators/z_polytopes_generators.h"

#include "random_walks/random_walks.hpp"

#include "volume/volume_sequence_of_balls.hpp"
#include "volume/volume_cooling_hpoly.hpp"
#include "volume/volume_cooling_gaussians.hpp"
#include "volume/volume_cooling_gaussians.hpp"

int main()
{
    typedef double NT;
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef boost::mt19937    RNGType;
    typedef HPolytope<Point> Hpolytope;
    typedef VPolytope<Point> Vpolytope;


    typedef BoostRandomNumberGenerator<boost::mt11213b, NT> RNG;
    Hpolytope HPoly = generate_cube<Hpolytope>(10, false);
    BallWalk BW(3);

    // Estimate the volume

    VPolytope<Point> VP2 = generate_cube<Vpolytope>(2, true);

    // VP2.init(VP2.dimension(), VP2.get_mat(), VP2.get_vec());
    double tstart = (double)clock()/(double)CLOCKS_PER_SEC;
    std::cout << "Cube-v cb = "
              << volume_cooling_balls<>(VP2) << " , ";
    std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;


    tstart = (double)clock()/(double)CLOCKS_PER_SEC;
    std::cout << "Ball SOB = "
              << volume_sequence_of_balls<>(HPoly) << " , ";
    std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;

    tstart = (double)clock()/(double)CLOCKS_PER_SEC;
    std::cout << "Ball CG = "
              << volume_cooling_gaussians<>(HPoly) << " , ";
    std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;

    tstart = (double)clock()/(double)CLOCKS_PER_SEC;
    std::cout << "Ball CB = "
              << volume_cooling_balls<>(HPoly) << " , ";
    std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;

    //
    tstart = (double)clock()/(double)CLOCKS_PER_SEC;
    std::cout << "Ball SOB = "
              << volume_sequence_of_balls<CDHRWalk, RNG>(HPoly) << " , ";
    std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;

    tstart = (double)clock()/(double)CLOCKS_PER_SEC;
    std::cout << "Ball CG = "
              << volume_cooling_gaussians<GaussianCDHRWalk, RNG>(HPoly) << " , ";
    std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;

    tstart = (double)clock()/(double)CLOCKS_PER_SEC;
    std::cout << "Ball CB = "
              << volume_cooling_balls<CDHRWalk, RNG>(HPoly) << " , ";
    std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;


    typedef Zonotope<Point> Zonotope;
    Zonotope Z = gen_zonotope_uniform<Zonotope, RNGType>(10, 15, 211);
    tstart = (double)clock()/(double)CLOCKS_PER_SEC;
    std::cout << "Zono CB = "
              << volume_cooling_hpoly<CDHRWalk, RNG, Hpolytope>(Z) << " , ";
    std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;


    return 0;
/*
    Vpolytope VP, VP2;

    VP.init(VP.dimension(), VP2.get_mat(), VP2.get_vec());
    tstart = (double)clock()/(double)CLOCKS_PER_SEC;
    std::cout << "RDHR (cross) = "
              << volume_sequence_of_balls<RDHRWalk>(VP) << " , ";
    std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;

    VP.init(VP.dimension(), VP2.get_mat(), VP2.get_vec());
    tstart = (double)clock()/(double)CLOCKS_PER_SEC;
    std::cout << "CDHR (cross) = "
              << volume_sequence_of_balls<CDHRWalk>(VP) << " , ";
    std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;

    VP.init(VP.dimension(), VP2.get_mat(), VP2.get_vec());
    tstart = (double)clock()/(double)CLOCKS_PER_SEC;
    std::cout << "Blrd (cross) = "
              << volume_sequence_of_balls<BilliardWalk>(VP) << " , ";
    std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;



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

    HP = generate_cube<Hpolytope>(10, false);

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


    VP.init(VP.dimension(), VP2.get_mat(), VP2.get_vec());
    auto VPCheBall = VP.ComputeInnerBall();
    {
        NT diameter;
        VP.comp_diam(diameter, VPCheBall.second);
        vars<NT, RNGType> var(rnum,n,walk_len,n_threads,err,e,0,0,0,
                              VPCheBall.second,diameter,rng,
                              urdist,urdist1,-1.0,false,false,false,
                              false,false,true,false,false,false);

        tstart = (double)clock()/(double)CLOCKS_PER_SEC;
        VP.init(VP.dimension(), VP2.get_mat(), VP2.get_vec());
        std::cout << "OLD Ball = " << volume(VP, var, VPCheBall)
                  << " , ";
        std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;
    }
    VP.init(VP.dimension(), VP2.get_mat(), VP2.get_vec());
    VPCheBall = VP.ComputeInnerBall();
    {
        NT diameter;
        VP.comp_diam(diameter, VPCheBall.second);
        vars<NT, RNGType> var(rnum,n,walk_len,n_threads,err,e,0,0,0,
                              VPCheBall.second,diameter,rng,
                              urdist,urdist1,-1.0,false,false,false,
                              false,false,false,true,false,false);

        tstart = (double)clock()/(double)CLOCKS_PER_SEC;
        VP.init(VP.dimension(), VP2.get_mat(), VP2.get_vec());
        std::cout << "OLD RDHR = " << volume(VP, var, VPCheBall)
                  << " , ";
        std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;
    }
    VP.init(VP.dimension(), VP2.get_mat(), VP2.get_vec());
    VPCheBall = VP.ComputeInnerBall();
    {
        NT diameter;
        VP.comp_diam(diameter, VPCheBall.second);
        vars<NT, RNGType> var(rnum,n,walk_len,n_threads,err,e,0,0,0,
                              VPCheBall.second,diameter,rng,
                              urdist,urdist1,-1.0,false,false,false,
                              false,false,false,false,true,false);

        tstart = (double)clock()/(double)CLOCKS_PER_SEC;
        VP.init(VP.dimension(), VP2.get_mat(), VP2.get_vec());
        std::cout << "OLD CDHR = " << volume(VP, var, VPCheBall)
                  << " , ";
        std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;
    }
    VP.init(VP.dimension(), VP2.get_mat(), VP2.get_vec());
    VPCheBall = VP.ComputeInnerBall();
    {
        NT diameter;
        VP.comp_diam(diameter, VPCheBall.second);
        vars<NT, RNGType> var(rnum,n,walk_len,n_threads,err,e,0,0,0,
                              VPCheBall.second,diameter,rng,
                              urdist,urdist1,-1.0,false,false,false,
                              false,false,false,false,false,true);

        tstart = (double)clock()/(double)CLOCKS_PER_SEC;
        VP.init(VP.dimension(), VP2.get_mat(), VP2.get_vec());
        std::cout << "OLD Blrd = " << volume(VP, var, VPCheBall)
                  << " , ";
        std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;
    }
    return 0;

    Hpolytope P = gen_skinny_cube<Hpolytope>(5);
    P.ComputeInnerBall();

    //Hpolytope P40 = generate_cube<Hpolytope>(40, false);
    //CheBall = P40.ComputeInnerBall();

    // Estimate the volume

    //std::cout << "Default (cube) = " << volume(HP) << std::endl;
    //std::cout << "Default (cube) = " << volume(HP, 0.5) << std::endl;
    //std::cout << "Default (cube) = " << volume(HP, 0.5, 2) << std::endl;

    //typedef BoostRandomNumberGenerator<boost::mt19937, NT> RNG;
    typedef BoostRandomNumberGenerator<boost::mt11213b, NT> RNG;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT, 3> RNGFixed;
    //RNG boost_rng(HP.dimension());
    //RNG3 boost_rng3(HP.dimension());

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
/*
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
              << volume_cooling_gaussians<GaussianBallWalk>(HP, e, walk_len) << " , ";
    std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;

    tstart = (double)clock()/(double)CLOCKS_PER_SEC;
    std::cout << "new GC RDHR (cube) = "
              << volume_cooling_gaussians<GaussianRDHRWalk>(HP, e, walk_len) << " , ";
    std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;

    tstart = (double)clock()/(double)CLOCKS_PER_SEC;
    std::cout << "new GC CDHR (cube) = "
              << volume_cooling_gaussians<GaussianCDHRWalk>(HP, e, walk_len) << " , ";
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
                  << volume_cooling_gaussians(HP, var1, var2, CheBall)
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
                  << volume_cooling_gaussians(HP, var1, var2, CheBall)
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
                  << volume_cooling_gaussians(HP, var1, var2, CheBall)
                  << " , ";
        std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;
    }
*/
    /*
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
*/
/*
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
*/
/*
    {
        NT diameter;
        HP.comp_diam(diameter, VPCheBall.second);
        vars<NT, RNGType> var(rnum,n,walk_len,n_threads,err,e,0,0,0,
                              VPCheBall.second,diameter,rng,
                              urdist,urdist1,-1.0,false,false,false,
                              false,false,false,false,true,false);

        tstart = (double)clock()/(double)CLOCKS_PER_SEC;
        //std::cout << "OLD CDHR = " << vol_cooling_balls(VP, var, var_ban, VPCheBall)
        //          << " , ";
        std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;
    }
*/
    /*
{
        VP.ComputeInnerBall();
    tstart = (double)clock()/(double)CLOCKS_PER_SEC;
    std::cout << "NEW Ball Vpoly = "
              << volume_cooling_balls<BallWalk>(VP, e, walk_len)
              << " , ";
    std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;
}
{
        VP.ComputeInnerBall();
    tstart = (double)clock()/(double)CLOCKS_PER_SEC;
    std::cout << "NEW CDHR Vpoly = "
              << volume_cooling_balls<CDHRWalk>(VP, e, walk_len)
              << " , ";
    std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;
}
    tstart = (double)clock()/(double)CLOCKS_PER_SEC;
    std::cout << "NEW RDHR Vpoly = "
              << volume_cooling_balls<RDHRWalk>(VP, e, walk_len)
              << " , ";
    std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;
    tstart = (double)clock()/(double)CLOCKS_PER_SEC;
    std::cout << "NEW Blrd Vpoly = "
              << volume_cooling_balls<BilliardWalk>(VP, e, walk_len)
              << " , ";
    std::cout << (double)clock()/(double)CLOCKS_PER_SEC - tstart << std::endl;
*/

    return 0;
}
