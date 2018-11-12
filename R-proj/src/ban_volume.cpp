// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

#define VOLESTI_DEBUG
#include <Rcpp.h>
#include <RcppEigen.h>
// [[Rcpp::depends(BH)]]
//#define VOLESTI_DEBUG
#include <iterator>
//#include <fstream>
#include <vector>
#include <list>
//#include <algorithm>
#include <math.h>
#include <chrono>
#include "cartesian_geom/cartesian_kernel.h"
#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include "vars.h"
#include "polytopes.h"
//#include "ellipsoids.h"
#include "ballintersectconvex.h"
//#include "vpolyintersectvpoly.h"
#include "samplers.h"
#include "rounding.h"
#include "ball_ann_vol.h"


// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
Rcpp::NumericVector ban_volume(Rcpp::Reference P, double e = 0.1, bool steps_only = false, bool rounding = false, bool verbose = false,
                                double lb_ratio=0.1, double ub_ratio=0.15) {

    typedef double NT;
    typedef Cartesian <NT> Kernel;
    typedef typename Kernel::Point Point;
    typedef boost::mt19937 RNGType;
    typedef HPolytope <Point> Hpolytope;
    typedef VPolytope <Point, RNGType> Vpolytope;
    typedef Zonotope <Point> zonotope;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;

    bool rand_only = false,
            NN = false,
                    ball_walk = false,
                            coordinate = true,
            birk = false;
    //verbose = false;
    unsigned int n_threads = 1;
    zonotope ZP;
    Vpolytope VP;
    Hpolytope HP;
    int n;

    int type = P.field("t");
    if (type==1) {
        //std::cout<<"H poly"<<std::endl;
        MT A = Rcpp::as<MT>(P.field("A"));
        n = A.cols();
        VT vec = Rcpp::as<VT>(P.field("b"));
        HP.init(n, A, vec);

    } else if(type==2) {
        MT V = Rcpp::as<MT>(P.field("V"));
        n = V.cols();
        VT vec = VT::Ones(V.rows());
        VP.init(n, V, vec);
        coordinate = false;
    } else {
        MT V = Rcpp::as<MT>(P.field("G"));
        n = V.cols();
        VT vec = VT::Ones(V.rows());
        ZP.init(n, V, vec);
        coordinate = false;
    }

    //Compute chebychev ball//
    std::pair<Point, NT> InnerBall;
    if (type==1) {
        InnerBall = HP.ComputeInnerBall();
    } else if(type==2) {
        InnerBall = VP.ComputeInnerBall();
    }else{
        InnerBall = ZP.ComputeInnerBall();
    }

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    // the random engine with this seed
    RNGType rng(seed);
    boost::random::uniform_real_distribution<>(urdist);
    boost::random::uniform_real_distribution<> urdist1(-1,1);

    vars<NT, RNGType> var(1, n, 1, n_threads, 0.0, e, 0, 0.0, 0, InnerBall.second, rng,
                           urdist, urdist1, -1, verbose, rand_only, rounding, NN, birk, ball_walk, coordinate);
    NT HnRsteps, nballs, MemLps, vol;

    if(type==2) {
        Point p = InnerBall.first;
        std::list<Point> randPoints;
        rand_point_generator(VP, p, 20*n, 1, randPoints, var);
        VP.get_vol_centroid(InnerBall, randPoints);
    }

    if(type==1) {
        vol = volesti_ball_ann(HP, InnerBall, lb_ratio, ub_ratio, var, HnRsteps, nballs, MemLps, steps_only);
    } else if(type==2) {
        vol = volesti_ball_ann(VP, InnerBall, lb_ratio, ub_ratio, var, HnRsteps, nballs, MemLps, steps_only);
    } else {
        vol = volesti_ball_ann(ZP, InnerBall, lb_ratio, ub_ratio, var, HnRsteps, nballs, MemLps, steps_only);
    }

    if (steps_only) {
        Rcpp::NumericVector res(1, vol);
        return res;
    }

    Rcpp::NumericVector res(4);
    res[0] = vol;
    res[1] = nballs;
    res[2] = HnRsteps;
    res[3] = MemLps;

    return res;


}