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
#include "gaussian_samplers.h"
#include "gaussian_annealing.h"
#include "rounding.h"
#include "zonovol_heads/cg_vol_test.h"


// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
Rcpp::NumericVector cg_volume(Rcpp::Reference P, double e = 0.1, bool steps_only = false, bool rounding = false, bool verbose = false) {

    typedef double NT;
    typedef Cartesian <NT> Kernel;
    typedef typename Kernel::Point Point;
    typedef boost::mt19937 RNGType;
    typedef HPolytope <Point> Hpolytope;
    typedef VPolytope <Point, RNGType> Vpolytope;
    typedef Zonotope <Point> zonotope;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
    typedef Eigen::Matrix <NT, Eigen::Dynamic, Eigen::Dynamic> MT;

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
    if (type == 1) {
        //std::cout << "H poly" << std::endl;
        MT A = Rcpp::as<MT>(P.field("A"));
        n = A.cols();
        VT vec = Rcpp::as<VT>(P.field("b"));
        HP.init(n, A, vec);

    } else if (type == 2) {
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
    std::pair <Point, NT> InnerBall;
    if (type == 1) {
        InnerBall = HP.ComputeInnerBall();
    } else if (type == 2) {
        InnerBall = VP.ComputeInnerBall();
    } else {
        InnerBall = ZP.ComputeInnerBall();
    }

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    // the random engine with this seed
    RNGType rng(seed);
    boost::random::uniform_real_distribution<>(urdist);
    boost::random::uniform_real_distribution<> urdist1(-1, 1);

    NT C = 2.0, frac=0.1, ratio2 = 1.0-1.0/(NT(n));
    int W = 4*n*n+500, N = 500 * ((int) C) + ((int) (n * n / 2));
    vars_g<NT, RNGType> var1(n, 1, N, W, 1, e, InnerBall.second, rng, C, frac, ratio2, -1.0, false, verbose,
                             rand_only, rounding, NN, birk, ball_walk, coordinate);

    vars <NT, RNGType> var2(1, n, 1, n_threads, 0.0, e, 0, 0.0, 0, InnerBall.second, rng,
                           urdist, urdist1, -1, verbose, rand_only, rounding, NN, birk, ball_walk, coordinate);

    NT HnRsteps, nGaussians, vol;
    if (type==1) {
        vol = params_gaussian_cooling(HP, var1, var2, InnerBall, HnRsteps, nGaussians, steps_only);
    } else if (type==2) { // if the input is a H-polytope
        vol = params_gaussian_cooling(VP, var1, var2, InnerBall, HnRsteps, nGaussians, steps_only);
    } else {  // if the input is a V-polytope
        vol = params_gaussian_cooling(ZP, var1, var2, InnerBall, HnRsteps, nGaussians, steps_only);
    }

    if (steps_only) {
        Rcpp::NumericVector res(1, vol);
        return res;
    }

    Rcpp::NumericVector res(3);
    res[0] = vol;
    res[1] = nGaussians;
    res[2] = HnRsteps;

    return res;

}
