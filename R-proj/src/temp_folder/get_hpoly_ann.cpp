

// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.
#define VOLESTI_DEBUG
#include <Rcpp.h>
#include <RcppEigen.h>
// [[Rcpp::depends(BH)]]
#define VOLESTI_DEBUG
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
#include "zonovol_heads/cg_hpoly.h"
//#include "gaussian_samplers.h"
//#include "gaussian_annealing.h"
#include "zonovol_heads/sampleTruncated.h"
#include "zonovol_heads/annealing_zono.h"
#include "zonovol_heads/outer_zono.h"
#include "zonovol_heads/cg_zonovol.h"
#include "zonovol_heads/ball_annealing.h"
#include "zonovol_heads/est_ratio1.h"



// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
Rcpp::NumericVector vol_zono (Rcpp::Reference P, double e, bool verbose=false, double delta_in=0.0,
                              double up_lim=0.15) {

    typedef double NT;
    typedef Cartesian <NT> Kernel;
    typedef typename Kernel::Point Point;
    typedef boost::mt19937 RNGType;
    typedef HPolytope <Point> Hpolytope;
    typedef VPolytope <Point, RNGType> Vpolytope;
    typedef Zonotope <Point> zonotope;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
    typedef Eigen::Matrix <NT, Eigen::Dynamic, Eigen::Dynamic> MT;
    //unsigned int n_threads=1,i,j;

    bool rand_only = false,
            NN = false,
            birk = false;
    //verbose = false;
    unsigned int n_threads = 1;
    zonotope ZP;
    Rcpp::NumericMatrix A = P.field("G");

    unsigned int m = A.nrow();
    unsigned int n = A.ncol();
    //std::cout<<"hello"<<std::endl;
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    // the random engine with this seed
    RNGType rng(seed);
    boost::random::uniform_real_distribution<>(urdist);
    boost::random::uniform_real_distribution<> urdist1(-1, 1);

    MT V = Rcpp::as<MT>(A);
    VT vec = VT::Ones(m);
    ZP.init(n, V, vec);
    MT G = V.transpose();
    std::list <Point> randPoints;
    NT delta, ratio;

    std::pair <Point, NT> InnerB = ZP.ComputeInnerBall();
    NT C = 2.0, frac = 0.1, ratio2 = 1.0 - 1.0 / (NT(n));
    int W = 4 * n * n + 500, N = 500 * ((int) C) + ((int) (n * n / 2));

    vars <NT, RNGType> var2(1, n, 10 + n / 10, 1, 0.0, e, 0, 0.0, 0, InnerB.second, rng,
                            urdist, urdist1, -1.0, verbose, false, false, NN, birk, false,
                            false);

    MT AA;
    VT b;
    AA.resize(2 * m, m);
    b.resize(2 * m);
    for (unsigned int i = 0; i < m; ++i) {
        b(i) = 1.0;
        for (unsigned int j = 0; j < m; ++j) {
            if (i == j) {
                AA(i, j) = 1.0;
            } else {
                AA(i, j) = 0.0;
            }
        }
    }
    for (unsigned int i = 0; i < m; ++i) {
        b(i + m) = 1.0;
        for (unsigned int j = 0; j < m; ++j) {
            if (i == j) {
                AA(i + m, j) = -1.0;
            } else {
                AA(i + m, j) = 0.0;
            }
        }
    }
    MT T = ZP.get_T();
    MT Tt = T.transpose();
    MT A2 = AA * Tt;
    MT B = G * Tt;
    MT A3 = A2 * B.inverse();
    Hpolytope HP;
    HP.init(n, A3, b);

    var2.coordinate=true;
    NT steps, Hsteps, HnRsteps = 0.0, MemLps = 0.0, cg_steps;
    get_hdelta<Point>(ZP, HP, delta_in, up_lim, ratio, randPoints, var2, Hsteps);
    var2.coordinate = true;

    typedef Ball <Point> ball;
    typedef BallIntersectPolytope <zonotope, ball> ZonoBall;
    typedef std::list <Point> PointList;
    std::vector <ZonoBall> ZonoBallSet;
    std::vector <PointList> PointSets;
    std::vector <NT> ratios;
    NT p_value = 0.1;
  //  ZonoBall zb1, zb2;

    var2.coordinate = false;
    var2.walk_steps = 1;
    get_sequence_of_zonoballs<ball>(ZP, HP2, ZonoBallSet, PointSets, ratios,
                                    p_value, var2, steps);

    return ZonoBallSet.size()+1;
}