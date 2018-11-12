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
#include "zonovol_heads/ball_annealing2.h"
#include "zonovol_heads/ball_ratios.h"

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
int get_ball_ann(Rcpp::Reference P, bool verbose = false,
                                double lb_ratio=0.1, double ub_ratio=0.15) {

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

    MT V = Rcpp::as<MT>(A);
    VT vec = VT::Ones(m);
    ZP.init(n, V, vec);

    typedef Ball <Point> ball;
    typedef BallIntersectPolytope <zonotope, ball> ZonoBall;
    typedef std::list <Point> PointList;
    std::vector <ball> BallSet;
    std::vector <PointList> PointSets;
    std::vector <NT> ratios;
    //NT p_value = 0.1;

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    // the random engine with this seed
    RNGType rng(seed);
    boost::random::uniform_real_distribution<>(urdist);
    boost::random::uniform_real_distribution<> urdist1(-1, 1);

    vars <NT, RNGType> var2(1, n, 10 + n / 10, 1, 0.0, 0.1, 0, 0.0, 0, 0.0, rng,
                            urdist, urdist1, -1.0, verbose, false, false, NN, birk, false,
                            false);
    std::pair <Point, NT> InnerB = ZP.ComputeInnerBall();
    ball B0;
    NT ratio0, steps, HnRSteps = 0.0, MemLps = 0.0, ballsteps;

    if (verbose) std::cout << "Computing ball annealing..." << std::endl;
    get_sequence_of_zonoballs<ZonoBall, RNGType>(ZP, BallSet, B0, ratio0, PointSets,
                                                 ratios, lb_ratio, ub_ratio, InnerB.second, var2,
                                                 ballsteps, steps);

    return BallSet.size()+1;

}