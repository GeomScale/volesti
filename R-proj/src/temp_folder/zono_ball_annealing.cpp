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
Rcpp::NumericVector ZonoBallAnn(Rcpp::Reference P, double e, bool verbose = false,
                                double lb_ratio=0.1, double ub_ratio=0.2) {

    typedef double NT;
    typedef Cartesian <NT> Kernel;
    typedef typename Kernel::Point Point;
    typedef boost::mt19937 RNGType;
    typedef HPolytope <Point> Hpolytope;
    typedef VPolytope <Point, RNGType> Vpolytope;
    typedef Zonotope <Point> zonotope;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    //unsigned int n_threads=1,i,j;

    bool rand_only = false,
            NN = false,
            birk = false;
    //verbose = false;
    unsigned int n_threads = 1;
    zonotope ZP;
    Rcpp::NumericMatrix A = P.field("G");

    unsigned int m=A.nrow();
    unsigned int n=A.ncol();

    MT V = Rcpp::as<MT>(A);
    VT vec = VT::Ones(m);
    ZP.init(n, V, vec);

    typedef Ball<Point> ball;
    typedef BallIntersectPolytope<zonotope , ball > ZonoBall;
    typedef std::list<Point> PointList;
    std::vector<ball> BallSet;
    std::vector<PointList> PointSets;
    std::vector<NT> ratios;
    NT p_value = 0.1;

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    // the random engine with this seed
    RNGType rng(seed);
    boost::random::uniform_real_distribution<>(urdist);
    boost::random::uniform_real_distribution<> urdist1(-1,1);

    vars<NT, RNGType> var2(1, n, 10 + n / 10, 1, 0.0, 0.1, 0, 0.0, 0, 0.0, rng,
                           urdist, urdist1, -1.0, verbose, false, false, NN, birk, false,
                           false);
    std::pair<Point,NT> InnerB = ZP.ComputeInnerBall();
    ball B0;
    NT ratio0, steps, HnRSteps = 0.0, MemLps=0.0, ballsteps;

    if(verbose) std::cout<<"Computing ball annealing..."<<std::endl;
    get_sequence_of_zonoballs<ZonoBall, RNGType>(ZP, BallSet, B0, ratio0, PointSets,
                                                 ratios, lb_ratio, ub_ratio, InnerB.second, var2,
                                                 ballsteps, steps);
    HnRSteps = steps;
    MemLps = ballsteps;

    typename std::vector<NT>::iterator rit = ratios.begin();
    for ( ; rit!=ratios.end(); ++rit) {
        if(verbose) std::cout << *rit << " ";
    }
    if(verbose) std::cout<<"\n";
    if(verbose) std::cout<<"size of ballSet = "<<BallSet.size()<<std::endl;
    if(verbose) std::cout<<"B0 radius = "<<B0.radius()<<std::endl;

    NT vol = (std::pow(M_PI,n/2.0)*(std::pow(B0.radius(), n) ) ) / (tgamma(n/2.0+1));
    vol = vol * esti_ratio2<RNGType>(B0, ZP, e/3.0, ratio0, steps);
    MemLps += steps;

    ball Biter;
    ZonoBall zb1, zb2;
    var2.coordinate = true;
    var2.walk_steps=1;
    NT tele_prod=1.0;
    NT er = e*0.942809;
    if (BallSet.size()>0) {
        er = er / std::sqrt(NT(BallSet.size()+1));
        tele_prod = tele_prod * esti_ratio(ZP, BallSet[0], ratios[0], er, var2, steps);
        HnRSteps += steps;
        for (int i = 0; i < BallSet.size() - 1; ++i) {
            zb1 = ZonoBall(ZP, BallSet[i]);
            //zb2 = ZonoBall(ZP, BallSet[i]);
            tele_prod = tele_prod * esti_ratio(zb1, BallSet[i+1], ratios[i+1], er, var2, steps);
            HnRSteps += steps;
        }


        zb1 = ZonoBall(ZP, BallSet[BallSet.size() - 1]);
        tele_prod = tele_prod * esti_ratio(zb1, B0, ratios[ratios.size() - 1], er, var2, steps);
        HnRSteps += steps;
        vol = vol / tele_prod;
    } else {
        if (ratios[0]!=1) {
            vol = vol / esti_ratio(ZP, B0, ratios[0], er, var2, steps);
            HnRSteps += steps;
        }
    }

    Rcpp::NumericVector res(4);
    res[0] = vol;
    res[1] = NT(BallSet.size()+1);
    res[2] = HnRSteps;
    res[3] = MemLps;

    return res;
}