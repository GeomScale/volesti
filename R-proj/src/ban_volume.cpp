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
#include <boost/math/special_functions/erf.hpp>
#include "vars.h"
#include "polytopes.h"

//#include "ellipsoids.h"
#include "ballintersectconvex.h"
//#include "vpolyintersectvpoly.h"
#include "samplers.h"
#include "vpolyintersectvpoly.h"
#include "rounding.h"
#include "ball_ann_vol.h"


// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
Rcpp::NumericVector ban_volume(Rcpp::Reference P, double e = 0.1, bool steps_only = false, bool const_win = true, bool rounding = false, bool verbose = false,
                                double lb_ratio=0.1, double ub_ratio=0.15, int len_subwin = 0, int len_tuple = 0) {

    typedef double NT;
    typedef Cartesian <NT> Kernel;
    typedef typename Kernel::Point Point;
    typedef boost::mt19937 RNGType;
    typedef HPolytope <Point> Hpolytope;
    typedef VPolytope <Point, RNGType> Vpolytope;
    typedef Zonotope <Point> zonotope;
    typedef IntersectionOfVpoly<Vpolytope> InterVP;
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
    Vpolytope VP2;
    InterVP VPcVP;
    int n;

    int type = P.field("t");
    if (type==1) {
        //std::cout<<"H poly"<<std::endl;
        MT A = Rcpp::as<MT>(P.field("A"));
        n = A.cols();
        VT vec = Rcpp::as<VT>(P.field("b"));
        HP.init(n, A, vec);
        //coordinate = false;
    } else if(type==2) {
        MT V = Rcpp::as<MT>(P.field("V"));
        n = V.cols();
        VT vec = VT::Ones(V.rows());
        VP.init(n, V, vec);
        coordinate = false;
    } else if(type==3) {
        MT V = Rcpp::as<MT>(P.field("G"));
        n = V.cols();
        VT vec = VT::Ones(V.rows());
        ZP.init(n, V, vec);
        coordinate = false;
    } else {
        MT V1 = Rcpp::as<MT>(P.field("V"));
        MT V2 = Rcpp::as<MT>(P.field("V2"));
        n = V1.cols();
        VT vec1 = VT::Ones(V1.rows());
        VT vec2 = VT::Ones(V2.rows());
        VP.init(n, V1, vec1);
        VP2.init(n, V2, vec2);
        VPcVP.init(VP, VP2);
        coordinate = false;
    }

    //Compute chebychev ball//
    std::pair<Point, NT> InnerBall;
    if (type==1) {
        InnerBall = HP.ComputeInnerBall();
    } else if(type==2) {
        InnerBall = VP.ComputeInnerBall();
    }else if(type==3){
        InnerBall = ZP.ComputeInnerBall();
    } else {
        bool empty=false;
        InnerBall = VPcVP.getInnerPoint_rad(empty);
        if(empty){
            Rcpp::NumericVector res;
            return res;
        }
    }

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    // the random engine with this seed
    RNGType rng(seed);
    boost::random::uniform_real_distribution<>(urdist);
    boost::random::uniform_real_distribution<> urdist1(-1,1);

    vars<NT, RNGType> var(1, n, 1, n_threads, 0.0, e, 0, 0.0, 0, InnerBall.second, rng,
                           urdist, urdist1, -1, verbose, rand_only, rounding, NN, birk, ball_walk, coordinate);
    NT HnRsteps, nballs, MemLps, vol;
    NT round_val = 1.0;

    if(type==2) {
        if (rounding) {
            std::pair<NT,NT> res_round = rounding_min_ellipsoid(VP,InnerBall,var);
            round_val=res_round.first;
            InnerBall = VP.ComputeInnerBall();
        }
        Point p = InnerBall.first;
        std::list<Point> randPoints;
        rand_point_generator(VP, p, 20*n, 1, randPoints, var);
        VP.get_vol_centroid(InnerBall, randPoints);
    } else if(type==4){
        //Point p = InnerBall.first;
        //std::list<Point> randPoints;
        //rand_point_generator(VP, p, 20*n, 1, randPoints, var);
        //VPcVP.get_vol_centroid(InnerBall, randPoints);
    }

    if(len_subwin==0) len_subwin = 30;// + int(std::log2(NT(n)));
    if(len_tuple==0) len_tuple = 150+n;
    if(type==1) {
        vol = volesti_ball_ann(HP, InnerBall, lb_ratio, ub_ratio, var, HnRsteps, nballs, MemLps, len_subwin, len_tuple, steps_only, const_win);
    } else if(type==2) {
        vol = volesti_ball_ann(VP, InnerBall, lb_ratio, ub_ratio, var, HnRsteps, nballs, MemLps, len_subwin, len_tuple, steps_only, const_win);
    } else if(type==3){
        vol = volesti_ball_ann(ZP, InnerBall, lb_ratio, ub_ratio, var, HnRsteps, nballs, MemLps, len_subwin, len_tuple, steps_only, const_win);
    } else {
        vol = volesti_ball_ann(VPcVP, InnerBall, lb_ratio, ub_ratio, var, HnRsteps, nballs, MemLps, len_subwin, len_tuple, steps_only, const_win);
    }

    if (steps_only) {
        Rcpp::NumericVector res(1, vol);
        return res;
    }

    Rcpp::NumericVector res(4);
    res[0] = vol*round_val;
    res[1] = nballs;
    res[2] = HnRsteps;
    res[3] = MemLps;

    return res;


}