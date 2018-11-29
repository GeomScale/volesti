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
#include <boost/math/special_functions/erf.hpp>
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
#include "ball_ann_vol.h"
#include "ZonoIntersectHPoly.h"



// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
Rcpp::NumericVector volume_zono (Rcpp::Reference P, double e=0.1, bool rounding = false, bool steps_only=false,
                               bool verbose=false, bool const_win=true, bool cg_hpol = false, bool PCA = false,
                               double lb_ratio=0.1, double ub_ratio=0.15, int len_subwin = 0, int len_tuple = 0) {

    typedef double NT;
    typedef Cartesian <NT> Kernel;
    typedef typename Kernel::Point Point;
    typedef boost::mt19937 RNGType;
    typedef HPolytope <Point> Hpolytope;
    //typedef VPolytope <Point, RNGType> Vpolytope;
    typedef Zonotope <Point> zonotope;
    //typedef IntersectionOfVpoly<Vpolytope> InterVP;
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
    //Hpolytope HP;

    int n;

    int type = P.field("t");
    MT V = Rcpp::as<MT>(P.field("G"));
    n = V.cols();
    VT vec = VT::Ones(V.rows());
    ZP.init(n, V, vec);
    coordinate = false;

    //NT HnRsteps, nballs, MemLps, vol;

    //Compute chebychev ball//
    std::pair<Point, NT> InnerBall;
    InnerBall = ZP.ComputeInnerBall();


    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    // the random engine with this seed
    RNGType rng(seed);
    boost::random::uniform_real_distribution<>(urdist);
    boost::random::uniform_real_distribution<> urdist1(-1,1);

    vars<NT, RNGType> var(1, n, 1, n_threads, 0.0, e, 0, 0.0, 0, InnerBall.second, rng,
                          urdist, urdist1, -1, verbose, rand_only, rounding, NN, birk, ball_walk, coordinate);
    NT HnRsteps=0.0, nballs, MemLps=0.0, vol, ratio;
    //------------------------------------------------------------------//
    typedef Ball<Point> ball;
    typedef BallIntersectPolytope<zonotope ,ball> ZonoBall;
    ball B0;
    get_first_ball<RNGType>(ZP, B0, ratio, InnerBall.second, var, MemLps);
    NT vol1= (std::pow(M_PI,n/2.0)*(std::pow(B0.radius(), n) ) ) / (tgamma(n/2.0+1));//*ratio;
    if(verbose)  std::cout<<"vol1 = "<<vol1*ratio<<std::endl;

    //------------------------------------------------------------------//
    MT G = V.transpose();
    int m = G.cols();
    MT AA; VT b;
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
    MT A3;
    NT volh;
    MT T = ZP.get_T();
    MT Tt = T.transpose();
    MT A2 = AA * Tt;
    MT B = G * Tt;
    A3 = A2 * B.inverse();
    Hpolytope HP;
    HP.init(n,A3,b);
    vars<NT, RNGType> var33(1, n, 10 + n / 10, 1, 0.0, e, 0, 0.0, 0, InnerBall.second, rng,
                            urdist, urdist1, -1.0, verbose, false, false, NN, birk, false,
                            true);
    NT steps, Hsteps, cg_steps, ratio2;
    VT Zs_max(2*m);
    std::list<Point> randPoints;
    get_hdelta<Point>(ZP, HP, Zs_max, ub_ratio, ratio2, randPoints, var33, Hsteps);
    std::pair<Point, NT> InnerBall2 = HP.ComputeInnerBall();
    var.error=e/5.0;
    var.walk_steps=1;
    NT HnRsteps2, nballs2, MemLps2;
    NT vol2 = volesti_ball_ann(HP, InnerBall2, lb_ratio, ub_ratio, var, HnRsteps2, nballs2, MemLps2, 0, 0, 0.75, false, false);//*ratio2;
    if(verbose) std::cout<<"vol2 = "<<vol2*ratio2<<std::endl;

    //------------------------------------------------------------//



    if(len_subwin==0) len_subwin = 2;// + int(std::log2(NT(n)));
    if(len_tuple==0) len_tuple = n*n+125;

    if(vol1*ratio>vol2*ratio2) {
        var.error=e;
        HnRsteps = 0.0;
        vol = volesti_ball_ann(ZP, InnerBall, lb_ratio, ub_ratio, var, HnRsteps, nballs, MemLps, len_subwin, len_tuple,
                               0.75,steps_only, const_win, B0.radius(), ratio);
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
    vol = vol2;

    typedef ZonoIntersectHPoly<zonotope , Hpolytope > ZonoHP;
    ZonoHP zb1, zb2;
    Hpolytope HP2 = HP;
    std::vector<Hpolytope> HPolySet;
    std::vector<NT> ratios;
    NT p_value=0.1;

    var.coordinate = false;
    get_sequence_of_zonoballs<ZonoHP>(ZP, HP2, HPolySet, Zs_max, ratios,
                                      p_value, var, steps);

    int mm=HPolySet.size()+2;
    NT prob = std::pow(0.75, 1.0/NT(mm));
    NT er0 = e/(2.0*std::sqrt(NT(mm)));
    NT er1 = (e*std::sqrt(2.0*NT(mm)-1))/(std::sqrt(2.0*NT(mm)));
    NT Her = e/(2.0*std::sqrt(NT(mm)));

    vars<NT, RNGType> var22(1, n, 10+n, 1, 0.0, e, 0, 0.0, 0, InnerBall.second, rng,
                            urdist, urdist1, -1.0, verbose, false, false, NN, birk, false,
                            false);
    NT ratio22;
    if(const_win){
        ratio22 = est_ratio_hzono_normal(ZP, HP2, er0, len_subwin, len_tuple, prob, ratio2, var22, Hsteps);
    } else {
        ratio22 = est_ratio_hzono(ZP, HP2, er0, ratio2, var22, Hsteps);
    }
    MemLps += Hsteps;

    //double tstop3 = (double)clock()/(double)CLOCKS_PER_SEC;
    //if(verbose) std::cout << "[3] rejection time = " << tstop3 - tstart3 << std::endl;
    if (verbose) std::cout<<"final ratio = "<<ratio22<<std::endl;
    vol = vol * ratio22;

    randPoints.clear();
    Point q(n);



    var.coordinate=false;
    var.walk_steps=1;
    //get_sequence_of_zonoballs<ball>(ZP, HP2, ZonoBallSet, PointSets, ratios,
    //p_value, var2, steps);



    Hpolytope b1, b2;
    //NT er = e*0.8602325;
    if (HPolySet.size()==0) {
        if (verbose) std::cout << "no hpoly | ratio = " << ratios[0] << std::endl;
        if (ratios[0]!=1) {
            //vol = vol * est_ratio_zball_sym<Point>(ZP, sigma, G2, Q0, l, u, delta_in, ratios[0], e*0.8602325, var2);
            if(const_win) {
                esti_ratio_interval<Point>(ZP, HP2, ratios[0], er1, len_subwin, len_tuple, prob, var, steps);
            } else {
                vol = vol * est_ratio_zonoballs<Point>(ZP, HP2, ratios[0], er1, var, steps);
            }
            HnRsteps += steps;
        }
    } else {
        //er = er/std::sqrt(HPolySet.size()+1);
        if(verbose) std::cout<<"number of hpolys = "<<HPolySet.size()<<std::endl;
        b1 = HPolySet[0];
        //b1 = zb1.second();
        if(const_win) {
            vol = vol / esti_ratio_interval<Point>(ZP, b1, ratios[0], er1, len_subwin, len_tuple, prob, var, steps);
        } else {
            vol = vol / est_ratio_zonoballs<Point>(ZP, b1, ratios[0], er1, var, steps);
        }
        HnRsteps += steps;

        for (int i = 0; i < HPolySet.size()-1; ++i) {
            zb1 = ZonoHP(ZP,HPolySet[i]);
            //b1 = zb1.second();
            b2 = HPolySet[i+1];
            //b2 = zb2.second();
            if(const_win) {
                vol = vol / esti_ratio_interval<Point>(zb1, b2, ratios[i], er1, len_subwin, len_tuple, prob, var, steps);
            } else {
                vol = vol / est_ratio_zonoballs<Point>(zb1, b2, ratios[i], er1, var, steps);
            }
            HnRsteps += steps;
        }

        zb1 = ZonoHP(ZP,HPolySet[HPolySet.size()-1]);
        //vol = vol * est_ratio_zball_sym<Point>(zb1, sigma, G2, Q0, l, u, delta_in, ratios[ratios.size()-1], e, var2);
        if (const_win) {
            vol = vol / esti_ratio_interval<Point>(zb1, HP2, ratios[ratios.size() - 1], er1, len_subwin, len_tuple, prob, var, steps);
        } else {
            vol = vol / est_ratio_zonoballs<Point>(zb1, HP2, ratios[ratios.size() - 1], er1, var, steps);
        }
        HnRsteps += steps;
    }

    Rcpp::NumericVector res(5);
    res[0] = vol;
    res[1] = NT(HPolySet.size()+1);
    res[2] = HnRsteps;
    res[3] = MemLps;
    res[4] = (10+10/n)*MemLps + cg_steps;

    //std::cout<<"final volume = "<<vol<<std::endl;
    return res;

}
