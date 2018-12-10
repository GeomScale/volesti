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
Rcpp::NumericVector vol_hzono (Rcpp::Reference P, double e=0.1, bool steps_only=false, bool verbose=false, bool const_win=true,
                               bool cg_hpol = false, bool PCA = false, bool pca_ratio = false, double delta_in=0.0, double up_lim=0.15,
                               int len_subwin = 0, int len_tuple = 0) {

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
    //std::cout<<"hello"<<std::endl;
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    // the random engine with this seed
    RNGType rng(seed);
    boost::random::uniform_real_distribution<>(urdist);
    boost::random::uniform_real_distribution<> urdist1(-1,1);

    MT V = Rcpp::as<MT>(A);
    VT vec = VT::Ones(m);
    ZP.init(n, V, vec);
    MT G = V.transpose();
    std::list<Point> randPoints;
    NT delta, ratio;

    std::pair<Point,NT> InnerB = ZP.ComputeInnerBall();
    NT C = 2.0, frac=0.1, ratio2 = 1.0-1.0/(NT(n));
    int W = 4*n*n+500, N = 500 * ((int) C) + ((int) (n * n / 2));

    vars<NT, RNGType> var2(1, n, 10 + n / 10, 1, 0.0, e, 0, 0.0, 0, InnerB.second, rng,
                           urdist, urdist1, -1.0, verbose, false, false, NN, birk, false,
                           false);
    vars_g<NT, RNGType> var1(n, 1, N, W, 1, e/10.0, InnerB.second, rng, C, frac,
                             ratio2, -1.0, false, verbose, false, false, NN, birk,
                             false, false);

    if(PCA) m = n;

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
    if(!PCA) {
        MT T = ZP.get_T();
        MT Tt = T.transpose();
        MT A2 = AA * Tt;
        MT B = G * Tt;
        A3 = A2 * B.inverse();
    } else {
        MT X(n, 2*G.cols());
        X << G, -G;
        MT X2 = X.transpose();
        MT C = X2.transpose()*X2;
        Eigen::JacobiSVD<MT> svd(C, Eigen::ComputeFullU | Eigen::ComputeFullV);
        volh = std::pow(2.0,NT(n))*svd.matrixU().determinant();
        A3 = AA*svd.matrixU().transpose();
    }

    MT A33;
    NT vol_pca;
    if(pca_ratio) {

        int mn=n;
        MT AA2; VT b2;
        AA2.resize(2 * mn, mn);
        b2.resize(2 * mn);
        for (unsigned int i = 0; i < mn; ++i) {
            b2(i) = 1.0;
            for (unsigned int j = 0; j < mn; ++j) {
                if (i == j) {
                    AA2(i, j) = 1.0;
                } else {
                    AA2(i, j) = 0.0;
                }
            }
        }
        for (unsigned int i = 0; i < mn; ++i) {
            b2(i + mn) = 1.0;
            for (unsigned int j = 0; j < mn; ++j) {
                if (i == j) {
                    AA2(i + mn, j) = -1.0;
                } else {
                    AA2(i + mn, j) = 0.0;
                }
            }
        }

        MT X(n, 2*G.cols());
        X << G, -G;
        MT X2 = X.transpose();
        MT C = X2.transpose()*X2;
        Eigen::JacobiSVD<MT> svd(C, Eigen::ComputeFullU | Eigen::ComputeFullV);
        //volh = std::pow(2.0,NT(n))*svd.matrixU().determinant();
        A33 = AA2*svd.matrixU().transpose();
        Hpolytope HP33;
        HP33.init(n,A33,b);
        std::pair<Point, NT> InnerBall33 = HP33.ComputeInnerBall();
        std::cout<<"radius = "<<InnerBall33.second<<std::endl;
        vars<NT, RNGType> var33(1, n, 1, 1, 0.0, e, 0, 0.0, 0, InnerBall33.second, rng,
                                urdist, urdist1, -1.0, verbose, false, false, NN, birk, false,
                                false);
        NT HnRsteps2, nballs, MemLps2;
        NT lb_ratio=0.1, up_ratio=0.15;
        VT Zs33_max(2*n);
        get_hdelta<Point>(ZP, HP33, Zs33_max, up_lim, ratio, randPoints, var33, HnRsteps2);
        HP33.set_vec(Zs33_max);
        vol_pca = volesti_ball_ann(HP33, InnerBall33, lb_ratio, up_ratio, var33, HnRsteps2, nballs, MemLps2, 0, 0, 0.75, false, false);
    }


    Hpolytope HP;
    HP.init(n,A3,b);
    vars<NT, RNGType> var33(1, n, 10 + n / 10, 1, 0.0, e, 0, 0.0, 0, InnerB.second, rng,
                            urdist, urdist1, -1.0, verbose, false, false, NN, birk, false,
                            true);
    NT steps, Hsteps, HnRsteps = 0.0, MemLps = 0.0, cg_steps;
    VT Zs_max(2*m);
    get_hdelta<Point>(ZP, HP, Zs_max, up_lim, ratio, randPoints, var33, Hsteps);
    //if(verbose) std::cout<<"delta = "<<delta_in<<std::endl;
    MemLps += Hsteps;
    Hpolytope HP2=HP;// = HP;
    //HP2.init(n,A3,b);
    std::pair<Point, NT> InnerBall = HP.ComputeInnerBall();

    typedef Ball<Point> ball;
    typedef ZonoIntersectHPoly<zonotope , Hpolytope > ZonoHP;
    typedef std::list<Point> PointList;
    std::vector<Hpolytope > HPolySet;
    std::vector<PointList> PointSets;
    std::vector<NT> ratios;
    NT p_value = 0.1;
    ZonoHP zb1, zb2;

    var2.coordinate=false;
    var2.walk_steps=1;
    get_sequence_of_zonoballs<ZonoHP>(ZP, HP2, HPolySet, Zs_max, ratios,
                                    p_value, var2, steps);
    if(steps_only) {
        Rcpp::NumericVector res(1);
        res[0] = NT(HPolySet.size() + 1);
        return res;
    }
    HnRsteps += steps;
    int mm=HPolySet.size()+2;
    int mm2=mm+1;
    NT prob = std::pow(0.75, 1.0/NT(mm2));
    NT er0 = e/(2.0*std::sqrt(NT(mm2)));
    NT er1 = (e*std::sqrt(2.0*NT(mm2)-1))/(std::sqrt(2.0*NT(mm2)));
    NT Her = e/(2.0*std::sqrt(NT(mm2)));

    var2.walk_steps=10 + n / 10;
    //InnerBall.first.print();
    //std::cout<<"radius = "<<InnerBall.second<<std::endl;
    vars_g<NT, RNGType> var111(n, 1, N, 2*W, 1, Her, InnerBall.second, rng, C, frac,
                               ratio2, -1.0, false, verbose, false, false, NN, birk,
                               false, true);
    var2.coordinate=true;
    NT vol;
    if( cg_hpol ) {
        vol = volume_gaussian_annealing(HP, var111, var2, InnerBall);
    } else {
        NT HnRsteps2, nballs, MemLps2;
        NT lb_ratio=0.1, up_ratio=0.15;
        var2.error=Her;
        var2.walk_steps=1;
        vol = volesti_ball_ann(HP, InnerBall, lb_ratio, up_ratio, var2, HnRsteps2, nballs, MemLps2, 0, 0, 0.75, false, false);
    }
    cg_steps = 0.0;
    if(verbose) std::cout<<"\n\nvol of h-polytope = "<<vol<<"volhPCA = "<<volh<<"\n\n"<<std::endl;
    double tstart3 = (double)clock()/(double)CLOCKS_PER_SEC;

    vars<NT, RNGType> var22(1, n, 10+n, 1, 0.0, e, 0, 0.0, 0, InnerB.second, rng,
                            urdist, urdist1, -1.0, verbose, false, false, NN, birk, false,
                            false);
    NT ratio22;
    if(len_subwin==0) len_subwin = 2;// + int(std::log2(NT(n)));
    if(len_tuple==0) len_tuple = n*n+125;
    if(const_win){
        ratio22 = est_ratio_hzono_normal(ZP, HP2, er0, len_subwin, len_tuple, prob, ratio, var22, Hsteps);
    } else {
        ratio22 = est_ratio_hzono(ZP, HP2, er0, ratio, var22, Hsteps);
    }
    MemLps += Hsteps;

    double tstop3 = (double)clock()/(double)CLOCKS_PER_SEC;
    if(verbose) std::cout << "[3] rejection time = " << tstop3 - tstart3 << std::endl;
    if (verbose) std::cout<<"final ratio = "<<ratio22<<std::endl;
    vol = vol * ratio22;

    randPoints.clear();
    Point q(n);



    var2.coordinate=false;
    var2.walk_steps=1;
    //get_sequence_of_zonoballs<ball>(ZP, HP2, ZonoBallSet, PointSets, ratios,
                                    //p_value, var2, steps);



    Hpolytope b1, b2;
    //NT er = e*0.8602325;
    if (HPolySet.size()==0) {
        if (verbose) std::cout << "no ball | ratio = " << ratios[0] << std::endl;
        if (ratios[0]!=1) {
            //vol = vol * est_ratio_zball_sym<Point>(ZP, sigma, G2, Q0, l, u, delta_in, ratios[0], e*0.8602325, var2);
            if(const_win) {
                vol = vol / esti_ratio_interval<Point>(ZP, HP2, ratios[0], er1, len_subwin, len_tuple, prob, var2, steps);
            } else {
                vol = vol / est_ratio_zonoballs<Point>(ZP, HP2, ratios[0], er1, var2, steps);
            }
            HnRsteps += steps;
        }
    } else {
        er1 = er1 / std::sqrt(NT(mm)-1.0);
        //er = er/std::sqrt(HPolySet.size()+1);
        if(verbose) std::cout<<"number of balls = "<<HPolySet.size()<<std::endl;
        b1 = HPolySet[0];
        //b1 = zb1.second();
        if(const_win) {
            vol = vol / esti_ratio_interval<Point>(ZP, b1, ratios[0], er1, len_subwin, len_tuple, prob, var2, steps);
        } else {
            vol = vol / est_ratio_zonoballs<Point>(ZP, b1, ratios[0], er1, var2, steps);
        }
        HnRsteps += steps;

        for (int i = 0; i < HPolySet.size()-1; ++i) {
            zb1 = ZonoHP(ZP,HPolySet[i]);
            //b1 = zb1.second();
            b2 = HPolySet[i+1];
            //b2 = zb2.second();
            if(const_win) {
                vol = vol / esti_ratio_interval<Point>(zb1, b2, ratios[i], er1, len_subwin, len_tuple, prob, var2, steps);
            } else {
                vol = vol / est_ratio_zonoballs<Point>(zb1, b2, ratios[i], er1, var2, steps);
            }
            HnRsteps += steps;
        }

        zb1 = ZonoHP(ZP,HPolySet[HPolySet.size()-1]);
        //vol = vol * est_ratio_zball_sym<Point>(zb1, sigma, G2, Q0, l, u, delta_in, ratios[ratios.size()-1], e, var2);
        if (const_win) {
            vol = vol / esti_ratio_interval<Point>(zb1, HP2, ratios[ratios.size() - 1], er1, len_subwin, len_tuple, prob, var2, steps);
        } else {
            vol = vol / est_ratio_zonoballs<Point>(zb1, HP2, ratios[ratios.size() - 1], er1, var2, steps);
        }
        HnRsteps += steps;
    }

    Rcpp::NumericVector res(5);
    res[0] = vol;
    res[1] = NT(HPolySet.size()+1);
    res[2] = HnRsteps;
    res[3] = MemLps;
    //res[4] = (10+10/n)*MemLps + cg_steps;
    res[4] = vol_pca;

    //std::cout<<"final volume = "<<vol<<std::endl;
    return res;
}