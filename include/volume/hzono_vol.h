// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

#ifndef HZONO_VOL_H
#define HZONO_VOL_H

#include "hpoly_utils/outer_zono.h"
#include "hpoly_utils/hpoly_annealing.h"
//#include "hpoly_utils/est_ratio1.h"
#include "esti_ratioGl.h"
#include "ZonoIntersectHPoly.h"


template <class Hpolytope, class Zonotope, class UParameters, class AParameters, class GParameters, class Point, typename NT>
NT vol_hzono (Zonotope &ZP, UParameters &var, AParameters &var_ban, GParameters &var_g,
              std::pair<Point,NT> &InnerB) {

    //typedef typename Zonotope::PolytopePoint Point;
    typedef typename Zonotope::VT VT;
    typedef typename Zonotope::MT MT;
    typedef typename UParameters::RNGType RNGType;

    NT lb = var_ban.lb, ub = var_ban.ub, prob = var_ban.p, B0_radius = var_ban.B0_radius, rmax = var_ban.rmax,
            radius = InnerB.second, round_value = 1.0, tele_prod = 1.0, e = var.error, ratio;
    int n = var.n, win_len = var_ban.win_len, N = var_ban.N, nu = var_ban.nu;

    var.verbose = true;
    bool verbose = var.verbose, round = var.round, window2 = var_ban.window2, cg_hpol = var_ban.hp_cg;
    //std::cout<<"error = "<<e<<std::endl;
    //int n= var.n;
    //NT ratio;


    //std::pair<Point,NT> InnerB = ZP.ComputeInnerBall();

    MT V = ZP.get_mat();
    MT G = V.transpose();
    int m = G.cols();
    std::list<Point> randPoints;

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
    //NT volh;
    MT T = ZP.get_T();
    MT Tt = T.transpose();
    MT A2 = AA * Tt;
    MT B = G * Tt;
    A3 = A2 * B.inverse();

    Hpolytope HP;
    HP.init(n,A3,b);

    VT Zs_max(2*m);
    get_hdelta(ZP, HP, Zs_max, lb, ub, ratio, randPoints, var);
    //if(verbose) std::cout<<"delta = "<<delta_in<<std::endl;
    Hpolytope HP2=HP;// = HP;
    //HP2.init(n,A3,b);
    std::pair<Point, NT> InnerBall = HP.ComputeInnerBall();

    //typedef Ball<Point> ball;
    typedef ZonoIntersectHPoly<Zonotope , Hpolytope > ZonoHP;
    //typedef std::list<Point> PointList;
    std::vector<Hpolytope > HPolySet;
    //std::vector<PointList> PointSets;
    std::vector<NT> ratios;
    NT p_value = 0.1;
    ZonoHP zb1, zb2;

    //var2.coordinate=false;
    //var2.walk_steps=1;
    get_sequence_of_zonopolys<ZonoHP>(ZP, HP2, HPolySet, Zs_max, ratios, N*nu, nu,
                                      lb, ub, var);

    //HnRsteps += steps;
    int mm=HPolySet.size()+2;
    int mm2=mm+1;
    prob = std::pow(prob, 1.0/NT(mm2));
    NT er0 = e/(2.0*std::sqrt(NT(mm2)));
    NT er1 = (e*std::sqrt(2.0*NT(mm2)-1))/(std::sqrt(2.0*NT(mm2)));
    NT Her = e/(2.0*std::sqrt(NT(mm2)));

    // var2.walk_steps=10 + n / 10;
    //InnerBall.first.print();
    //std::cout<<"radius = "<<InnerBall.second<<std::endl;
    //vars_g<NT, RNGType> var111(n, 1, N, 2*W, 1, Her, InnerBall.second, var.rng, C, frac,
                               //ratio2, -1.0, false, verbose, false, false, NN, birk,
                              // false, true);
    //var2.coordinate=true;
    NT vol;
    var_g.error = Her;
    var.error = Her;
    var.cdhr_walk = false;
    var.rdhr_walk = true;
    if( cg_hpol ) {
        vol = volume_gaussian_annealing(HP, var_g, var, InnerBall);
    } else {
        vol = volesti_ball_ann(HP, var, var_ban, InnerBall);
    }
    var.cdhr_walk = true;
    var.rdhr_walk = false;
    if(verbose) std::cout<<"\nvol of h-polytope = "<<vol<<"\n"<<std::endl;
    //double tstart3 = (double)clock()/(double)CLOCKS_PER_SEC;

    //NT ratio22;
    //if(len_subwin==0) len_subwin = 2;// + int(std::log2(NT(n)));
    //if(len_tuple==0) len_tuple = n*n+125;
    if (!window2) {
        //vol *= est_ratio_hzono_normal(HP2, ZP, er0, Win_len, prob, ratio, var22, Hsteps);
        vol *= esti_ratio_interval<RNGType, Point>(HP2, ZP, ratio, er0, win_len, N*nu, prob, var);
    } else {
        //vol *= est_ratio_hzono(ZP, HP2, er0, ratio, var22);
    }
    //MemLps += Hsteps;

    //double tstop3 = (double)clock()/(double)CLOCKS_PER_SEC;
    //if(verbose) std::cout << "[3] rejection time = " << tstop3 - tstart3 << std::endl;
    //if (verbose) std::cout<<"final ratio = "<<ratio22<<std::endl;
    //vol = vol * ratio22;

    randPoints.clear();
    Point q(n);

    Hpolytope b1, b2;
    //std::cout<<"coordinate: "<<var.coordinate<<std::endl;
    //NT er = e*0.8602325;
    if (HPolySet.size()==0) {
        if (verbose) std::cout << "no hpoly | ratio = " << ratios[0] << std::endl;
        if (ratios[0]!=1) {
            //vol = vol * est_ratio_zball_sym<Point>(ZP, sigma, G2, Q0, l, u, delta_in, ratios[0], e*0.8602325, var2);
            if(!window2) {
                //vol = vol / esti_ratio_interval<RNGType, Point>(ZP, HP2, ratios[0], er1, Win_len, prob, var, steps);
                vol = vol / esti_ratio_interval<RNGType, Point>(ZP, HP2, ratios[0], er1, win_len, N*nu, prob, var);
            } else {
                //vol = vol / est_ratio_zonoballs<Point>(ZP, HP2, ratios[0], er1, var2, steps);
            }
            //HnRsteps += steps;
        }
    } else {
        er1 = er1 / std::sqrt(NT(mm)-1.0);
        //er = er/std::sqrt(HPolySet.size()+1);
        if(verbose) std::cout<<"number of hpolys = "<<HPolySet.size()<<std::endl;
        b1 = HPolySet[0];
        //b1 = zb1.second();
        if(!window2) {
            vol = vol / esti_ratio_interval<RNGType, Point>(ZP, b1, ratios[0], er1, win_len, N*nu, prob, var);
        } else {
            //vol = vol / est_ratio_zonoballs<Point>(ZP, b1, ratios[0], er1, var2, steps);
        }
        //HnRsteps += steps;

        for (int i = 0; i < HPolySet.size()-1; ++i) {
            zb1 = ZonoHP(ZP,HPolySet[i]);
            //b1 = zb1.second();
            b2 = HPolySet[i+1];
            //b2 = zb2.second();
            if(!window2) {
                vol = vol / esti_ratio_interval<RNGType, Point>(zb1, b2, ratios[i], er1, win_len, N*nu, prob, var);
            } else {
                //vol = vol / est_ratio_zonoballs<Point>(zb1, b2, ratios[i], er1, var2, steps);
            }
            //HnRsteps += steps;
        }

        zb1 = ZonoHP(ZP,HPolySet[HPolySet.size()-1]);
        //vol = vol * est_ratio_zball_sym<Point>(zb1, sigma, G2, Q0, l, u, delta_in, ratios[ratios.size()-1], e, var2);
        if (!window2) {
            vol = vol / esti_ratio_interval<RNGType, Point>(zb1, HP2, ratios[ratios.size() - 1], er1, win_len, N*nu, prob, var);
        } else {
            //vol = vol / est_ratio_zonoballs<Point>(zb1, HP2, ratios[ratios.size() - 1], er1, var2, steps);
        }
        //HnRsteps += steps;
    }

    //Rcpp::NumericVector res(5);
   // res[0] = vol;
    //nHpoly = HPolySet.size()+1;
    //res[4] = (10+10/n)*MemLps + cg_steps;
    //res[4] = vol_pca;

    //std::cout<<"final volume = "<<vol<<std::endl;
    //std::cout<<"\nNumber of phases (h-polytopes) = "<<nHpoly<<"\nNumber of steps = "<<HnRsteps<<"\nVolume = "<<vol<<std::endl;
    return vol;

}

#endif
