// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

#ifndef COOLING_HPOLY_H
#define COOLING_HPOLY_H

#include "ball_annealing.h"
#include "hpoly_annealing.h"
#include "ratio_estimation.h"
#include "zonoIntersecthpoly.h"


template <class Hpolytope, class Zonotope, class UParameters, class AParameters, class GParameters, class Point, typename NT>
NT vol_cooling_hpoly (Zonotope &ZP, UParameters &var, AParameters &var_ban, GParameters &var_g,
              std::pair<Point,NT> &InnerB) {

    typedef typename Zonotope::VT VT;
    typedef typename Zonotope::MT MT;
    typedef typename UParameters::RNGType RNGType;

    NT lb = var_ban.lb, ub = var_ban.ub, prob = var_ban.p, rmax = var_ban.rmax, diam0 = var.diameter,
            radius = InnerB.second, round_value = 1.0, e = var.error, ratio, vol, alpha = var_ban.alpha;
    int n = var.n, win_len = var_ban.win_len, N = var_ban.N, nu = var_ban.nu;

  
    bool verbose = var.verbose, round = var.round, window2 = var_ban.window2;

    MT V = ZP.get_mat();
    MT G = V.transpose();
    int m = G.cols();
    std::list<Point> randPoints;

    MT XX(m, 2*m);
    XX << MT::Identity(m,m), -MT::Identity(m,m);
    MT AA = XX.transpose(); VT b = VT::Ones(2*m);
    MT T = ZP.get_T();
    MT Tt = T.transpose();
    MT A2 = AA * Tt, B = G * Tt;
    MT A3 = A2 * B.inverse();

    NT row_norm;
    for (int i = 0; i < A3.rows(); ++i) {
        row_norm = A3.row(i).norm();
        A3.row(i) = A3.row(i) / row_norm;
        b(i) = b(i) / row_norm;
    }
    MT A = A3*G;

    Hpolytope HP;
    HP.init(n,A3,b);

    VT Zs_max(2*m);
    UParameters var3 = var;
    var3.cdhr_walk = true;
    var3.ball_walk = var3.rdhr_walk = false;
    get_hdelta(ZP, HP, Zs_max, lb, ub, ratio, var3);
    //Hpolytope HP2 = HP;
    HP.normalize();

    std::pair<Point, NT> InnerBall = HP.ComputeInnerBall();

    typedef ZonoIntersectHPoly<Zonotope , Hpolytope > ZonoHP;
    std::vector<Hpolytope > HPolySet;
    std::vector<NT> ratios;
    NT p_value = 0.1;
    ZonoHP zb1, zb2;
    std::vector<NT> diams_inter;

    get_sequence_of_zonopolys<ZonoHP>(ZP, HP, HPolySet, Zs_max, ratios, N*nu, nu, lb, ub, alpha, var, var3, diams_inter);
    //var.diameter = diam0;

    int mm=HPolySet.size()+2;
    int mm2=mm+1;
    prob = std::pow(prob, 1.0/NT(mm2));
    NT er0 = e/(2.0*std::sqrt(NT(mm2))), er1 = (e*std::sqrt(2.0*NT(mm2)-1))/(std::sqrt(2.0*NT(mm2))),
        Her = e/(2.0*std::sqrt(NT(mm2)));

    var_g.error = Her;
    vol = volume_gaussian_annealing(HP, var_g, var, InnerBall);

    if (!window2) {
        UParameters var2 = var;
        var2.cdhr_walk = true;
        var2.ball_walk = var2.rdhr_walk = false;
        var2.walk_steps = 10+2*n;
        vol *= esti_ratio_interval<RNGType, Point>(HP, ZP, ratio, er0, win_len, N*nu, prob, var2);
    } else {
        vol *= esti_ratio<RNGType, Point>(HP, ZP, ratio, er0, var_g.W, N*nu, var);
    }

    Hpolytope b1, b2;
    if (HPolySet.size()==0) {
        if (ratios[0]!=1) {
            if(!window2) {
                vol = vol / esti_ratio_interval<RNGType, Point>(ZP, HP, ratios[0], er1, win_len, N*nu, prob, var);
            } else {
                vol = vol / esti_ratio<RNGType, Point>(ZP, HP, ratios[0], er1, var_g.W, N*nu, var);
            }
        }
    } else {
        er1 = er1 / std::sqrt(NT(mm)-1.0);
        b1 = HPolySet[0];
        if(!window2) {
            vol = vol / esti_ratio_interval<RNGType, Point>(ZP, b1, ratios[0], er1, win_len, N*nu, prob, var);
        } else {
            vol = vol / esti_ratio<RNGType, Point>(ZP, b1, ratios[0], er1, var_g.W, N*nu, var);
        }

        for (int i = 0; i < HPolySet.size()-1; ++i) {
            zb1 = ZonoHP(ZP,HPolySet[i]);
            b2 = HPolySet[i+1];
            //var.diameter = diams_inter[i];
            if(!window2) {
                vol = vol / esti_ratio_interval<RNGType, Point>(zb1, b2, ratios[i], er1, win_len, N*nu, prob, var);
            } else {
                vol = vol / esti_ratio<RNGType, Point>(zb1, b2, ratios[i], er1, var_g.W, N*nu, var);
            }
        }

        zb1 = ZonoHP(ZP,HPolySet[HPolySet.size()-1]);
        //var.diameter = diams_inter[diams_inter.size()-1];
        if (!window2) {
            vol = vol / esti_ratio_interval<RNGType, Point>(zb1, HP, ratios[ratios.size() - 1], er1, win_len, N*nu, prob, var);
        } else {
            vol = vol / esti_ratio<RNGType, Point>(zb1, HP, ratios[ratios.size() - 1], er1, var_g.W, N*nu, var);
        }
    }

    ZP.free_them_all();

    return vol;

}

#endif
