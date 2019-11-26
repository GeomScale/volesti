// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

#ifndef HZONO_VOL_H
#define HZONO_VOL_H

#include "ball_annealingGl.h"
#include "hpoly_annealing.h"
#include "esti_ratioGl.h"
#include "ZonoIntersectHPoly.h"


template <class Hpolytope, class Zonotope, class UParameters, class AParameters, class GParameters, class Point, typename NT>
NT vol_hzono (Zonotope &ZP, UParameters &var, AParameters &var_ban, GParameters &var_g,
              std::pair<Point,NT> &InnerB, NT &nballs, bool only_phases = false) {

    typedef typename Zonotope::VT VT;
    typedef typename Zonotope::MT MT;
    typedef typename UParameters::RNGType RNGType;

    NT lb = var_ban.lb, ub = var_ban.ub, prob = var_ban.p, rmax = var_ban.rmax, diam0 = var.diameter,
            radius = InnerB.second, round_value = 1.0, e = var.error, ratio, vol, alpha = var_ban.alpha;
    int n = var.n, win_len = var_ban.win_len, N = var_ban.N, nu = var_ban.nu;

    //var.verbose = true;
    bool verbose = var.verbose, round = var.round, window2 = var_ban.window2;

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

    //MT A3;
    //NT volh;
    MT T = ZP.get_T();
    MT Tt = T.transpose();
    MT A2 = AA * Tt;
    MT B = G * Tt;
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
    if(verbose) std::cout<<"get first hpoly.. = "<<vol<<"\n"<<std::endl;
    UParameters var3 = var;
    var3.cdhr_walk = true;
    var3.ball_walk = false;
    var3.rdhr_walk = false;
    var3.bill_walk = false;
    get_hdelta(ZP, HP, Zs_max, lb, ub, ratio, var3);
    var.MemLps = var.MemLps + var3.MemLps;
    Hpolytope HP2=HP;
    HP2.normalize();

    std::pair<Point, NT> InnerBall = HP.ComputeInnerBall();

    typedef ZonoIntersectHPoly<Zonotope , Hpolytope > ZonoHP;
    std::vector<Hpolytope > HPolySet;
    std::vector<NT> ratios;
    NT p_value = 0.1;
    ZonoHP zb1, zb2;

    std::vector<NT> diams_inter;

    //var3.walk_steps = 1;

    if(verbose) std::cout<<"computing hpoly annealing.. = "<<vol<<"\n"<<std::endl;
    std::cout<<"N = "<<N<<" nu = "<< nu<<std::endl;
    get_sequence_of_zonopolys<ZonoHP>(ZP, HP2, HPolySet, Zs_max, ratios, N*nu, nu, lb, ub, alpha, var, var3, diams_inter);
    var.diameter = diam0;
    nballs = NT(HPolySet.size()+1);
    if (only_phases){
        ZP.free_them_all();
        return vol;
    }

    int mm=HPolySet.size()+2;
    int mm2=mm+1;
    prob = std::pow(prob, 1.0/NT(mm2));
    NT er0 = e/(2.0*std::sqrt(NT(mm2)));
    NT er1 = (e*std::sqrt(2.0*NT(mm2)-1))/(std::sqrt(2.0*NT(mm2)));
    NT Her = e/(2.0*std::sqrt(NT(mm2)));

    var_g.error = Her;
    std::cout<<"computing vol of h-polytope = "<<vol<<std::endl;
    NT fake_nballs;
    vol = volume_gaussian_annealing(HP, var_g, var, InnerBall, fake_nballs);

    if(verbose) std::cout<<"\nvol of h-polytope = "<<vol<<"\n"<<std::endl;
    if (!window2) {
        UParameters var2 = var;
        var2.cdhr_walk = true;
        var2.ball_walk = false;
        var2.rdhr_walk = false;
        var2.bill_walk = false;
        var2.walk_steps = 10+2*n;
        vol *= esti_ratio_interval<RNGType, Point>(HP2, ZP, ratio, er0, win_len, N*nu, prob, var2);
    } else {
        vol *= esti_ratio<RNGType, Point>(HP2, ZP, ratio, er0, var_g.W, N*nu, var);
    }

    Hpolytope b1, b2;
    if (HPolySet.size()==0) {
        if (verbose) std::cout << "no hpoly | ratio = " << ratios[0] << std::endl;
        if (ratios[0]!=1) {
            if(!window2) {
                vol = vol / esti_ratio_interval<RNGType, Point>(ZP, HP2, ratios[0], er1, win_len, N*nu, prob, var);
            } else {
                vol = vol / esti_ratio<RNGType, Point>(ZP, HP2, ratios[0], er1, var_g.W, N*nu, var);
            }
        }
    } else {
        er1 = er1 / std::sqrt(NT(mm)-1.0);
        if(verbose) std::cout<<"number of hpolys = "<<HPolySet.size()<<std::endl;
        b1 = HPolySet[0];
        if(!window2) {
            vol = vol / esti_ratio_interval<RNGType, Point>(ZP, b1, ratios[0], er1, win_len, N*nu, prob, var);
        } else {
            vol = vol / esti_ratio<RNGType, Point>(ZP, b1, ratios[0], er1, var_g.W, N*nu, var);
        }

        for (int i = 0; i < HPolySet.size()-1; ++i) {
            zb1 = ZonoHP(ZP,HPolySet[i]);
            b2 = HPolySet[i+1];
            var.diameter = diams_inter[i];
            if(!window2) {
                vol = vol / esti_ratio_interval<RNGType, Point>(zb1, b2, ratios[i], er1, win_len, N*nu, prob, var);
            } else {
                vol = vol / esti_ratio<RNGType, Point>(zb1, b2, ratios[i], er1, var_g.W, N*nu, var);
            }
        }

        zb1 = ZonoHP(ZP,HPolySet[HPolySet.size()-1]);
        var.diameter = diams_inter[diams_inter.size()-1];
        if (!window2) {
            vol = vol / esti_ratio_interval<RNGType, Point>(zb1, HP2, ratios[ratios.size() - 1], er1, win_len, N*nu, prob, var);
        } else {
            vol = vol / esti_ratio<RNGType, Point>(zb1, HP2, ratios[ratios.size() - 1], er1, var_g.W, N*nu, var);
        }
    }

    ZP.free_them_all();

    return vol;

}

#endif
