// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

#ifndef HVOL_VPOLY_H
#define HVOL_VPOLY_H

//#include "low_dimensional_sampling.h"
#include "ball_annealingGl.h"
#include "hpoly_annealing.h"
#include "esti_ratioGl.h"
#include "ZonoIntersectHPoly.h"
#include "hpoly_mmc_vp.h"


template <class Hpolytope, class Vpolytope, class UParameters, class AParameters, class GParameters, class Point, typename NT>
NT hvol_vpoly (Vpolytope &VP, UParameters &var, AParameters &var_ban, GParameters &var_g,
              std::pair<Point,NT> &InnerB, NT &nballs, int k = 0, bool only_balls = false) {

    typedef typename Vpolytope::VT VT;
    typedef typename Vpolytope::MT MT;
    typedef typename UParameters::RNGType RNGType;
    typedef Ball<Point> ball1;
    typedef BallIntersectPolytope<Hpolytope,ball1> BallPoly;

    NT lb = var_ban.lb, ub = var_ban.ub, prob = var_ban.p, rmax = var_ban.rmax, diam0 = var.diameter,
            radius = InnerB.second, round_value = 1.0, e = var.error, ratio, vol, alpha = var_ban.alpha;
    int n = var.n, win_len = var_ban.win_len, N = var_ban.N, nu = var_ban.nu, walkL = var.walk_steps;

    //var.verbose = true;
    bool verbose = var.verbose, round = var.round, window2 = var_ban.window2;

    if (k == 0) k=2*VP.num_of_vertices();

    UParameters var3 = var;
    var3.cdhr_walk = true;
    var3.ball_walk = false;
    var3.rdhr_walk = false;
    var3.bill_walk = false;

    Hpolytope HP(n);
    ball1 B0;
    enclosing_ball(VP, B0, var3);
    std::cout<<"B0.rad = "<<B0.radius()<<std::endl;
    BallPoly BP(HP, B0);
    construct_hpoly(VP, HP, BP, 10+n, k, var3);
    Hpolytope HP3;
    HP3.init(n,HP.get_mat(),HP.get_vec());
    //HP3.print();
    get_first_poly(VP, HP3, lb, ub, ratio, var3);
    HP3.normalize();

    var.MemLps = var.MemLps + var3.MemLps;
   // Hpolytope HP2=HP;

    //std::cout<<"facets of HP = "<<HP.num_of_hyperplanes()<<"\nA="<<HP.get_mat()<<"\n b="<<HP.get_vec()<<std::endl;
    std::pair<Point, NT> InnerBall = HP3.ComputeInnerBall();
    //std::cout<<"Cheb center = "<<std::endl;
    //for (int j = 0; j < n; ++j) {
        //std::cout<<InnerBall.first[j]<<" ";
    //}
    //std::cout<<"\n";

    typedef ZonoIntersectHPoly<Vpolytope, Hpolytope > ZonoHP;
    std::vector<Hpolytope > HPolySet;
    std::vector<NT> ratios;
    ZonoHP zb1, zb2;
    std::vector<NT> diams_inter;

    //var3.walk_steps = 1;

    if(verbose) std::cout<<"computing hpoly annealing.. = "<<vol<<"\n"<<std::endl;
    std::cout<<"N = "<<N<<" nu = "<< nu<<std::endl;
    get_sequence_of_vpoly_hpolys<ZonoHP>(VP, HP3, HPolySet, ratios, N*nu, nu, lb, ub, alpha, var, var3, diams_inter);
    var.diameter = diam0;
    nballs = NT(HPolySet.size()+1);
    if (only_balls){
        VP.free_them_all();
        return vol;
    }

    int mm=HPolySet.size()+2;
    int mm2=mm+1;
    prob = std::pow(prob, 1.0/NT(mm2));
    NT er0 = e/(2.0*std::sqrt(NT(mm2)));
    NT er1 = (e*std::sqrt(2.0*NT(mm2)-1))/(std::sqrt(2.0*NT(mm2)));
    NT Her = e/(2.0*std::sqrt(NT(mm2)));

    var_g.error = Her/2.0;
    std::cout<<"computing vol of h-polytope = "<<vol<<std::endl;
    NT fake_nballs;
    //BallPoly BP2(HP, B0);
    //std::pair<Point, NT> InnerBall2 = BP2.ComputeInnerBall();
    vol = volume_gaussian_annealing(HP3, var_g, var, InnerBall, fake_nballs);

    if(verbose) std::cout<<"\nvol of h-polytope = "<<vol<<"\n"<<std::endl;
    if (!window2) {
        UParameters var2 = var;
        var2.cdhr_walk = true;
        var2.ball_walk = false;
        var2.rdhr_walk = false;
        var2.bill_walk = false;
        var2.walk_steps = 10+2*n;
        vol *= esti_ratio_interval<RNGType, Point>(HP3, VP, ratio, er0, win_len, N*nu, prob, var2);
    } else {
        vol *= esti_ratio<RNGType, Point>(HP3, VP, ratio, er0, var_g.W, N*nu, var);
    }

    Hpolytope b1, b2;
    if (HPolySet.size()==0) {
        if (verbose) std::cout << "no hpoly | ratio = " << ratios[0] << std::endl;
        if (ratios[0]!=1) {
            if(!window2) {
                vol = vol / esti_ratio_interval<RNGType, Point>(VP, HP3, ratios[0], er1, win_len, N*nu, prob, var);
            } else {
                vol = vol / esti_ratio<RNGType, Point>(VP, HP3, ratios[0], er1, var_g.W, N*nu, var);
            }
        }
    } else {
        er1 = er1 / std::sqrt(NT(mm)-1.0);
        if(verbose) std::cout<<"number of hpolys = "<<HPolySet.size()<<std::endl;
        b1 = HPolySet[0];
        if(!window2) {
            vol = vol / esti_ratio_interval<RNGType, Point>(VP, b1, ratios[0], er1, win_len, N*nu, prob, var);
        } else {
            vol = vol / esti_ratio<RNGType, Point>(VP, b1, ratios[0], er1, var_g.W, N*nu, var);
        }

        for (int i = 0; i < HPolySet.size()-1; ++i) {
            zb1 = ZonoHP(VP,HPolySet[i]);
            b2 = HPolySet[i+1];
            var.diameter = diams_inter[i];
            if(!window2) {
                vol = vol / esti_ratio_interval<RNGType, Point>(zb1, b2, ratios[i], er1, win_len, N*nu, prob, var);
            } else {
                vol = vol / esti_ratio<RNGType, Point>(zb1, b2, ratios[i], er1, var_g.W, N*nu, var);
            }
        }

        zb1 = ZonoHP(VP,HPolySet[HPolySet.size()-1]);
        var.diameter = diams_inter[diams_inter.size()-1];
        if (!window2) {
            vol = vol / esti_ratio_interval<RNGType, Point>(zb1, HP3, ratios[ratios.size() - 1], er1, win_len, N*nu, prob, var);
        } else {
            vol = vol / esti_ratio<RNGType, Point>(zb1, HP3, ratios[ratios.size() - 1], er1, var_g.W, N*nu, var);
        }
    }

    VP.free_them_all();

    return vol;

}

#endif
