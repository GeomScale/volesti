// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

#ifndef BALL_ANN_VOL_H
#define BALL_ANN_VOL_H

#include <boost/math/distributions/students_t.hpp>
#include <boost/math/special_functions/erf.hpp>
#include "ball_annealingGl.h"
#include "esti_ratioGl.h"

template <class Polytope, class Point, class UParameters, class AParameters, typename NT>
NT volesti_ball_ann(Polytope &P, std::pair<Point,NT> &InnerBall, UParameters &var, AParameters &var_ban){

    typedef Ball<Point> ball;
    typedef BallIntersectPolytope<Polytope,ball> PolyBall;
    typedef typename UParameters::RNGType RNGType;
    typedef typename Polytope::VT VT;
    typedef std::list<Point> PointList;
    NT e = var.error;
    int n = var.n;
    bool verbose = var.verbose, round=var.round;

    NT lb = var_ban.lb, ub = var_ban.ub, PR = var_ban.p, B0_radius = var_ban.Bo_radius,
            ratio_B0 = var_ban.ratio_B0, rmax = var_ban.rmax;
    int win_len = var_ban.win_len, N = var_ban.N, nu = var_ban.nu;
    bool window2 = var_ban.window2;


    std::vector<ball> BallSet;
    std::vector<PointList> PointSets;
    std::vector<NT> ratios;
    NT p_value = 0.1;
    Point c = InnerBall.first;
    NT radius = InnerBall.second, round_value = 1.0;

    if(round){
#ifdef VOLESTI_DEBUG
        if(verbose) std::cout<<"\nRounding.."<<std::endl;
#endif
        double tstart1 = (double)clock()/(double)CLOCKS_PER_SEC;
        std::pair<NT,NT> res_round = rounding_min_ellipsoid(P,InnerBall,var);
        double tstop1 = (double)clock()/(double)CLOCKS_PER_SEC;
#ifdef VOLESTI_DEBUG
        if(verbose) std::cout << "Rounding time = " << tstop1 - tstart1 << std::endl;
#endif
       round_value=res_round.first;
       std::pair<Point,NT> res=P.ComputeInnerBall();
       c=res.first; radius=res.second;
    }

    // Save the radius of the Chebychev ball
    var.che_rad = radius;

    VT c_e(n);
    for(unsigned int i=0; i<n; i++){
        c_e(i)=c[i];  // write chebychev center in an eigen vector
    }
    P.shift(c_e);

    NT ratio0;//, steps, HnRSteps = 0.0, MemLps=0.0, ballsteps;
    ball B0;

    if(verbose) std::cout<<"Computing ball annealing..."<<std::endl;
    if(ratio_B0!=0.0) ratio0 = ratio_B0;
    NT ballsteps, steps;

    std::vector<std::vector<NT> > all_ratios;
    get_sequence_of_zonoballs<PolyBall, RNGType>(P, BallSet, B0, ratio0,
                                                 ratios, N, nu, lb, ub, radius, var,
                                                 ballsteps, steps, all_ratios, B0_radius, rmax);

    NT vol = (std::pow(M_PI,n/2.0)*(std::pow(B0.radius(), n) ) ) / (tgamma(n/2.0+1));

    int mm=BallSet.size()+2;
    NT prob = std::pow(PR, 1.0/NT(mm));
    NT er0 = e/(2.0*std::sqrt(NT(mm)));
    NT er1 = (e*std::sqrt(4.0*NT(mm)-1))/(2.0*std::sqrt(NT(mm)));

    int WW = 4*n*n+500;
    vol *= (window2) ? esti_ratio2<RNGType>(B0, P, er0, WW, ratio0, steps) :
            esti_ratio2_const<RNGType>(B0, P, er0, win_len, ratio0, prob, steps);
    //if(window2){
        //vol = vol * esti_ratio2<RNGType>(B0, P, er0, WW, ratio0, steps);
    //} else {
        //vol = vol * esti_ratio2_const<RNGType>(B0, P, er0, win_len, ratio0, prob, steps);
    //}

    ball Biter;
    PolyBall zb1, zb2;

    NT tele_prod=1.0;
    if (BallSet.size()>0) {
        er1 = er1 / std::sqrt(NT(mm)-1.0);
        vol *= (window2) ? 1/esti_ratio_interval<Point>(P, BallSet[0], ratios[0], er1, win_len, prob, var, steps) :
                1/esti_ratio(P, BallSet[0], ratios[0], er1, WW, var, steps);
        //if(window2) {
            //vol = vol / esti_ratio_interval<Point>(P, BallSet[0], ratios[0], er1, win_len, prob, var, steps);
        //} else {
            //vol = vol / esti_ratio(P, BallSet[0], ratios[0], er1, WW, var, steps);
        //}

        //HnRSteps += steps;
        for (int i = 0; i < BallSet.size() - 1; ++i) {
            zb1 = PolyBall(P, BallSet[i]);
            vol *= (window2) ? 1/esti_ratio_interval<Point>(zb1, BallSet[i+1], ratios[i+1], er1, win_len, prob, var, steps) :
                    1/esti_ratio(zb1, BallSet[i+1], ratios[i+1], er1, WW, var, steps);
            //if(window2) {
                //vol = vol / esti_ratio_interval<Point>(zb1, BallSet[i+1], ratios[i+1], er1, win_len, prob, var, steps);
            //} else {
                //vol = vol / esti_ratio(zb1, BallSet[i+1], ratios[i+1], er1, WW, var, steps);
            //}
        }


        zb1 = PolyBall(P, BallSet[BallSet.size() - 1]);
        vol *= (window2) ? 1/esti_ratio_interval<Point>(zb1, B0, ratios[ratios.size() - 1], er1, win_len, prob, var, steps) :
                1/esti_ratio(zb1, B0, ratios[ratios.size() - 1], er1, WW, var, steps);
        //if(window2) {
            //vol = vol / esti_ratio_interval<Point>(zb1, B0, ratios[ratios.size() - 1], er1, win_len, prob, var, steps);
        //} else {
            //vol = vol / esti_ratio(zb1, B0, ratios[ratios.size() - 1], er1, WW, var, steps);
        //}
    } else {
        if (ratios[0]!=1) {
            vol *= (window2) ? 1/esti_ratio_interval<Point>(P, B0, ratios[0], er1, win_len, prob, var, steps) :
                    1/esti_ratio(P, B0, ratios[0], er1, WW, var, steps);
            //if(window2) {
                //vol = vol / esti_ratio_interval<Point>(P, B0, ratios[0], er1, win_len, prob, var, steps);
            //} else {
                //vol = vol / esti_ratio(P, B0, ratios[0], er1, WW, var, steps);
            //}
        }
    }

    return vol*round_value;

}


#endif
