

// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

#ifndef BALL_ANN_VOL_H
#define BALL_ANN_VOL_H

#include "ball_annealingGl.h"
#include "esti_ratioGl.h"

template <class Polytope, class Point, class Parameters, typename NT>
NT volesti_ball_ann(Polytope &P, std::pair<Point,NT> &InnerBall, NT &lb, NT &up, Parameters &var, NT &hnrst, NT &nballs, NT &memball, bool steps_only = false){

    typedef Ball<Point> ball;
    typedef BallIntersectPolytope<Polytope,ball> ZonoBall;
    typedef typename Parameters::RNGType RNGType;
    typedef typename Polytope::VT VT;
    NT e = var.error;
    int n = var.n;
    bool verbose = var.verbose, round=var.round;

    typedef std::list<Point> PointList;
    std::vector<ball> BallSet;
    std::vector<PointList> PointSets;
    std::vector<NT> ratios;
    NT p_value = 0.1;
    Point c = InnerBall.first;
    NT radius, round_value = 1.0;

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

    //Point c = InnerBall.first;
    VT c_e(n);
    for(unsigned int i=0; i<n; i++){
        c_e(i)=c[i];  // write chebychev center in an eigen vector
    }
    P.shift(c_e);

    NT ratio0, steps, HnRSteps = 0.0, MemLps=0.0, ballsteps;
   // var.coordinate=true;
    ball B0;

    if(verbose) std::cout<<"Computing ball annealing..."<<std::endl;
    get_sequence_of_zonoballs<ZonoBall, RNGType>(P, BallSet, B0, ratio0, PointSets,
                                                 ratios, lb, up, InnerBall.second, var,
                                                 ballsteps, steps);
    if(steps_only) {
        return NT(BallSet.size()+1);
    }
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
    vol = vol * esti_ratio2<RNGType>(B0, P, e/3.0, ratio0, steps);
    MemLps += steps;

    ball Biter;
    ZonoBall zb1, zb2;
    //var.coordinate = true;
    var.walk_steps=1;
    NT tele_prod=1.0;
    NT er = e*0.942809;
    if (BallSet.size()>0) {
        er = er / std::sqrt(NT(BallSet.size()+1));
        tele_prod = tele_prod * esti_ratio(P, BallSet[0], ratios[0], er, var, steps);
        HnRSteps += steps;
        for (int i = 0; i < BallSet.size() - 1; ++i) {
            zb1 = ZonoBall(P, BallSet[i]);
            //zb2 = ZonoBall(ZP, BallSet[i]);
            tele_prod = tele_prod * esti_ratio(zb1, BallSet[i+1], ratios[i+1], er, var, steps);
            HnRSteps += steps;
        }


        zb1 = ZonoBall(P, BallSet[BallSet.size() - 1]);
        tele_prod = tele_prod * esti_ratio(zb1, B0, ratios[ratios.size() - 1], er, var, steps);
        HnRSteps += steps;
        vol = vol / tele_prod;
    } else {
        if (ratios[0]!=1) {
            vol = vol / esti_ratio(P, B0, ratios[0], er, var, steps);
            HnRSteps += steps;
        }
    }

    //Rcpp::NumericVector res(4);
    //res[0] = vol;
    //res[1] = NT(BallSet.size()+1);
    //res[2] = HnRSteps;
    //res[3] = MemLps;
    if(verbose) {
        std::cout<<"number of balls = "<<BallSet.size()+1<<std::endl;
        std::cout<<"number of HnR steps = "<<HnRSteps<<std::endl;
        std::cout<<"number of memberships (ball random points) = "<<MemLps<<std::endl;
    }
    hnrst = HnRSteps;
    memball = MemLps;
    nballs = NT(BallSet.size()+1);

    return vol*round_value;


}


#endif
