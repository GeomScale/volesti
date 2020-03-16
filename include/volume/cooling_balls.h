// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2019 Vissarion Fisikopoulos
// Copyright (c) 2018-2019 Apostolos Chalkis

//Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.

#ifndef COOLING_BALLS_H
#define COOLING_BALLS_H

#include "cartesian_geom/cartesian_kernel.h"
#include "vars.h"
#include "hpolytope.h"
#include "vpolytope.h"
#include "zpolytope.h"
#include "ballintersectconvex.h"
#include "vpolyintersectvpoly.h"
#include "samplers.h"
#include "rounding.h"
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/special_functions/erf.hpp>
#include "ball_annealing.h"
#include "ratio_estimation.h"

template <typename Polytope, typename Point, typename UParameters, typename AParameters, typename NT>
NT vol_cooling_balls(Polytope &P, UParameters &var, AParameters &var_ban, std::pair<Point,NT> &InnerBall) {

    typedef Ball <Point> ball;
    typedef BallIntersectPolytope <Polytope, ball> PolyBall;
    typedef typename UParameters::RNGType RNGType;
    typedef typename Polytope::VT VT;
    typedef std::list <Point> PointList;

    int n = var.n, win_len = var_ban.win_len, N = var_ban.N, nu = var_ban.nu;
    bool verbose = var.verbose, round = var.round, window2 = var_ban.window2;
    NT lb = var_ban.lb, ub = var_ban.ub, prob = var_ban.p, rmax = var_ban.rmax, radius = InnerBall.second,
            round_value = 1.0, e = var.error, alpha = var_ban.alpha, diam = var.diameter;

    std::vector <ball> BallSet;
    std::vector <NT> ratios;
    Point c = InnerBall.first;
    P.normalize();

    if (round) {
#ifdef VOLESTI_DEBUG
        if(verbose) std::cout<<"\nRounding.."<<std::endl;
#endif
        double tstart1 = (double) clock() / (double) CLOCKS_PER_SEC;
        std::pair <NT, NT> res_round = rounding_min_ellipsoid(P, InnerBall, var);
        double tstop1 = (double) clock() / (double) CLOCKS_PER_SEC;
#ifdef VOLESTI_DEBUG
        if(verbose) std::cout << "Rounding time = " << tstop1 - tstart1 << std::endl;
#endif
        round_value = res_round.first;
        std::pair <Point, NT> res = P.ComputeInnerBall();
        c = res.first;
        radius = res.second;
        P.normalize();
        if (var.bill_walk) {
            P.comp_diam(var.diameter, radius);
        }
        if (var.ball_walk){
            var.delta = 4.0 * radius / NT(n);
        }
    }

    // Save the radius of the Chebychev ball
    var.che_rad = radius;
    // Move the chebychev center to the origin and apply the same shifting to the polytope

    VT c_e = c.getCoefficients();
    P.shift(c_e);


    if ( !get_sequence_of_polyballs<PolyBall, RNGType>(P, BallSet, ratios, N * nu, nu, lb, ub, radius, alpha, var, rmax) ){
        return -1.0;
    }
    var.diameter = diam;

    NT vol = (std::pow(M_PI, n / 2.0) * (std::pow((*(BallSet.end() - 1)).radius(), n))) / (tgamma(n / 2.0 + 1));

    int mm = BallSet.size() + 1;
    prob = std::pow(prob, 1.0 / NT(mm));
    NT er0 = e / (2.0 * std::sqrt(NT(mm))), er1 = (e * std::sqrt(4.0 * NT(mm) - 1)) / (2.0 * std::sqrt(NT(mm)));

    vol *= (window2) ? esti_ratio<RNGType, Point>(*(BallSet.end() - 1), P, *(ratios.end() - 1), er0, win_len, 1200, var,
            true, (*(BallSet.end() - 1)).radius()) :
           esti_ratio_interval<RNGType, Point>(*(BallSet.end() - 1), P, *(ratios.end() - 1), er0, win_len, 1200, prob,
                                               var, true, (*(BallSet.end() - 1)).radius());

    PolyBall Pb;
    typename std::vector<ball>::iterator balliter = BallSet.begin();
    typename std::vector<NT>::iterator ratioiter = ratios.begin();

    er1 = er1 / std::sqrt(NT(mm) - 1.0);

    if (*ratioiter != 1) vol *= (!window2) ? 1 / esti_ratio_interval<RNGType, Point>(P, *balliter, *ratioiter, er1,
            win_len, N * nu, prob, var) : 1 / esti_ratio<RNGType, Point>(P, *balliter, *ratioiter, er1, win_len, N * nu,
                                                                         var);
    for ( ; balliter < BallSet.end() - 1; ++balliter, ++ratioiter) {
        Pb = PolyBall(P, *balliter);
        Pb.comp_diam(var.diameter, 0.0);
        vol *= (!window2) ? 1 / esti_ratio_interval<RNGType, Point>(Pb, *(balliter + 1), *(ratioiter + 1), er1,
                win_len, N * nu, prob, var) : 1 / esti_ratio<RNGType, Point>(Pb, *balliter, *ratioiter, er1,
                                                                             win_len, N * nu, var);
    }

    P.free_them_all();
    return vol * round_value;

}

#endif
