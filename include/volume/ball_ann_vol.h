// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

#ifndef BALL_ANN_VOL_H
#define BALL_ANN_VOL_H

//#include <iterator>
//#include <fstream>
//#include <vector>
//#include <list>
//#include <algorithm>
//#include <math.h>
//#include <chrono>
#include "cartesian_geom/cartesian_kernel.h"
//#include "random.hpp"
//#include "random/uniform_int.hpp"
//#include "random/normal_distribution.hpp"
//#include "random/uniform_real_distribution.hpp"
#include "vars.h"
#include "polytopes.h"
//#include "ellipsoids.h"
#include "ballintersectconvex.h"
#include "vpolyintersectvpoly.h"
#include "samplers.h"
#include "rounding.h"
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/special_functions/erf.hpp>
#include "ball_annealingGl.h"
#include "esti_ratioGl.h"

template <class Polytope, class Point, class UParameters, class AParameters, typename NT>
NT volesti_ball_ann(Polytope &P, UParameters &var, AParameters &var_ban, std::pair<Point,NT> &InnerBall) {

    typedef Ball <Point> ball;
    typedef BallIntersectPolytope <Polytope, ball> PolyBall;
    typedef typename UParameters::RNGType RNGType;
    typedef typename Polytope::VT VT;
    typedef std::list <Point> PointList;
    NT e = var.error;
    int n = var.n;
    var.verbose = true;
    bool verbose = var.verbose, round = var.round;

    NT lb = var_ban.lb, ub = var_ban.ub, PR = var_ban.p, B0_radius = var_ban.B0_radius,
            ratio_B0 = var_ban.ratio_B0, rmax = var_ban.rmax;
    int win_len = var_ban.win_len, N = var_ban.N, nu = var_ban.nu;
    bool window2 = var_ban.window2;


    std::vector <ball> BallSet;
    std::vector <PointList> PointSets;
    std::vector <NT> ratios;
    NT p_value = 0.1;
    Point c = InnerBall.first;
    NT radius = InnerBall.second, round_value = 1.0;

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
    }

    // Save the radius of the Chebychev ball
    var.che_rad = radius;

    VT c_e(n);
    for (unsigned int i = 0; i < n; i++) {
        c_e(i) = c[i];  // write chebychev center in an eigen vector
    }
    P.shift(c_e);

    NT ratio0;//, steps, HnRSteps = 0.0, MemLps=0.0, ballsteps;
    ball B0;

    if (verbose) std::cout << "Computing ball annealing..." << std::endl;
    if (ratio_B0 != 0.0) ratio0 = ratio_B0;

    get_sequence_of_zonoballs<PolyBall, RNGType>(P, BallSet, B0, ratio0,
                                                 ratios, N * nu, nu, lb, ub, radius, var);

    NT vol = (std::pow(M_PI, n / 2.0) * (std::pow(B0.radius(), n))) / (tgamma(n / 2.0 + 1));

    int mm = BallSet.size() + 2;
    NT prob = std::pow(PR, 1.0 / NT(mm));
    NT er0 = e / (2.0 * std::sqrt(NT(mm)));
    NT er1 = (e * std::sqrt(4.0 * NT(mm) - 1)) / (2.0 * std::sqrt(NT(mm)));

    int WW = 4 * n * n + 500;
    vol *= (window2) ? esti_ratio2<RNGType>(B0, P, er0, WW, ratio0) :
           esti_ratio_interval<RNGType, Point>(B0, P,  ratio0, er1, win_len, 1200, prob, var, true, B0.radius()); //esti_ratio2_const<RNGType>(B0, P, er0, win_len, ratio0, prob);

    //ball Biter;
    PolyBall Pb;

    NT tele_prod = 1.0;
    typename std::vector<ball> balliter = BallSet.begin();
    typename std::vector<NT> ratioiter = ratios.begin();
    if (BallSet.size() > 0) {
        er1 = er1 / std::sqrt(NT(mm) - 1.0);
        vol *= (!window2) ? 1 /
                            esti_ratio_interval<RNGType, Point>(P, BallSet[0], ratios[0], er1, win_len, N * nu, prob,
                                                                var) : 1 / esti_ratio(P, BallSet[0], ratios[0], er1, WW,
                                                                                      var);

        for (int i = 0; i < BallSet.size() - 1; ++i) {
            Pb = PolyBall(P, BallSet[i]);
            vol *= (!window2) ? 1 / esti_ratio_interval<RNGType, Point>(Pb, BallSet[i + 1], ratios[i + 1], er1, win_len,
                                                                        N * nu, prob, var) :
                   1 / esti_ratio(zb1, BallSet[i + 1], ratios[i + 1], er1, WW, var);
        }


        Pb = PolyBall(P, BallSet[BallSet.size() - 1]);
        vol *= (!window2) ? 1 /
                            esti_ratio_interval<RNGType, Point>(Pb, B0, ratios[ratios.size() - 1], er1, win_len, N * nu,
                                                                prob, var) :
               1 / esti_ratio(zb1, B0, ratios[ratios.size() - 1], er1, WW, var);
    } else {
        if (ratios[0] != 1) {
            vol *= (!window2) ? 1 /
                                esti_ratio_interval<RNGType, Point>(P, B0, ratios[0], er1, win_len, N * nu, prob, var) :
                   1 / esti_ratio(P, B0, ratios[0], er1, WW, var);
        }
    }

    return vol * round_value;

}


#endif
