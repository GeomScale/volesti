// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

#ifndef BALL_ANN_VOL_H
#define BALL_ANN_VOL_H

#include "cartesian_geom/cartesian_kernel.h"
#include "vars.h"
#include "polytopes.h"
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

    int n = var.n, win_len = var_ban.win_len, N = var_ban.N, nu = var_ban.nu;
    var.verbose = true;
    bool verbose = var.verbose, round = var.round, window2 = var_ban.window2;
    NT lb = var_ban.lb, ub = var_ban.ub, prob = var_ban.p, B0_radius = var_ban.B0_radius, rmax = var_ban.rmax,
            radius = InnerBall.second, round_value = 1.0, tele_prod = 1.0, e = var.error, ratio0;

    std::vector <ball> BallSet;
    std::vector <PointList> PointSets;
    std::vector <NT> ratios;
    ball B0;
    Point c = InnerBall.first;

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

    if (verbose) std::cout << "Computing ball annealing..." << std::endl;
    get_sequence_of_zonoballs<PolyBall, RNGType>(P, BallSet, B0, ratio0,
                                                 ratios, N * nu, nu, lb, ub, radius, var);

    NT vol = (std::pow(M_PI, n / 2.0) * (std::pow(B0.radius(), n))) / (tgamma(n / 2.0 + 1));

    int mm = BallSet.size() + 2;
    prob = std::pow(prob, 1.0 / NT(mm));
    NT er0 = e / (2.0 * std::sqrt(NT(mm))), er1 = (e * std::sqrt(4.0 * NT(mm) - 1)) / (2.0 * std::sqrt(NT(mm)));

    vol *= (window2) ? esti_ratio<RNGType, Point>(B0, P, ratio0, er0, win_len, 1200, var, true, B0.radius()) :
           esti_ratio_interval<RNGType, Point>(B0, P,  ratio0, er0, win_len, 1200, prob, var, true, B0.radius());

    PolyBall Pb;
    typename std::vector<ball>::iterator balliter = BallSet.begin();
    typename std::vector<NT>::iterator ratioiter = ratios.begin();

    if (BallSet.size() > 0) {
        er1 = er1 / std::sqrt(NT(mm) - 1.0);
        vol *= (!window2) ? 1 / esti_ratio_interval<RNGType, Point>(P, *balliter, *ratioiter, er1, win_len, N * nu, prob,
                                                                var) :
               1 / esti_ratio<RNGType, Point>(P, BallSet[0], ratios[0], er1, win_len, N * nu, var);

        ++balliter;
        ++ratioiter;
        for ( ; balliter < BallSet.end(); ++balliter, ++ratioiter) {
            Pb = PolyBall(P, *(balliter-1));
            vol *= (!window2) ? 1 / esti_ratio_interval<RNGType, Point>(Pb, *balliter, *ratioiter, er1, win_len,
                                                                        N * nu, prob, var) :
                   1 / esti_ratio<RNGType, Point>(Pb, *balliter, *ratioiter, er1, win_len, N * nu, var);
        }

        Pb = PolyBall(P, *(BallSet.end() - 1));
        vol *= (!window2) ? 1 /
                            esti_ratio_interval<RNGType, Point>(Pb, B0, *(ratios.end() - 1), er1, win_len, N * nu, prob, var) :
               1 / esti_ratio<RNGType, Point>(Pb, B0, *(ratios.end() - 1), er1, win_len, N * nu, var);
    } else {
        if (*ratioiter != 1) {
            vol *= (!window2) ? 1 /
                                esti_ratio_interval<RNGType, Point>(P, B0, *ratioiter, er1, win_len, N * nu, prob, var) :
                   1 / esti_ratio<RNGType, Point>(P, B0, *ratioiter, er1, win_len, N * nu, var);
        }
    }

    return vol * round_value;

}

#endif
