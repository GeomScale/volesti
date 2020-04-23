#ifndef BALL_ANN_VOL_H
#define BALL_ANN_VOL_H

//#include "cartesian_geom/cartesian_kernel.h"
//#include "vars.h"
//#include "hpolytope.h"
//#include "vpolytope.h"
//#include "zpolytope.h"
//#include "spectrahedron.h"
#include "ballintersectspec.h"
//#include "ballintersectconvex.h"
//#include "vpolyintersectvpoly.h"
//#include "samplers.h"
//#include "rounding.h"
#include "boost/math/distributions/students_t.hpp"
#include "boost/math/special_functions/erf.hpp"
#include "ball_annealingGl.h"
#include "esti_ratioGl.h"

template <class Spectrahedron, class Point, class UParameters, class AParameters, class SpecSettings, typename NT>
NT volesti_ball_ann(Spectrahedron &P, UParameters &var, AParameters &var_ban, SpecSettings &settings, std::pair<Point,NT> &InnerBall, NT &nballs,
                     bool only_phases = false) {

    typedef Ball <Point> ball;
    typedef BallIntersectSpectra <Spectrahedron, ball> PolyBall;
    typedef typename UParameters::RNGType RNGType;
    typedef typename Spectrahedron::VT VT;
    typedef typename Spectrahedron::MT MT;
    typedef std::list <Point> PointList;

    int n = var.n, win_len = var_ban.win_len, N = var_ban.N, nu = var_ban.nu;
    var.verbose = false;
    bool verbose = var.verbose, round = var.round, window2 = var_ban.window2;
    //verbose =true, var.verbose = true;
    NT lb = var_ban.lb, ub = var_ban.ub, prob = var_ban.p, rmax = var_ban.rmax, radius = InnerBall.second,
            round_value = 1.0, e = var.error, alpha = var_ban.alpha, diam = var.diameter;

    std::vector <ball> BallSet;
    std::vector <NT> ratios;
    Point c = InnerBall.first;

    /*
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
        var.diameter = 2.0 * std::sqrt(NT(n))*radius;
        P.comp_diam(var.diameter);
        diam = var.diameter;
    }*/

    // Save the radius of the Chebychev ball
    var.che_rad = radius;
    //VT c_e(n);
    //for (unsigned int i = 0; i < n; i++) {
    //    c_e(i) = c[i];  // write chebychev center in an eigen vector
    //}
    P.shift(c.getCoefficients()); // TODO we need to shift the spectrahedron s.t. the center of balls to be the origin
    //P.normalize();

    std::cout << "\nComputing ball annealing..." << std::endl;
    std::cout<<"N = "<<N<<" nu = "<< nu<<" e = "<<e<<std::endl;

    //MT LMIatP(P.getLMI().getMatricesDim(), P.getLMI().getMatricesDim());
    // TODO CHECK the settings!!!
    get_sequence_of_polyballs<PolyBall, RNGType>(P, BallSet, ratios, N * nu, nu, lb, ub, radius, alpha, var, settings, rmax); // we get the sequence of balls
    var.diameter = diam;

    NT vol = (std::pow(M_PI, n / 2.0) * (std::pow((*(BallSet.end() - 1)).radius(), n))) / (tgamma(n / 2.0 + 1));
    std::cout<<"rad of C_m = "<<(*(BallSet.end() - 1)).radius()<<", vol of B_m = "<<vol<<std::endl;

    int mm = BallSet.size() + 1;
    nballs = NT(mm - 1);
    if (only_phases) {
        //P.free_them_all();
        return vol;
    }
    prob = std::pow(prob, 1.0 / NT(mm));
    NT er0 = e / (2.0 * std::sqrt(NT(mm))), er1 = (e * std::sqrt(4.0 * NT(mm) - 1)) / (2.0 * std::sqrt(NT(mm)));

    std::cout<<"Estimating the ratios...\n"<<std::endl;
    P.set_LMIatP_A0(settings);
    std::cout<<"ratio_m = "<<*(ratios.end() - 1)<<std::endl;
    vol *= (window2) ? esti_ratio<RNGType, Point>(*(BallSet.end() - 1), P, *(ratios.end() - 1), er0, win_len, 1200, var,
            true, (*(BallSet.end() - 1)).radius()) :
           esti_ratio_interval<RNGType, Point>(*(BallSet.end() - 1), P, *(ratios.end() - 1), er0, win_len, 1200, prob,
                                               var, settings, true, (*(BallSet.end() - 1)).radius());

    PolyBall Pb;
    typename std::vector<ball>::iterator balliter = BallSet.begin();
    typename std::vector<NT>::iterator ratioiter = ratios.begin();

    er1 = er1 / std::sqrt(NT(mm) - 1.0);

    if (*ratioiter != 1) vol *= (!window2) ? 1 / esti_ratio_interval<RNGType, Point>(P, *balliter, *ratioiter, er1,
            win_len, N * nu, prob, var, settings) : 1 / esti_ratio<RNGType, Point>(P, *balliter, *ratioiter, er1, win_len, N * nu,
                                                                         var);
    for ( ; balliter < BallSet.end() - 1; ++balliter, ++ratioiter) {
        Pb = PolyBall(P, *balliter);
        Pb.comp_diam(var.diameter);
        Pb.set_LMIatP_A0(settings);
        vol *= (!window2) ? 1 / esti_ratio_interval<RNGType, Point>(Pb, *(balliter + 1), *(ratioiter + 1), er1,
                win_len, N * nu, prob, var, settings) : 1 / esti_ratio<RNGType, Point>(Pb, *balliter, *ratioiter, er1,
                                                                             win_len, N * nu, var);
    }

    //P.free_them_all();
    return vol * round_value;

}

#endif
