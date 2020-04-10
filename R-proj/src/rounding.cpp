// [[Rcpp::depends(BH)]]

// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

#include <Rcpp.h>
#include <RcppEigen.h>
#include <chrono>
#include "cartesian_geom/cartesian_kernel.h"
#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include "vars.h"
#include "hpolytope.h"
#include "vpolytope.h"
#include "zpolytope.h"
#include "samplers.h"
#include "rounding.h"
#include "vpolyintersectvpoly.h"
#include "extractMatPoly.h"

//' Internal rcpp function for the rounding of a convex polytope
//'
//' @param P A convex polytope (H- or V-representation or zonotope).
//'
//' @section warning:
//' Do not use this function.
//'
//' @return A numerical matrix that describes the rounded polytope and contains the round value.
// [[Rcpp::export]]
Rcpp::List rounding (Rcpp::Reference P){

    typedef double NT;
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef boost::mt19937    RNGType;
    typedef HPolytope<Point> Hpolytope;
    typedef VPolytope<Point, RNGType > Vpolytope;
    typedef Zonotope<Point> zonotope;
    typedef IntersectionOfVpoly<Vpolytope> InterVP;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;

    Hpolytope HP;
    Vpolytope VP;
    zonotope ZP;
    InterVP VPcVP;

    bool rand_only=false,
            NN=false,
            birk=false,
            verbose=false,
            cdhr=false, rdhr = false, ball_walk = false, billiard = false;
    NT delta = -1.0, diam = -1.0;

    unsigned int n = P.field("dimension");
    unsigned int rnum = std::pow(1.0,-2.0) * 400 * n * std::log(n);
    unsigned int walkL = 10+n/10;

    std::pair <Point, NT> InnerBall;
    Rcpp::NumericMatrix Mat;
    int type = P.field("type");

    if (type == 1) {
        walkL = 10 + 10/n;
        cdhr = true;
    } else {
        walkL = 5;
        billiard = true;
    }

    switch (type) {
        case 1: {
            // Hpolytope
            HP.init(n, Rcpp::as<MT>(P.field("A")), Rcpp::as<VT>(P.field("b")));
            InnerBall = HP.ComputeInnerBall();
            //if (billiard && diam < 0.0) HP.comp_diam(diam, InnerBall.second);
            break;
        }
        case 2: {
            VP.init(n, Rcpp::as<MT>(P.field("V")), VT::Ones(Rcpp::as<MT>(P.field("V")).rows()));
            InnerBall = VP.ComputeInnerBall();
            VP.comp_diam(diam, 0.0);
            break;
        }
        case 3: {
            // Zonotope
            ZP.init(n, Rcpp::as<MT>(P.field("G")), VT::Ones(Rcpp::as<MT>(P.field("G")).rows()));
            InnerBall = ZP.ComputeInnerBall();
            ZP.comp_diam(diam, 0.0);
            break;
        }
        case 4: {
            throw Rcpp::exception("volesti does not support rounding for this representation currently.");
            /*
            Vpolytope VP1;
            Vpolytope VP2;
            VP1.init(n, Rcpp::as<MT>(Rcpp::as<Rcpp::Reference>(P).field("V1")),
                     VT::Ones(Rcpp::as<MT>(Rcpp::as<Rcpp::Reference>(P).field("V1")).rows()));
            VP2.init(n, Rcpp::as<MT>(Rcpp::as<Rcpp::Reference>(P).field("V2")),
                     VT::Ones(Rcpp::as<MT>(Rcpp::as<Rcpp::Reference>(P).field("V2")).rows()));
            VPcVP.init(VP1, VP2);

            if (!VPcVP.is_feasible()) throw Rcpp::exception("Empty set!");
            InnerBall = VPcVP.ComputeInnerBall();

            diam = 2.0 * std::sqrt(NT(n)) * InnerBall.second;
            VPcVP.comp_diam(diam);
            delta = 4.0 * InnerBall.second / std::sqrt(NT(n));*/

        }
        //default: throw Rcpp::exception("Wrong polytope input");
    }

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    // the random engine with this seed
    RNGType rng(seed);
    boost::random::uniform_real_distribution<>(urdist);
    boost::random::uniform_real_distribution<> urdist1(-1,1);

    // initialization
    vars<NT, RNGType> var(rnum,n,walkL,1,0.0,0.0,0,0.0,0,InnerBall.second,diam,rng,urdist,urdist1,
                          delta,verbose,rand_only,false,NN,birk,ball_walk,cdhr,rdhr,billiard);
    std::pair< std::pair<MT, VT>, NT > round_res;

    switch (type) {
        case 1: {
            round_res = rounding_min_ellipsoid<MT, VT>(HP, InnerBall, var);
            Mat = extractMatPoly(HP);
            break;
        }
        case 2: {
            round_res = rounding_min_ellipsoid<MT, VT>(VP, InnerBall, var);
            Mat = extractMatPoly(VP);
            break;
        }
        case 3: {
            round_res = rounding_min_ellipsoid<MT, VT>(ZP, InnerBall, var);
            Mat = extractMatPoly(ZP);
            break;
        }
    }

    return Rcpp::List::create(Rcpp::Named("Mat") = Mat, Rcpp::Named("T") = Rcpp::wrap(round_res.first.first),
                              Rcpp::Named("shift") = Rcpp::wrap(round_res.first.second),
                              Rcpp::Named("round_value") = round_res.second);

}
