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
//' @param random_walk Optional. A string that declares the random walk.
//' @param walk_length Optional. The number of the steps for the random walk.
//' @param parameters Optional. A list for the parameters of the methods:
//'
//' @section warning:
//' Do not use this function.
//'
//' @return A numerical matrix that describes the rounded polytope and contains the round value.
// [[Rcpp::export]]
Rcpp::List rounding (Rcpp::Reference P,
                              Rcpp::Nullable<std::string> random_walk = R_NilValue,
                              Rcpp::Nullable<unsigned int> walk_length = R_NilValue,
                     Rcpp::Nullable<Rcpp::List> parameters = R_NilValue){

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
            cdhr=true, rdhr = false, ball_walk = false, billiard = false;
    NT delta = -1.0, diam = -1.0;

    unsigned int n = P.field("dimension");
    unsigned int rnum = std::pow(1.0,-2.0) * 400 * n * std::log(n);
    unsigned int walkL = 10+n/10;

    std::pair <Point, NT> InnerBall;
    Rcpp::NumericMatrix Mat;
    int type = P.field("type");

    switch (type) {
        case 1: {
            // Hpolytope
            HP.init(n, Rcpp::as<MT>(P.field("A")), Rcpp::as<VT>(P.field("b")));
            InnerBall = HP.ComputeInnerBall();
            diam = 2.0 * std::sqrt(NT(n)) * InnerBall.second;
            delta = 4.0 * InnerBall.second / std::sqrt(NT(n));
            break;
        }
        case 2: {
            VP.init(n, Rcpp::as<MT>(P.field("V")), VT::Ones(Rcpp::as<MT>(P.field("V")).rows()));
            InnerBall = VP.ComputeInnerBall();
            VP.comp_diam(diam);
            delta = 4.0 * InnerBall.second / std::sqrt(NT(n));
            break;
        }
        case 3: {
            // Zonotope
            ZP.init(n, Rcpp::as<MT>(P.field("G")), VT::Ones(Rcpp::as<MT>(P.field("G")).rows()));
            InnerBall = ZP.ComputeInnerBall();
            ZP.comp_diam(diam);
            delta = 4.0 * InnerBall.second / std::sqrt(NT(n));
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

    if(!random_walk.isNotNull()) {
        if (type == 1) {
            cdhr = true;
        } else {
            billiard = true;
        }
    } else if(Rcpp::as<std::string>(random_walk).compare(std::string("CDHR"))==0) {
        cdhr = true;
    } else if (Rcpp::as<std::string>(random_walk).compare(std::string("RDHR"))==0) {
        rdhr = true;
    } else if (Rcpp::as<std::string>(random_walk).compare(std::string("BW"))==0) {
        if (Rcpp::as<Rcpp::List>(parameters).containsElementNamed("BW_rad")) {
            delta = Rcpp::as<NT>(Rcpp::as<Rcpp::List>(parameters)["BW_rad"]);
        }
        ball_walk = true;
    } else if (Rcpp::as<std::string>(random_walk).compare(std::string("BiW")) == 0) {
        billiard = true;
        walkL = 1;
        if (Rcpp::as<Rcpp::List>(parameters).containsElementNamed("diameter"))
            diam = Rcpp::as<NT>(Rcpp::as<Rcpp::List>(parameters)["diameter"]);
    } else {
        throw Rcpp::exception("Unknown walk type!");
    }

    if(walk_length.isNotNull()) walkL = Rcpp::as<unsigned int>(walk_length);

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    // the random engine with this seed
    RNGType rng(seed);
    boost::random::uniform_real_distribution<>(urdist);
    boost::random::uniform_real_distribution<> urdist1(-1,1);

    // initialization
    vars<NT, RNGType> var(rnum,n,walkL,1,0.0,0.0,0,0.0,0,InnerBall.second,diam,rng,urdist,urdist1,
                          delta,verbose,rand_only,false,NN,birk,ball_walk,cdhr,rdhr,billiard);
    std::pair <NT, NT> round_res;

    switch (type) {
        case 1: {
            round_res = rounding_min_ellipsoid(HP, InnerBall, var);
            Mat = extractMatPoly(HP);
            break;
        }
        case 2: {
            round_res = rounding_min_ellipsoid(VP, InnerBall, var);
            Mat = extractMatPoly(VP);
            break;
        }
        case 3: {
            round_res = rounding_min_ellipsoid(ZP, InnerBall, var);
            Mat = extractMatPoly(ZP);
            break;
        }
    }

    return Rcpp::List::create(Rcpp::Named("Mat") = Mat , Rcpp::Named("round_value") = round_res.first);

}
