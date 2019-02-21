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
#include "polytopes.h"
#include "samplers.h"
#include "rounding.h"
#include "extractMatPoly.h"

//' Internal rcpp function for the rounding of a convex polytope
//'
//' @param P A convex polytope (H- or V-representation or zonotope).
//' @param WalkType Optional. A string that declares the random walk.
//' @param walk_step Optional. The number of the steps for the random walk.
//' @param radius Optional. The radius for the ball walk.
//'
//' @section warning:
//' Do not use this function.
//'
//' @return A numerical matrix that describes the rounded polytope and contains the round value.
// [[Rcpp::export]]
Rcpp::NumericMatrix rounding (Rcpp::Reference P,
                              Rcpp::Nullable<std::string> WalkType = R_NilValue,
                              Rcpp::Nullable<unsigned int> walk_step = R_NilValue,
                              Rcpp::Nullable<double> radius = R_NilValue) {

    typedef double NT;
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef boost::mt19937    RNGType;
    typedef HPolytope<Point> Hpolytope;
    typedef VPolytope<Point, RNGType > Vpolytope;
    typedef Zonotope<Point> zonotope;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;

    Hpolytope HP;
    Vpolytope VP;
    zonotope ZP;

    bool rand_only=false,
            NN=false,
            birk=false,
            verbose=false,
            coordinate=true, ball_walk = false;
    NT delta = -1.0;

    unsigned int n = P.field("dimension");
    unsigned int rnum = std::pow(1.0,-2.0) * 400 * n * std::log(n);
    unsigned int walkL = 10+n/10;
    if(walk_step.isNotNull()) walkL = Rcpp::as<unsigned int>(walk_step);

    std::pair <Point, NT> InnerBall;
    Rcpp::NumericMatrix Mat;
    int type = P.field("type");

    switch (type) {
        case 1: {
            // Hpolytope
            HP.init(n, Rcpp::as<MT>(P.field("A")), Rcpp::as<VT>(P.field("b")));
            InnerBall = HP.ComputeInnerBall();
            break;
        }
        case 2: {
            VP.init(n, Rcpp::as<MT>(P.field("V")), VT::Ones(Rcpp::as<MT>(P.field("V")).rows()));
            InnerBall = VP.ComputeInnerBall();
            break;
        }
        case 3: {
            // Zonotope
            ZP.init(n, Rcpp::as<MT>(P.field("G")), VT::Ones(Rcpp::as<MT>(P.field("G")).rows()));
            InnerBall = ZP.ComputeInnerBall();
            break;
        }
        //default: throw Rcpp::exception("Wrong polytope input");
    }

    if(!WalkType.isNotNull() || Rcpp::as<std::string>(WalkType).compare(std::string("CDHR"))==0){
        coordinate = true;
        ball_walk = false;
    } else if (Rcpp::as<std::string>(WalkType).compare(std::string("RDHR"))==0) {
        coordinate = false;
        ball_walk = false;
    } else if (Rcpp::as<std::string>(WalkType).compare(std::string("BW"))==0) {
        if(radius.isNotNull()){
            delta = Rcpp::as<NT>(radius);
        } else {
            delta = 4.0 * InnerBall.second / std::sqrt(NT(n));
        }
        coordinate = false;
        ball_walk = true;
    } else {
        throw Rcpp::exception("Unknown walk type!");
    }

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    // the random engine with this seed
    RNGType rng(seed);
    boost::random::uniform_real_distribution<>(urdist);
    boost::random::uniform_real_distribution<> urdist1(-1,1);

    // initialization
    vars<NT, RNGType> var(rnum,n,walkL,1,0.0,0.0,0,0.0,0,InnerBall.second,rng,urdist,urdist1,
                          delta,verbose,rand_only,false,NN,birk,ball_walk,coordinate);
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

    /*
    if (type == 3) {
        round_res = rounding_min_ellipsoid(ZP, InnerBall, var);
        Mat = extractMatPoly(ZP);
    } else if (type == 1) {
        round_res = rounding_min_ellipsoid(HP, InnerBall, var);
        Mat = extractMatPoly(HP);
    } else {
        round_res = rounding_min_ellipsoid(VP, InnerBall, var);
        Mat = extractMatPoly(VP);
    }*/
    // store rounding value and the ratio between min and max axe in the first row
    // the matrix is in ine format so the first row is useless and is going to be removed by R function modifyMat()
    Mat(0,0) = round_res.first;
    Mat(0,1) = round_res.second;
    return Mat;

}
