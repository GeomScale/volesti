
// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis


#include <Rcpp.h>
#include <RcppEigen.h>
#include <chrono>
#include "cartesian_geom/cartesian_kernel.h"
#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include "polytopes.h"
#include "vpolyintersectvpoly.h"

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
Rcpp::NumericVector InnerBall(Rcpp::Reference P) {

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
    unsigned int n = P.field("dimension");

    std::pair <Point, NT> InnerBall;

    int type = P.field("type");
    if (type==1) {
        // Hpolytope
        Hpolytope HP;
        HP.init(n, Rcpp::as<MT>(P.field("A")), Rcpp::as<VT>(P.field("b")));
        InnerBall = HP.ComputeInnerBall();
    } else if(type==2) {
        // Vpolytope
        Vpolytope VP;
        VP.init(n, Rcpp::as<MT>(P.field("V")), VT::Ones(Rcpp::as<MT>(P.field("V")).rows()));
        InnerBall = VP.ComputeInnerBall();

    } else if(type==3){
        // Zonotope
        zonotope ZP;
        InnerBall = ZP.ComputeInnerBall();
        ZP.init(n, Rcpp::as<MT>(P.field("G")), VT::Ones(Rcpp::as<MT>(P.field("G")).rows()));
        InnerBall = ZP.ComputeInnerBall();
    } else {
        // Intersection of two V-polytopes
        Vpolytope VP1;
        Vpolytope VP2;
        InterVP VPcVP;
        VP1.init(n, Rcpp::as<MT>(P.field("V1")), VT::Ones(Rcpp::as<MT>(P.field("V1")).rows()));
        VP2.init(n, Rcpp::as<MT>(P.field("V2")), VT::Ones(Rcpp::as<MT>(P.field("V2")).rows()));
        VPcVP.init(VP1, VP2);
        InnerBall = VPcVP.ComputeInnerBall();
    }


    Rcpp::NumericVector vec(n + 1);
    for (unsigned int k = 0; k < n; ++k) {
        vec[k] = InnerBall.first[k];
    }
    vec[n] = InnerBall.second;
    return vec;

}