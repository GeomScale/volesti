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
#include "hpolytope.h"
#include "vpolytope.h"
#include "zpolytope.h"
#include "vpolyintersectvpoly.h"
#include "samplers.h"
#include "rotating.h"
#include "extractMatPoly.h"

//'  An internal Rccp function for the random rotation of a convex polytope
//'
//' @param P A convex polytope (H-, V-polytope or a zonotope).
//' @param T Optional. A rotation matrix.
//' @param seed Optional. A fixed seed for the random linear map generator.
//'
//' @section warning:
//' Do not use this function.
//'
//' @return A matrix that describes the rotated polytope
// [[Rcpp::export]]
Rcpp::NumericMatrix rotating (Rcpp::Reference P, Rcpp::Nullable<Rcpp::NumericMatrix> T = R_NilValue,
                              Rcpp::Nullable<double> seed = R_NilValue){

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

    MT TransorfMat;
    Rcpp::NumericMatrix Mat;
    unsigned int n = P.field("dimension");
    int type = P.field("type");

    double seed2 = (!seed.isNotNull()) ? std::numeric_limits<double>::signaling_NaN() : Rcpp::as<double>(seed);

    switch (type) {
        case 1: {
            // Hpolytope
            Hpolytope HP;
            HP.init(n, Rcpp::as<MT>(P.field("A")), Rcpp::as<VT>(P.field("b")));
            if (T.isNotNull()) {
                TransorfMat = Rcpp::as<MT>(T);
                HP.linear_transformIt(TransorfMat.inverse());
            } else {
                TransorfMat = rotating < MT > (HP, seed2);
            }
            Mat = extractMatPoly(HP);
            break;
        }
        case 2: {
            // Vpolytope
            Vpolytope VP;
            VP.init(n, Rcpp::as<MT>(P.field("V")), VT::Ones(Rcpp::as<MT>(P.field("V")).rows()));
            if (T.isNotNull()) {
                TransorfMat = Rcpp::as<MT>(T);
                VP.linear_transformIt(TransorfMat.inverse());
            } else {
                TransorfMat = rotating < MT > (VP, seed2);
            }
            Mat = extractMatPoly(VP);
            break;
        }
        case 3: {
            // Zonotope
            zonotope ZP;
            ZP.init(n, Rcpp::as<MT>(P.field("G")), VT::Ones(Rcpp::as<MT>(P.field("G")).rows()));
            if (T.isNotNull()) {
                TransorfMat = Rcpp::as<MT>(T);
                ZP.linear_transformIt(TransorfMat.inverse());
            } else {
                TransorfMat = rotating < MT > (ZP, seed2);
            }
            Mat = extractMatPoly(ZP);
            break;
        }
        case 4: {
            Vpolytope VP1;
            Vpolytope VP2;
            InterVP VPcVP;
            VP1.init(n, Rcpp::as<MT>(P.field("V1")), VT::Ones(Rcpp::as<MT>(P.field("V1")).rows()));
            VP2.init(n, Rcpp::as<MT>(P.field("V2")), VT::Ones(Rcpp::as<MT>(P.field("V2")).rows()));
            VPcVP.init(VP1, VP2);
            if (T.isNotNull()) {
                TransorfMat = Rcpp::as<MT>(T);
                VPcVP.linear_transformIt(TransorfMat.inverse());
            } else {
                TransorfMat = rotating < MT > (VPcVP, seed2);
            }
        }
    }



    TransorfMat.conservativeResize(n+1, n);
    TransorfMat.row(n) = VT::Ones(n);
    MT res(TransorfMat.rows(), Rcpp::as<MT>(Mat).rows()+n);
    res << Rcpp::as<MT>(Mat).transpose(), TransorfMat;
    return Rcpp::wrap(res);

}
