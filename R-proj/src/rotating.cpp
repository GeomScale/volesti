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
#include "polytopes.h"
#include "samplers.h"
#include "rotating.h"
#include "extractMatPoly.h"

//'  An internal Rccp function for the random rotation of a convex polytope
//'
//' @param P A convex polytope (H-, V-polytope or a zonotope).
//'
//' @section warning:
//' Do not use this function.
//'
//' @return A matrix that describes the rotated polytope
// [[Rcpp::export]]
Rcpp::NumericMatrix rotating (Rcpp::Reference P){

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

    MT TransorfMat;
    Rcpp::NumericMatrix Mat;
    unsigned int n = P.field("dimension");
    int type = P.field("type");

    switch (type) {
        case 1: {
            // Hpolytope
            Hpolytope HP;
            HP.init(n, Rcpp::as<MT>(P.field("A")), Rcpp::as<VT>(P.field("b")));
            TransorfMat = rotating < MT > (HP);
            Mat = extractMatPoly(HP);
            break;
        }
        case 2: {
            // Vpolytope
            Vpolytope VP;
            VP.init(n, Rcpp::as<MT>(P.field("V")), VT::Ones(Rcpp::as<MT>(P.field("V")).rows()));
            TransorfMat = rotating< MT >(VP);
            Mat = extractMatPoly(VP);
            break;
        }
        case 3: {
            // Zonotope
            zonotope ZP;
            ZP.init(n, Rcpp::as<MT>(P.field("G")), VT::Ones(Rcpp::as<MT>(P.field("G")).rows()));
            TransorfMat = rotating < MT > (ZP);
            Mat = extractMatPoly(ZP);
            break;
        }
    }

    TransorfMat.conservativeResize(n+1, n);
    TransorfMat.row(n) = VT::Ones(n);
    MT res(TransorfMat.rows(), Rcpp::as<MT>(Mat).rows()+n);
    res << Rcpp::as<MT>(Mat).transpose(), TransorfMat;
    return Rcpp::wrap(res);

}
