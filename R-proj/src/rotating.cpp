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
#include "volume/volume_sequence_of_balls.hpp"
#include "volume/rotating.hpp"
#include "extractMatPoly.h"

//'  An internal Rccp function for the random rotation of a convex polytope
//'
//' @param P A convex polytope (H-, V-polytope or a zonotope).
//' @param T Optional. A rotation matrix.
//' @param seed Optional. A fixed seed for the random linear map generator.
//'
//' @keywords internal
//'
//' @return A matrix that describes the rotated polytope
// [[Rcpp::export]]
Rcpp::NumericMatrix rotating (Rcpp::Reference P, Rcpp::Nullable<Rcpp::NumericMatrix> T = R_NilValue,
                              Rcpp::Nullable<int> seed = R_NilValue){

    typedef double NT;
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT> RNGType;
    typedef HPolytope<Point> Hpolytope;
    typedef VPolytope<Point> Vpolytope;
    typedef Zonotope<Point> zonotope;
    typedef IntersectionOfVpoly<Vpolytope, RNGType> InterVP;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;

    MT TransorfMat;
    Rcpp::NumericMatrix Mat;
    unsigned int n, type_num;

    std::string type = Rcpp::as<std::string>(P.slot("type"));

    if (type.compare(std::string("Hpolytope")) == 0) {
        n = Rcpp::as<MT>(P.slot("A")).cols();
        type_num = 1;
    } else if (type.compare(std::string("Vpolytope")) == 0) {
        n = Rcpp::as<MT>(P.slot("V")).cols();
        type_num = 2;
    } else if (type.compare(std::string("Zonotope")) == 0) {
        n = Rcpp::as<MT>(P.slot("G")).cols();
        type_num = 3;
    } else if (type.compare(std::string("VpolytopeIntersection")) == 0) {
        throw Rcpp::exception("volesti does not support roatation of this kind of representation.");
    } else {
        throw Rcpp::exception("Unknown polytope representation!");
    }

    int seed2 = (!seed.isNotNull()) ? std::chrono::system_clock::now()
                                      .time_since_epoch().count()
                                    : Rcpp::as<int>(seed);

    switch (type_num) {
        case 1: {
            // Hpolytope
            Hpolytope HP;
            HP.init(n, Rcpp::as<MT>(P.slot("A")), Rcpp::as<VT>(P.slot("b")));
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
            VP.init(n, Rcpp::as<MT>(P.slot("V")), VT::Ones(Rcpp::as<MT>(P.slot("V")).rows()));
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
            ZP.init(n, Rcpp::as<MT>(P.slot("G")), VT::Ones(Rcpp::as<MT>(P.slot("G")).rows()));
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
            throw Rcpp::exception("volesti does not support rotation for this representation currently.");
        }
    }

    TransorfMat.conservativeResize(n+1, n);
    TransorfMat.row(n) = VT::Ones(n);
    MT res(TransorfMat.rows(), Rcpp::as<MT>(Mat).rows()+n);
    res << Rcpp::as<MT>(Mat).transpose(), TransorfMat;
    return Rcpp::wrap(res);
}
