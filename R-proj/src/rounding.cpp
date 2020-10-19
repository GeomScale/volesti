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
#include "random_walks/random_walks.hpp"
#include "volume/volume_sequence_of_balls.hpp"
#include "volume/volume_cooling_gaussians.hpp"
#include "extractMatPoly.h"

//' Internal rcpp function for the rounding of a convex polytope
//'
//' @param P A convex polytope (H- or V-representation or zonotope).
//' @param settings A list to set the random walk and its walk length
//' @param seed Optional. A fixed seed for the number generator.
//'
//' @keywords internal
//'
//' @return A numerical matrix that describes the rounded polytope, a numerical matrix of the inverse linear transofmation that is applied on the input polytope, the numerical vector the the input polytope is shifted and the determinant of the matrix of the linear transformation that is applied on the input polytope.
// [[Rcpp::export]]
Rcpp::List rounding (Rcpp::Reference P,
                     Rcpp::Nullable<Rcpp::List> settings = R_NilValue,
                     Rcpp::Nullable<double> seed = R_NilValue){

    typedef double NT;
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT> RNGType;
    typedef HPolytope<Point> Hpolytope;
    typedef VPolytope<Point> Vpolytope;
    typedef Zonotope<Point> zonotope;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;

    bool cdhr = false, rdhr = false;
    unsigned int n, walkL, type_num;

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

    RNGType rng(n);
    if (seed.isNotNull()) {
        unsigned seed2 = Rcpp::as<double>(seed);
        rng.set_seed(seed2);
    }

    Hpolytope HP;
    Vpolytope VP;
    zonotope ZP;

    std::pair <Point, NT> InnerBall;
    Rcpp::NumericMatrix Mat;

    if (!settings.isNotNull() || !Rcpp::as<Rcpp::List>(settings).containsElementNamed("random_walk")) {
        if (type_num == 1) {
            walkL = 10 + 10 * n;
            cdhr = true;
        } else {
            walkL = 2;
        }
    } else if (Rcpp::as<std::string>(Rcpp::as<Rcpp::List>(settings)["random_walk"]).compare(std::string("CDHR")) == 0) {
        walkL = 10 + 10 * n;
        cdhr = true;
    } else if (Rcpp::as<std::string>(Rcpp::as<Rcpp::List>(settings)["random_walk"]).compare(std::string("BiW")) == 0) {
        walkL = 2;
    } else if (Rcpp::as<std::string>(Rcpp::as<Rcpp::List>(settings)["random_walk"]).compare(std::string("RDHR")) == 0) {
        walkL = 10 + 10 * n;
        rdhr = true;
    } else {
        throw Rcpp::exception("Unknown walk type!");
    }

    if (Rcpp::as<Rcpp::List>(settings).containsElementNamed("walk_length")) {
        walkL = Rcpp::as<int>(Rcpp::as<Rcpp::List>(settings)["walk_length"]);
        if (walkL <= 0) {
            throw Rcpp::exception("The walk length has to be a positive integer!");
        }
    }

    switch (type_num) {
        case 1: {
            // Hpolytope
            HP.init(n, Rcpp::as<MT>(P.slot("A")), Rcpp::as<VT>(P.slot("b")));
            HP.normalize();
            InnerBall = HP.ComputeInnerBall();
            break;
        }
        case 2: {
            // Vpolytope
            VP.init(n, Rcpp::as<MT>(P.slot("V")), VT::Ones(Rcpp::as<MT>(P.slot("V")).rows()));
            InnerBall = VP.ComputeInnerBall();
            break;
        }
        case 3: {
            // Zonotope
            ZP.init(n, Rcpp::as<MT>(P.slot("G")), VT::Ones(Rcpp::as<MT>(P.slot("G")).rows()));
            InnerBall = ZP.ComputeInnerBall();
            break;
        }
        case 4: {
            throw Rcpp::exception("volesti does not support rounding for this representation currently.");
        }
    }

    std::pair< std::pair<MT, VT>, NT > round_res;
    switch (type_num) {
        case 1: {
            if (cdhr) {
                round_res = round_polytope<CDHRWalk, MT, VT>(HP, InnerBall, walkL, rng);
            } else if(rdhr) {
                round_res = round_polytope<RDHRWalk, MT, VT>(HP, InnerBall, walkL, rng);
            } else {
                round_res = round_polytope<BilliardWalk, MT, VT>(HP, InnerBall, walkL, rng);
            }
            Mat = extractMatPoly(HP);
            break;
        }
        case 2: {
            if (cdhr) {
                round_res = round_polytope<CDHRWalk, MT, VT>(VP, InnerBall, walkL, rng);
            } else if(rdhr) {
                round_res = round_polytope<RDHRWalk, MT, VT>(VP, InnerBall, walkL, rng);
            } else {
                round_res = round_polytope<BilliardWalk, MT, VT>(VP, InnerBall, walkL, rng);
            }
            Mat = extractMatPoly(VP);
            break;
        }
        case 3: {
            if (cdhr) {
                round_res = round_polytope<CDHRWalk, MT, VT>(ZP, InnerBall, walkL, rng);
            } else if (rdhr) {
                round_res = round_polytope<RDHRWalk, MT, VT>(ZP, InnerBall, walkL, rng);
            } else {
                round_res = round_polytope<BilliardWalk, MT, VT>(ZP, InnerBall, walkL, rng);
            }
            Mat = extractMatPoly(ZP);
            break;
        }
    }

    return Rcpp::List::create(Rcpp::Named("Mat") = Mat, Rcpp::Named("T") = Rcpp::wrap(round_res.first.first),
                              Rcpp::Named("shift") = Rcpp::wrap(round_res.first.second),
                              Rcpp::Named("round_value") = round_res.second);
}
