// [[Rcpp::depends(BH)]]

// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2020 Vissarion Fisikopoulos
// Copyright (c) 2020 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.
//Contributed and/or modified by Alexandros Manochis, as part of Google Summer of Code 2020 program.

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
#include "preprocess/get_full_dimensional_polytope.hpp"
#include "extractMatPoly.h"

//' Internal rcpp function for the rounding of a convex polytope
//'
//' @param P A convex polytope (H- or V-representation or zonotope).
//' @param method Optional. The method to use for rounding, a) \code{'mee'} for the method based on mimimmum volume enclosing ellipsoid of a dataset, b) \code{'mve'} for the method based on maximum volume enclosed ellipsoid, (c) \code{'svd'} for the method based on svd decomposition.  The default method is \code{'mee'} for all the representations.
//' @param seed Optional. A fixed seed for the number generator.
//'
//' @keywords internal
//'
//' @return A numerical matrix that describes the rounded polytope, a numerical matrix of the inverse linear transofmation that is applied on the input polytope, the numerical vector the the input polytope is shifted and the determinant of the matrix of the linear transformation that is applied on the input polytope.
// [[Rcpp::export]]
Rcpp::List rounding (Rcpp::Reference P, Rcpp::Nullable<std::string> method = R_NilValue,
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

    bool cdhr = false;
    unsigned int n = P.field("dimension"), walkL, type = P.field("type");
    std::string mthd = std::string("mee");
    if(method.isNotNull()) {
        mthd =  Rcpp::as<std::string>(method);
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

    if (type == 1) {
        walkL = 10 + 10*n;
        cdhr = true;
    } else {
        walkL = 2;
    }

    MT N = MT::Identity(n,n);
    VT N_shift = VT::Zero(n);
    NT svd_prod = 1.0;
    switch (type) {
        case 1: {
            // Hpolytope
            if (Rcpp::as<MT>(P.field("Aeq")).rows() == 0) {
                // full dimensional
                HP.init(n, Rcpp::as<MT>(P.field("A")), Rcpp::as<VT>(P.field("b")));
            } else {
                // low dimensional
                std::pair<Hpolytope, std::pair<MT, VT> > temp_res = get_full_dimensional_polytope<Hpolytope>(Rcpp::as<MT>(P.field("A")), Rcpp::as<VT>(P.field("b")),
                            Rcpp::as<MT>(P.field("Aeq")), Rcpp::as<VT>(P.field("beq")));

                HP = temp_res.first;
                N = temp_res.second.first;
                N_shift = temp_res.second.second;

                Eigen::JacobiSVD<MT> svd(N, Eigen::ComputeThinU | Eigen::ComputeThinV);
                svd_prod = svd.singularValues().prod();
                walkL = 10 + 10*HP.dimension();
              
                rng = RNGType(HP.dimension());
                if (seed.isNotNull()) {
                    unsigned seed2 = Rcpp::as<double>(seed);
                    rng.set_seed(seed2);
                }
            }
            HP.normalize();
            InnerBall = HP.ComputeInnerBall();
            break;
        }
        case 2: {
            // Vpolytope
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
        case 4: {
            throw Rcpp::exception("volesti does not support rounding for this representation currently.");
        }
    }

    std::pair< std::pair<MT, VT>, NT > round_res;
    switch (type) {
        case 1: {
            if (mthd.compare(std::string("mve")) == 0) {
                round_res = mve_rounding<MT, VT>(HP, InnerBall);
            } else if (mthd.compare(std::string("mee")) == 0) {
                if (cdhr) {
                    round_res = round_polytope<CDHRWalk, MT, VT>(HP, InnerBall, walkL, rng);
                } else {
                    round_res = round_polytope<BilliardWalk, MT, VT>(HP, InnerBall, walkL, rng);
                }
            } else if (mthd.compare(std::string("svd")) == 0) {
                if (cdhr) {
                    round_res = round_isotropy<CDHRWalk, MT, VT>(HP, InnerBall, walkL, rng);
                } else {
                    round_res = round_isotropy<BilliardWalk, MT, VT>(HP, InnerBall, walkL, rng);
                }
            } else {
                throw Rcpp::exception("Unknown method!");
            }
            Mat = extractMatPoly(HP);
            break;
        }
        case 2: {
            if (mthd.compare(std::string("mee")) == 0){
                if (cdhr) {
                    round_res = round_polytope<CDHRWalk, MT, VT>(VP, InnerBall, walkL, rng);
                } else {
                    round_res = round_polytope<BilliardWalk, MT, VT>(VP, InnerBall, walkL, rng);
                }
            } else if (mthd.compare(std::string("svd")) == 0) {
                if (cdhr) {
                    round_res = round_isotropy<CDHRWalk, MT, VT>(VP, InnerBall, walkL, rng);
                } else {
                    round_res = round_isotropy<BilliardWalk, MT, VT>(VP, InnerBall, walkL, rng);
                }
            } else {
                throw Rcpp::exception("Unknown method!");
            }
            Mat = extractMatPoly(VP);
            break;
        }
        case 3: {
            if (mthd.compare(std::string("mee")) == 0){
                if (cdhr) {
                    round_res = round_polytope<CDHRWalk, MT, VT>(ZP, InnerBall, walkL, rng);
                } else {
                    round_res = round_polytope<BilliardWalk, MT, VT>(ZP, InnerBall, walkL, rng);
                }
            } else if (mthd.compare(std::string("svd")) == 0) {
                if (cdhr) {
                    round_res = round_isotropy<CDHRWalk, MT, VT>(ZP, InnerBall, walkL, rng);
                } else {
                    round_res = round_isotropy<BilliardWalk, MT, VT>(ZP, InnerBall, walkL, rng);
                }
            } else {
                throw Rcpp::exception("Unknown method!");
            }
            Mat = extractMatPoly(ZP);
            break;
        }
    }

    return Rcpp::List::create(Rcpp::Named("Mat") = Mat, Rcpp::Named("T") = Rcpp::wrap(round_res.first.first),
                              Rcpp::Named("shift") = Rcpp::wrap(round_res.first.second),
                              Rcpp::Named("round_value") = round_res.second,
                              Rcpp::Named("N") = Rcpp::wrap(N),
                              Rcpp::Named("N_shift") = Rcpp::wrap(N_shift),
                              Rcpp::Named("svd_prod") = svd_prod);
}
