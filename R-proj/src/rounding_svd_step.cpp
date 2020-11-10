// [[Rcpp::depends(BH)]]

// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

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
#include "preprocess/min_sampling_covering_ellipsoid_rounding.hpp"
#include "preprocess/svd_rounding.hpp"
#include "preprocess/max_inscribed_ellipsoid_rounding.hpp"
#include "extractMatPoly.h"


//' Internal rcpp function for the rounding of a convex polytope
//'
//' @param P A convex polytope (H- or V-representation or zonotope).
//' @param method Optional. The method to use for rounding, a) \code{'min_ellipsoid'} for the method based on mimimmum volume enclosing ellipsoid of a uniform sample from P, b) \code{'max_ellipsoid'} for the method based on maximum volume enclosed ellipsoid in P, (c) \code{'svd'} for the method based on svd decomposition. The default method is \code{'min_ellipsoid'} for all the representations.
//' @param seed Optional. A fixed seed for the number generator.
//'
//' @keywords internal
//'
//' @return A numerical matrix that describes the rounded polytope, a numerical matrix of the inverse linear transofmation that is applied on the input polytope, the numerical vector the the input polytope is shifted and the determinant of the matrix of the linear transformation that is applied on the input polytope.

//MT T;
 //   VT T_shift;
//    unsigned int num_rounding_steps;                                
//    bool fail;
//    bool converged;
//    bool last_round_under_p;
//    NT max_s;
//    NT prev_max_s;
//    unsigned int round_it;

// [[Rcpp::export]]
Rcpp::List rounding_svd_step (Rcpp::NumericMatrix Ar, Rcpp::NumericVector br,
                              Rcpp::NumericVector center, double radius, int walk_length,
                              Rcpp::List parameters){

    typedef double NT;
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT> RNGType;
    typedef HPolytope<Point> Hpolytope;
    typedef VPolytope<Point> Vpolytope;
    typedef Zonotope<Point> zonotope;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;

    MT T = Rcpp::as<MT>(parameters["T"]), A = Rcpp::as<MT>(Ar);
    VT b = Rcpp::as<VT>(br), T_shift = Rcpp::as<VT>(parameters["T_shift"]);

    int round_it = parameters["round_it"], num_rounding_steps = parameters["num_rounding_steps"];
    NT max_s = parameters["max_s"], prev_max_s = parameters["prev_max_s"];
    bool fail = parameters["fail"], converged = parameters["converged"],
         last_round_under_p = parameters["last_round_under_p"];
 
    std::pair<Point, NT> InnerBall;
    InnerBall.first = Point(Rcpp::as<VT>(center));
    InnerBall.second = radius;

    Hpolytope P(A.cols(), A, b);
    P.set_InnerBall(InnerBall);

    int n = P.dimension(); 

    //typedef BoostRandomNumberGenerator<boost::mt19937, NT> RNGType;
    RNGType rng(n);

    MT round_mat, r_inv;
    VT shift(n), s(n);

    Point p(n);

    NT s_cutoff, p_cutoff;
    MT V(n,n), S(n,n);

    round_it = 1;
    max_s = std::numeric_limits<NT>::max();
    s_cutoff = 2.3;
    p_cutoff = 10.0;
    int num_its = 20;

    p = InnerBall.first;
    svd_on_sample<AcceleratedBilliardWalk>(P, p, num_rounding_steps, V, s,
                                  shift, walk_length, rng);

    max_s = s.maxCoeff();

    if (max_s <= p_cutoff && max_s > s_cutoff) {
        if (last_round_under_p) {
            num_rounding_steps = num_rounding_steps * 2;
            p = InnerBall.first;
            svd_on_sample<AcceleratedBilliardWalk>(P, p, num_rounding_steps, V, s,
                                          shift, walk_length, rng);
            max_s = s.maxCoeff();
        } else {
            last_round_under_p = true;
        }
    } else {
        last_round_under_p = false;
    }
    S = s.asDiagonal();
    round_mat = V * S;
    r_inv = VT::Ones(n).cwiseProduct(s.cwiseInverse()).asDiagonal() * V.transpose();

    if (round_it != 1 && max_s >= NT(4) * prev_max_s) {
        fail = true;
    }

    round_it++;
    prev_max_s = max_s;

    P.shift(shift);
    P.linear_transformIt(round_mat);
    P.normalize();
    T_shift += T * shift;
    T = T * round_mat;

    if (max_s <= s_cutoff || round_it > num_its) {
        converged = true;
    }

    return Rcpp::List::create(Rcpp::Named("A") = Rcpp::wrap(A),
                              Rcpp::Named("b") = Rcpp::wrap(b),
                              Rcpp::Named("T") = Rcpp::wrap(T),
                              Rcpp::Named("T_shift") = Rcpp::wrap(T_shift),
                              Rcpp::Named("round_it") = round_it,
                              Rcpp::Named("num_rounding_steps") = num_rounding_steps,
                              Rcpp::Named("max_s") = max_s,
                              Rcpp::Named("prev_max_s") = prev_max_s,
                              Rcpp::Named("fail") = fail,
                              Rcpp::Named("converged") = converged,
                              Rcpp::Named("last_round_under_p") = last_round_under_p);

}

