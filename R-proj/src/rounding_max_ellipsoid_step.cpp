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
#include "preprocess/max_inscribed_ellipsoid_rounding.hpp"


//' Internal rcpp function for the rounding of a convex polytope
//'
//' @param P A convex polytope (H- or V-representation or zonotope).
//' @param method Optional. The method to use for rounding, a) \code{'min_ellipsoid'} for the method based on mimimmum volume enclosing ellipsoid of a uniform sample from P, b) \code{'max_ellipsoid'} for the method based on maximum volume enclosed ellipsoid in P, (c) \code{'svd'} for the method based on svd decomposition. The default method is \code{'min_ellipsoid'} for all the representations.
//' @param seed Optional. A fixed seed for the number generator.
//'
//' @keywords internal
//'
//' @return A numerical matrix that describes the rounded polytope, a numerical matrix of the inverse linear transofmation that is applied on the input polytope, the numerical vector the the input polytope is shifted and the determinant of the matrix of the linear transformation that is applied on the input polytope.
// [[Rcpp::export]]
Rcpp::List rounding_max_ellipsoid_step(Rcpp::NumericMatrix Ar, Rcpp::NumericVector br, 
                                        Rcpp::NumericVector center, double radius){

    typedef double NT;
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT> RNGType;
    typedef HPolytope<Point> Hpolytope;
    typedef VPolytope<Point> Vpolytope;
    typedef Zonotope<Point> zonotope;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;

    MT A = Rcpp::as<MT>(Ar);
    MT b = Rcpp::as<VT>(br);
    std::pair <Point, NT> InnerBall;

    InnerBall.first = Point(Rcpp::as<VT>(center));
    InnerBall.second = radius;

    Hpolytope HP(A.cols(), A, b);
    HP.set_InnerBall(InnerBall);

    unsigned int n = HP.dimension();
    RNGType rng(n);

    std::tuple<MT, VT, NT> round_res = max_inscribed_ellipsoid_rounding<MT, VT, NT>(HP, InnerBall.first, 1);
            
    MT A2 = HP.get_mat();
    VT b2 = HP.get_vec();

    return Rcpp::List::create(Rcpp::Named("A") = A2, Rcpp::Named("b") = b2, Rcpp::Named("T") = Rcpp::wrap(std::get<0>(round_res)),
                              Rcpp::Named("shift") = Rcpp::wrap(std::get<1>(round_res)),
                              Rcpp::Named("round_value") = std::get<2>(round_res));

}

