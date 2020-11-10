// [[Rcpp::depends(BH)]]
// [[ Rcpp :: depends ( RcppArmadillo )]]

// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2020 Vissarion Fisikopoulos
// Copyright (c) 2020 Apostolos Chalkis

//Contributed and/or modified by Alexandros Manochis, as part of Google Summer of Code 2020 program.

//#include <Rcpp.h>
#include <RcppArmadillo.h>
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

//' Internal rcpp function to compute the full dimensional polytope when a low dimensional is given
//'
//' @param P A low dimensional convex polytope in H-representation.
//'
//' @keywords internal
//'
//' @return A numerical matrix that describes the full dimensional polytope, a numerical matrix of the inverse
//'         linear transformation that is applied on the input polytope, the numerical vector - point that the
//'         input polytope is shifted and the product of the singular values of the matrix of the linear map 
//'         applied on the input polytope.
//'
// [[Rcpp::export]]
Rcpp::List full_dimensional_polytope (Rcpp::Reference P)
{
    typedef double NT;
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef HPolytope<Point> Hpolytope;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;

    NT svd_prod;
    MT A = Rcpp::as<MT>(P.field("A")), Aeq = Rcpp::as<MT>(P.field("Aeq"));
    VT b = Rcpp::as<VT>(P.field("b")), beq = Rcpp::as<VT>(P.field("beq"));

    std::pair<Hpolytope, std::pair<MT, VT> > result = get_full_dimensional_polytope<Hpolytope>(A, b, Aeq, beq);

    Hpolytope HP = result.first;
    MT N = result.second.first;
    Rcpp::NumericMatrix Mat = extractMatPoly(HP);

    Eigen::JacobiSVD<MT> svd(N, Eigen::ComputeThinU | Eigen::ComputeThinV);
    svd_prod = svd.singularValues().prod();

    return Rcpp::List::create(Rcpp::Named("Mat") = Mat, Rcpp::Named("N") = Rcpp::wrap(N),
                              Rcpp::Named("shift") = Rcpp::wrap(result.second.second),
                              Rcpp::Named("svd_prod") = svd_prod);
}


// [[Rcpp::export]]
Rcpp::List full_dimensional_polytope_with_arma (Rcpp::NumericMatrix Ar, Rcpp::NumericVector br)
{
    typedef double NT;
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef HPolytope<Point> Hpolytope;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;

    MT Aeq = Rcpp::as<MT>(Ar);
    VT beq = Rcpp::as<VT>(br);

    VT p = Aeq.fullPivLu().solve(beq);
    //b = b - A * p;

    arma::mat Arm = Rcpp::as<arma::mat>(Ar);
    arma::mat N = arma::null(Arm);

    return Rcpp::List::create(Rcpp::Named("N") = Rcpp::wrap(N),
                              Rcpp::Named("N_shift") = Rcpp::wrap(p));
}

/*Rcpp::List full_dimensional_polytope_with_arma_sparse (Rcpp::NumericMatrix Ar, Rcpp::NumericVector br)
{
    typedef double NT;
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef HPolytope<Point> Hpolytope;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;

    MT Aeq = Rcpp::as<MT>(Ar);
    VT beq = Rcpp::as<VT>(br);

    VT p = Aeq.colPivHouseholderQr().solve(beq);
    //b = b - A * p;

    arma::mat Arm = Rcpp::as<arma::mat>(Ar);
    arma::SpMat<double> Arm_sp(Arm);
    arma::SpMat<double> NN = arma::null(Arm_sp);
    arma::mat N(NN);

    return Rcpp::List::create(Rcpp::Named("N") = Rcpp::wrap(N),
                              Rcpp::Named("N_shift") = Rcpp::wrap(p));
}*/
