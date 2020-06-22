// [[Rcpp::depends(BH)]]

// VolEsti (volume computation and sampling library)

// Copyright (c) 20014-2020 Vissarion Fisikopoulos
// Copyright (c) 2020 Apostolos Chalkis

//Contributed and/or modified by Alexandros Manochis, as part of Google Summer of Code 2020 program.

#include <Rcpp.h>
#include <RcppEigen.h>
#include <chrono>
#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include "diagnostics/psrf.hpp"
#include "diagnostics/geweke.hpp"

//' Internal rcpp function for the rounding of a convex polytope
//'
//' @param samples The point set.
//' @param method Optional. The statistical method to use.
//'
//' @return A numerical matrix that describes the rounded polytope, a numerical matrix of the inverse linear transofmation that is applied on the input polytope, the numerical vector the the input polytope is shifted and the determinant of the matrix of the linear transformation that is applied on the input polytope.
// [[Rcpp::export]]
double mc_diagnostics(Rcpp::NumericMatrix samples, Rcpp::Nullable<std::string> method = R_NilValue)
{
    typedef double NT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;

    NT frac1, frac2;

    

    return perform_psrf(Rcpp::as<MT>(samples));

}
