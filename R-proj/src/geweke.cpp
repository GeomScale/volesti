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
#include "diagnostics/geweke.hpp"


//' Internal rcpp function for the rounding of a convex polytope
//'
//' @param samples The point set.
//' @param frac1 Optional.
//' @param frac2 Optional. 
//'
//' @return A numerical matrix that describes the rounded polytope, a numerical matrix of the inverse linear transofmation that is applied on the input polytope, the numerical vector the the input polytope is shifted and the determinant of the matrix of the linear transformation that is applied on the input polytope.
// [[Rcpp::export]]
bool geweke(Rcpp::NumericMatrix samples, Rcpp::Nullable<double> frac1 = R_NilValue, 
              Rcpp::Nullable<double> frac2 = R_NilValue)
{
    typedef double NT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;

    NT frac_1, frac_2;

    frac_1 = (!frac1.isNotNull()) ? NT(0.1) : Rcpp::as<NT>(frac1);
    frac_2 = (!frac1.isNotNull()) ? NT(0.5) : Rcpp::as<NT>(frac2);

    return perform_geweke<VT>(Rcpp::as<MT>(samples), frac_1, frac_2);

}
