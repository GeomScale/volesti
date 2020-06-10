// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis


#include <Rcpp.h>
#include <RcppEigen.h>
#include "volume/exact_vols.h"

//' Compute the percentage of the volume of the simplex that is contained in the intersection of a half-space and the simplex.
//'
//' A half-space \eqn{H} is given as a pair of a vector \eqn{a\in R^d} and a scalar \eqn{z0\in R} s.t.: \eqn{a^Tx\leq z0}. This function calls the Ali's version of the Varsi formula to compute a frustum of the simplex.
//'
//' @param a A \eqn{d}-dimensional vector that defines the direction of the hyperplane.
//' @param z0 The scalar that defines the half-space.
//'
//' @references \cite{Varsi, Giulio,
//' \dQuote{The multidimensional content of the frustum of the simplex,} \emph{Pacific J. Math. 46, no. 1, 303--314,} 1973.}
//'
//' @references \cite{Ali, Mir M.,
//' \dQuote{Content of the frustum of a simplex,} \emph{ Pacific J. Math. 48, no. 2, 313--322,} 1973.}
//'
//' @return The percentage of the volume of the simplex that is contained in the intersection of a given half-space and the simplex.
//'
//' @examples
//' # compute the frustum of H: -x1+x2<=0
//' a=c(-1,1)
//' z0=0
//' frustum = frustum_of_simplex(a, z0)
//' @export
// [[Rcpp::export]]
double frustum_of_simplex(Rcpp::NumericVector a, double z0){

    unsigned int dim = a.size();
    if (dim < 2) {
        throw Rcpp::exception("Dimension has to be greater than 2");
    }
    std::vector<double> hyp = Rcpp::as<std::vector<double> >(a);

    return vol_Ali(hyp, -z0, dim);

}
