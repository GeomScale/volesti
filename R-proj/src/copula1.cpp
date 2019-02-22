// [[Rcpp::depends(BH)]]

// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

#include <Rcpp.h>
#include <chrono>
#include "cartesian_geom/cartesian_kernel.h"
#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include "simplex_samplers.h"
#include "copulas.h"

//' Construct a copula using uniform sampling from the unit simplex
//'
//' Given two families of parallel hyperplanes intersecting the canonical simplex, this function uniformly samples from the canonical simplex and construct an approximation of the bivariate probability distribution, called copula.
//'
//' @param h1 A \eqn{d}-dimensional vector that describes the direction of the first family of parallel hyperplanes.
//' @param h2 A \eqn{d}-dimensional vector that describes the direction of the second family of parallel hyperplanes.
//' @param numSlices The number of the slices for the copula. Default value is 100.
//' @param N The number of points to sample. Default value is \eqn{4\cdot 10^6}.
//'
//' @references \cite{L. Cales, A. Chalkis, I.Z. Emiris, V. Fisikopoulos,
//' \dQuote{Practical volume computation of structured convex bodies, and an application to modeling portfolio dependencies and financial crises,} \emph{Proc. of Symposium on Computational Geometry, Budapest, Hungary,} 2018.}
//'
//' @return A \eqn{numSlices\times numSlices} numerical matrix that corresponds to a copula.
//' @examples
//' # compute a copula for two random families of parallel hyperplanes
//' h1 = runif(n = 10, min = 1, max = 1000)
//' h1 = h1 / 1000
//' h2=runif(n = 10, min = 1, max = 1000)
//' h2 = h2 / 1000
//' cop = copula1(h1=h1, h2=h2, numSlices = 10, N = 100000)
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix copula1 (Rcpp::NumericVector h1, Rcpp::NumericVector h2,
                             Rcpp::Nullable<unsigned int> numSlices, Rcpp::Nullable<unsigned int> N){

    typedef double NT;
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef boost::mt19937    RNGType;
    unsigned int num_slices = 100, numpoints = 4000000;

    if (numSlices.isNotNull()) {
        num_slices = Rcpp::as<unsigned int>(numSlices);
    }

    if (N.isNotNull()) {
        numpoints = Rcpp::as<unsigned int>(N);
    }

    Rcpp::NumericMatrix copula(num_slices, num_slices);
    std::vector<std::vector<NT> > StdCopula;
    unsigned int dim = h1.size(), i, j;

    std::vector<NT> hyp1 = Rcpp::as<std::vector<NT> >(h1);
    std::vector<NT> hyp2 = Rcpp::as<std::vector<NT> >(h2);

    StdCopula = twoParHypFam<Point, RNGType >(dim, numpoints, num_slices, hyp1, hyp2);

    for(i=0; i<num_slices; i++) {
        for(j=0; j<num_slices; j++){
            copula(i,j) = StdCopula[i][j];
        }
    }

    return copula;
}
