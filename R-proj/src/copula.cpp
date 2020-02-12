// [[Rcpp::depends(BH)]]

// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

#include <Rcpp.h>
#include <chrono>
#include "ellipsoids.h"
#include "cartesian_geom/cartesian_kernel.h"
#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include "simplex_samplers.h"
#include "copulas.h"

//' Construct a copula using uniform sampling from the unit simplex
//'
//' Given two families of parallel hyperplanes or a family of parallel hyperplanes and a family of concentric ellispoids centered at the origin intersecting the canonical simplex, this function uniformly samples from the canonical simplex and construct an approximation of the bivariate probability distribution, called copula.
//'
//' @param h1 A \eqn{d}-dimensional vector that describes the direction of the first family of parallel hyperplanes.
//' @param h2 Optional. A \eqn{d}-dimensional vector that describes the direction of the second family of parallel hyperplanes.
//' @param E Optional. The \eqn{d\times d} symmetric positive semidefine matrix that describes the family of concentric ellipsoids centered at the origin.
//' @param numSlices The number of the slices for the copula. The default value is 100.
//' @param N The number of points to sample. The default value is \eqn{5\cdot 10^5}.
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
//' cop = copula(h1=h1, h2=h2, numSlices = 10, N = 100000)
//'
//' # compute a copula for a family of parallel hyperplanes and a family of conentric ellipsoids
//' h = runif(n = 10, min = 1, max = 1000)
//' h = h / 1000
//' E = replicate(10, rnorm(20))
//' E = cov(E)
//' cop = copula(h1=h, E=E, numSlices=10, N=100000)
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix copula (Rcpp::NumericVector h1 = R_NilValue, Rcpp::NumericVector h2 = R_NilValue,
        Rcpp::NumericMatrix E = R_NilValue, Rcpp::Nullable<unsigned int> numSlices = R_NilValue,
        Rcpp::Nullable<unsigned int> N = R_NilValue){

    typedef double NT;
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef boost::mt19937    RNGType;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    unsigned int num_slices = 100, numpoints = 500000;

    if (numSlices.isNotNull()) {
        num_slices = Rcpp::as<unsigned int>(numSlices);
    }

    if (N.isNotNull()) {
        numpoints = Rcpp::as<unsigned int>(N);
    }

    Rcpp::NumericMatrix copula(num_slices, num_slices);
    std::vector<std::vector<NT> > StdCopula;
    unsigned int dim = Rcpp::as<std::vector<NT> >(h1).size(), i, j;

    if(!h1.isNotNull()) {

        throw Rcpp::exception("You have to give at least one normal of a hyperplane!");

    }

    std::vector<NT> hyp1 = Rcpp::as<std::vector<NT> >(h1);
    if (h2.isNotNull()) {

        std::vector <NT> hyp2 = Rcpp::as < std::vector < NT > > (h2);
        StdCopula = twoParHypFam<Point, RNGType>(dim, numpoints, num_slices, hyp1, hyp2);

    } else if (E.isNotNull()) {

        std::vector<std::vector<NT> > Gin(dim, std::vector<NT>(dim));
        MT EE = Rcpp::as<MT>(E);
        for (i=0; i<dim; i++) {
            for (j = 0; j < dim; j++) {
                Gin[i][j] = EE(i, j);
            }
        }
        CopEll Ell(Gin);
        StdCopula = hypfam_ellfam<Point, RNGType >(dim, numpoints, num_slices, hyp1, Ell);
    } else {

        throw Rcpp::exception("Wrong inputs");

    }

    for(i=0; i<num_slices; i++) {
        for(j=0; j<num_slices; j++){
            copula(i,j) = StdCopula[i][j];
        }
    }

    return copula;
}
