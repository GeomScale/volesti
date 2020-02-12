// [[Rcpp::depends(BH)]]

// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

#include <Rcpp.h>
#include <RcppEigen.h>
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
//' @param M The number of the slices for the copula. The default value is 100.
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
//' cop = copula(R1 = h1, R2 = h2, M = 10, N = 100000)
//'
//' # compute a copula for a family of parallel hyperplanes and a family of conentric ellipsoids
//' h = runif(n = 10, min = 1, max = 1000)
//' h = h / 1000
//' E = replicate(10, rnorm(20))
//' E = cov(E)
//' cop = copula(R1 = h, Sigma = E, M = 10, N = 100000)
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix copula (Rcpp::Nullable<Rcpp::NumericVector> R1 = R_NilValue,
                            Rcpp::Nullable<Rcpp::NumericVector> R2 = R_NilValue,
                            Rcpp::Nullable<Rcpp::NumericMatrix> Sigma = R_NilValue,
                            Rcpp::Nullable<unsigned int> M = R_NilValue,
                            Rcpp::Nullable<unsigned int> N = R_NilValue){

    typedef double NT;
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef boost::mt19937    RNGType;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    typedef copula_ellipsoid<Point, MT, VT> CopEll;
    unsigned int num_slices = 100, numpoints = 500000;

    if (M.isNotNull()) {
        num_slices = Rcpp::as<unsigned int>(M);
    }

    if (N.isNotNull()) {
        numpoints = Rcpp::as<unsigned int>(N);
    }

    Rcpp::NumericMatrix copula(num_slices, num_slices);
    std::vector<std::vector<NT> > StdCopula;
    unsigned int dim = Rcpp::as<std::vector<NT> >(R1).size(), i, j;

    if(!R1.isNotNull()) {

        throw Rcpp::exception("You have to give at least one normal of a hyperplane!");

    }

    std::vector<NT> hyp1 = Rcpp::as<std::vector<NT> >(R1);
    if (R2.isNotNull()) {

        std::vector <NT> hyp2 = Rcpp::as < std::vector < NT > > (R2);
        StdCopula = twoParHypFam<Point, RNGType>(dim, numpoints, num_slices, hyp1, hyp2);

    } else if (Sigma.isNotNull()) {

        std::vector<std::vector<NT> > Gin(dim, std::vector<NT>(dim));
        MT EE = Rcpp::as<MT>(Sigma);
        for (int i=0; i<dim; i++) {
            for (int j = 0; j < dim; j++) {
                Gin[i][j] = EE(i, j);
            }
        }
        CopEll Ell(Gin);
        StdCopula = hypfam_ellfam<Point, RNGType >(dim, numpoints, num_slices, hyp1, Ell);
    } else {

        throw Rcpp::exception("Wrong inputs");

    }

    for(int i=0; i<num_slices; i++) {
        for(int j=0; j<num_slices; j++){
            copula(i,j) = StdCopula[i][j];
        }
    }

    return copula;
}
