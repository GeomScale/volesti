// [[Rcpp::depends(BH)]]

// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2019 Vissarion Fisikopoulos
// Copyright (c) 2018-2019 Apostolos Chalkis

#include <Rcpp.h>
#include <RcppEigen.h>
#include <chrono>
#include "cartesian_geom/cartesian_kernel.h"
#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include "polytopes.h"
#include "exact_vols.h"

template <typename FT>
FT factorial(FT n)
{
    return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

//' Compute the exact volume of (a) a zonotope (b) an arbitrary simplex (c) a unit simplex (d) a cross polytope (e) a hypercube
//'
//' Given a zonotope (as an object of class Zonotope), this function computes the sum of the absolute values of the determinants of all the \eqn{d \times d} submatrices of the \eqn{m\times d} matrix \eqn{G} that contains row-wise the \eqn{m} \eqn{d}-dimensional segments that define the zonotope.
//' For an arbitrary simplex that is given in V-representation this function computes the absolute value of the determinant formed by the simplex's points assuming it is shifted to the origin.
//' For a \eqn{d}-dimensional unit simplex, hypercube or cross polytope this function computes the exact well known formulas.
//'
//' @param P A zonotope or a simplex in V-representation.
//' @param body A string that declares the type of the body for the exact sampling: a) \code{'simplex'} for the unit simplex, b) \code{'cross'} for the cross polytope, c) \code{'hypersphere'} for the hypersphere, d) \code{'cube'} for the unit cube.
//' @param Parameters A list for the parameters of the methods:
//' \itemize{
//' \item{\code{dimension} }{ An integer that declares the dimension when exact sampling is enabled for a simplex or a hypersphere.}
//' \item{\code{radius} }{ The radius of the \eqn{d}-dimensional hypersphere. Default value is \eqn{1}.}
//' }
//'
//' @return The exact volume of the zonotope
//' @examples
//'
//' # compute the exact volume of a 5-dimensional zonotope defined by the Minkowski sum of 10 segments
//' Z = GenZonotope(5, 10)
//' vol = exact_vol(Z)
//'
//' \donttest{# compute the exact volume of a 2-d arbitrary simplex
//' V = matrix(c(2,3,-1,7,0,0),ncol = 2, nrow = 3, byrow = TRUE)
//' P = Vpolytope$new(V)
//' vol = exact_vol(P)
//' }
//'
//' # compute the exact volume the 10-dimensional cross polytope
//' vol = exact_vol(body = "cross", Parameters = list("dimension" = 10))
//' @export
// [[Rcpp::export]]
double exact_vol(Rcpp::Nullable<Rcpp::Reference> P = R_NilValue, Rcpp::Nullable<std::string> body = R_NilValue,
                 Rcpp::Nullable<Rcpp::List> Parameters = R_NilValue){

    typedef double NT;
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef boost::mt19937    RNGType;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;

    typedef Zonotope<Point> zonotope;
    typedef VPolytope <Point, RNGType> Vpolytope;
    zonotope ZP;
    Vpolytope VP;
    int type, dim;
    NT vol, rad;

    if (P.isNotNull()) {
        type = Rcpp::as<Rcpp::Reference>(P).field("type");
        if (!body.isNotNull()) {
            dim = Rcpp::as<Rcpp::Reference>(P).field("dimension");
            if (type == 3) {
                ZP.init(dim, Rcpp::as<MT>(Rcpp::as<Rcpp::Reference>(P).field("G")),
                        VT::Ones(Rcpp::as<MT>(Rcpp::as<Rcpp::Reference>(P).field("G")).rows()));
                vol = exact_zonotope_vol<NT>(ZP);
            } else if (type == 2) {
                if (Rcpp::as<MT>(Rcpp::as<Rcpp::Reference>(P).field("V")).rows() == Rcpp::as<MT>(Rcpp::as<Rcpp::Reference>(P).field("V")).cols()+1) {
                    MT V = Rcpp::as<MT>(Rcpp::as<Rcpp::Reference>(P).field("V")).transpose();
                    VT v0 = V.col(dim);
                    VT V2 = V.block(0, 0, dim, dim);
                    V2 = V2.colwise() - v0;

                    vol = std::abs(V2.determinant());
                    vol = vol / factorial(NT(dim));
                } else {
                    throw Rcpp::exception("Not a simplex!");
                }
            } else {
                throw Rcpp::exception("Wrong input!");
            }
        } else {
            throw Rcpp::exception("When a polytope is given, input for body has to be null.");
        }
    } else {
        if (body.isNotNull()) {
            if (Parameters.isNotNull()) {
                if (Rcpp::as<Rcpp::List>(Parameters).containsElementNamed("dimension")) {
                    dim = Rcpp::as<int>(Rcpp::as<Rcpp::List>(Parameters)["dimension"]);
                } else {
                    throw Rcpp::exception("You have to declare the dimension!");
                }
            } else {
                throw Rcpp::exception("You have to declare the dimension!");
            }
            if (Rcpp::as<std::string>(body).compare(std::string("simplex"))==0) {
                vol = 1.0 / factorial(NT(dim));
            } else if (Rcpp::as<std::string>(body).compare(std::string("cross"))==0) {
                vol = std::pow(2.0, NT(dim));
                vol = vol / factorial(NT(dim));
            } else if (Rcpp::as<std::string>(body).compare(std::string("cube"))==0) {
                vol = std::pow(2.0, NT(dim));
            } else if (Rcpp::as<std::string>(body).compare(std::string("hypersphere"))==0) {
                if (Rcpp::as<Rcpp::List>(Parameters).containsElementNamed("radius")) {
                    rad = Rcpp::as<NT>(Rcpp::as<Rcpp::List>(Parameters)["radius"]);
                } else {
                    throw Rcpp::exception("You have to declare the dimension!");
                }
                vol = (std::pow(M_PI,dim/2.0)*(std::pow(rad, dim) ) ) / (tgamma(dim/2.0+1));
            } else {
                throw Rcpp::exception("Unknown body!");
            }

        } else {
            throw Rcpp::exception("When a polytope is null, input for specific body required.");
        }
    }

    return vol;
}
