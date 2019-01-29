// [[Rcpp::depends(BH)]]

// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

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

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
double exact_vol(Rcpp::Nullable<Rcpp::Reference> P, Rcpp::Nullable<std::string> body = R_NilValue,
                 Rcpp::Nullable<Rcpp::List> Parameters = R_NilValue){

    typedef double NT;
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef boost::mt19937    RNGType;

    typedef Zonotope<Point> zonotope;
    typedef VPolytope <Point, RNGType> Vpolytope;
    zonotope ZP;
    Vpolytope VP;
    int type, dim;
    NT vol, rad;

    if (P.isNotNull()) {
        if (!body.isNotNull()) {
            dim = P.field("dimension");
            if (type == 3) {
                ZP.init(dim, Rcpp::as<MT>(Rcpp::as<Rcpp::Reference>(P).field("G")),
                        VT::Ones(Rcpp::as<MT>(Rcpp::as<Rcpp::Reference>(P).field("G")).rows()));
                vol = exact_zonotope_vol<NT>(ZP);
            } else if (type == 2) {
                MT V = Rcpp::as<MT>(Rcpp::as<Rcpp::Reference>(P).field("V")).transpose();
                VT v0 = V.col(dim);
                VT V2 = V.block(0,0,dim,dim);
                V2 = V2.colwise() - v0;

                vol = V2.determinant();
                vol = vol / factorial(NT(dim));
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
