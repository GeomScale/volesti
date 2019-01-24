// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// rcpp_module.cpp: Rcpp R/C++ interface class library -- Rcpp Module examples
//
// Copyright (C) 2010 - 2012  Dirk Eddelbuettel and Romain Francois
//
// This file is part of Rcpp.
//
// Rcpp is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Rcpp is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Rcpp.  If not, see <http://www.gnu.org/licenses/>.

#include <Rcpp.h>

class Hpolytope {
public:
    Hpolytope() {}
    Hpolytope(Rcpp::NumericMatrix _A, Rcpp::NumericVector _b) : A(_A), b(_b) {
        dimension = _A.ncol();
    }
    int type = 1;
    unsigned int dimension;
    Rcpp::NumericMatrix A;
    Rcpp::NumericVector b;

};

class Vpolytope {
public:
    Vpolytope() {}
    Vpolytope(Rcpp::NumericMatrix _V) : V(_V) {
        dimension = _V.ncol();
    }
    int type = 2;
    unsigned int dimension;
    Rcpp::NumericMatrix V;

};

class Zonotope {
public:
    Zonotope() {}
    Zonotope(Rcpp::NumericMatrix _G) : G(_G) {
        dimension = _G.ncol();
    }
    int type = 3;
    unsigned int dimension;
    Rcpp::NumericMatrix G;

};

class VPinterVP {
public:
    VPinterVP() {}
    VPinterVP(Rcpp::NumericMatrix _V1, Rcpp::NumericMatrix _V2) : V1(_V1), V2(_V2) {
        dimension = _V1.ncol();
    }
    int type = 4;
    unsigned int dimension;
    Rcpp::NumericMatrix V1;
    Rcpp::NumericMatrix V2;

};



RCPP_MODULE(yada){
    using namespace Rcpp ;

    class_<Hpolytope>("Hpolytope")
    // expose the default constructor
    .constructor()
    .constructor<Rcpp::NumericMatrix, Rcpp::NumericVector>()

    .field( "type", &Hpolytope::type )
    .field( "dimension", &Hpolytope::dimension )
    .field( "b", &Hpolytope::b )
    .field( "A", &Hpolytope::A );

    class_<Vpolytope>("Vpolytope")
    // expose the default constructor
    .constructor()
    .constructor<Rcpp::NumericMatrix>()

    .field( "type", &Vpolytope::type )
    .field( "dimension", &Vpolytope::dimension )
    .field( "V", &Vpolytope::V );

    class_<Zonotope>("Zonotope")
    // expose the default constructor
    .constructor()
    .constructor<Rcpp::NumericMatrix>()

    .field( "type", &Zonotope::type )
    .field( "dimension", &Zonotope::dimension )
    .field( "G", &Zonotope::G );

    class_<VPinterVP>("IntVP")
    // expose the default constructor
    .constructor()
    .constructor<Rcpp::NumericMatrix, Rcpp::NumericMatrix>()

    .field( "type", &VPinterVP::type )
    .field( "dimension", &VPinterVP::dimension )
    .field( "V1", &VPinterVP::V1 )
    .field( "V2", &VPinterVP::V2 );
}
