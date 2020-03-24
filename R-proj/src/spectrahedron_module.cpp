//
// Created by panagiotis on 3/24/20.
//

#include <Rcpp.h>

class Spectrahedron {
public:
    /// A list with the matrices A0, ..., An
    Rcpp::List matrices;

    Spectrahedron() {}

    Spectrahedron(Rcpp::List _matrices) : matrices(_matrices) {}
};


RCPP_MODULE(spec){
        using namespace Rcpp ;

        //' An exposed class to represent a Spectrahedron
        //'
        //' @description A Spectrahedron is a convex polytope defined by a linear matrix inequality of the form \eqn{A_0 + x_1 A_1 + ... + x_n A_n \prec 0}.
        //'
        //' @field matrices The matrices \eqn{A_0, A_1, ..., A_n}
        //' @export
        class_<Spectrahedron>("Spectrahedron")
        // expose the default constructor
        .constructor()
        .constructor<Rcpp::List>()
        .field( "matrices", &Spectrahedron::matrices);
}

extern SEXP _rcpp_module_boot_spec(void);
