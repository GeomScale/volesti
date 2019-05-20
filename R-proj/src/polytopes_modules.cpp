// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2019 Vissarion Fisikopoulos
// Copyright (c) 2018-2019 Apostolos Chalkis

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


    //' An exposed class to represent a H-polytope
    //'
    //' @description A H-polytope is a convex polytope defined by a set of linear inequalities or equivalently a \eqn{d}-dimensional H-polytope with \eqn{m} facets is defined by a \eqn{m\times d} matrix A and a \eqn{m}-dimensional vector b, s.t.: \eqn{Ax\leq b}.
    //'
    //' @field A \eqn{m\times d} numerical matrix A
    //' @field b \eqn{m}-dimensional vector b
    //' @field type An integer that declares the representation of the polytope. For H-representation the default value is 1. It has not be given to the constructor.
    //' @field An integer that declares the dimension of the polytope. It has not be given to the constructor.
    //'
    //' @example
    //' # create a 2-d unit simplex
    //' A = matrix(c(-1,0,0,-1,1,1), ncol=2, nrow=3, byrow=TRUE)
    //' b = c(0,0,1)
    //' P = Hpolytope$new(A,b)
    //' @export
    class_<Hpolytope>("Hpolytope")
    // expose the default constructor
    .constructor()
    .constructor<Rcpp::NumericMatrix, Rcpp::NumericVector>()

    .field( "type", &Hpolytope::type )
    .field( "dimension", &Hpolytope::dimension )
    .field( "b", &Hpolytope::b )
    .field( "A", &Hpolytope::A );

    //' An exposed C++ class to represent a V-polytope
    //'
    //' @description A V-polytope is a convex polytope defined by the set of its vertices.
    //'
    //' @field V \eqn{m\times d} numerical matrix that contains the vertices row-wise
    //' @field type An integer that declares the representation of the polytope. For H-representation the default value is 2. It has not be given to the constructor.
    //' @field An integer that declares the dimension of the polytope. It has not be given to the constructor.
    //'
    //' @example
    //' # Create a 2-d cube in V-representation
    //' V = matrix(c(-1,1,1,1,1,-1,-1,-1), ncol=3, byrow=TRUE)
    //' P = Vpolytope$new(V)
    //' @export
    class_<Vpolytope>("Vpolytope")
    // expose the default constructor
    .constructor()
    .constructor<Rcpp::NumericMatrix>()

    .field( "type", &Vpolytope::type )
    .field( "dimension", &Vpolytope::dimension )
    .field( "V", &Vpolytope::V );

    //' An exposed C++ class to represent a zonotope
    //'
    //' @description A zonotope is a convex polytope defined by the Minkowski sum of \eqn{m} \eqn{d}-dimensional segments.
    //'
    //' @field G \eqn{m\times d} numerical matrix that contains the segments (or generators) row-wise
    //' @field type An integer that declares the representation of the polytope. For H-representation the default value is 3. It has not be given to the constructor.
    //' @field An integer that declares the dimension of the polytope. It has not be given to the constructor.
    //'
    //' @example
    //' # Create a 2-d zonotope with 4 generators
    //' G = matrix(c(1,0,0,1,-0.73,0.67,-0.25,0.96), ncol = 2, nrow = 4, byrow = TRUE)
    //' P = Zonotope$new(G)
    //' @export
    class_<Zonotope>("Zonotope")
    // expose the default constructor
    .constructor()
    .constructor<Rcpp::NumericMatrix>()

    .field( "type", &Zonotope::type )
    .field( "dimension", &Zonotope::dimension )
    .field( "G", &Zonotope::G );

    //' An exposed C++ class to represent an intersection of two V-polytopes
    //'
    //' @description An intersection of two V-polytopes is defined by the intersection of the two coresponding convex hulls.
    //'
    //' @field G \eqn{m\times d} numerical matrix that contains the segments (or generators) row-wise
    //' @field type An integer that declares the representation of the polytope. For H-representation the default value is 4. It has not be given to the constructor.
    //' @field An integer that declares the dimension of the polytope. It has not be given to the constructor.
    //'
    //' @example
    //' # define the intwrsection of a 2-d simplex with a 2-d cross polytope
    //' P1 = GenSimplex(2,'V')
    //' P2 = GenCross(2,'V')
    //' P = IntP$new(P1$V, P2$V)
    //' @export
    class_<VPinterVP>("IntVP")
    // expose the default constructor
    .constructor()
    .constructor<Rcpp::NumericMatrix, Rcpp::NumericMatrix>()

    .field( "type", &VPinterVP::type )
    .field( "dimension", &VPinterVP::dimension )
    .field( "V1", &VPinterVP::V1 )
    .field( "V2", &VPinterVP::V2 );
}

extern SEXP _rcpp_module_boot_yada(void);
