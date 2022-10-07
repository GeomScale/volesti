// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2019 Vissarion Fisikopoulos
// Copyright (c) 2018-2019 Apostolos Chalkis


#include <Rcpp.h>

class Hpolytope {
public:
    Hpolytope() {}
    Hpolytope(Rcpp::NumericMatrix _A, Rcpp::NumericVector _b) : A(_A), b(_b), Aeq(Rcpp::NumericMatrix(0,0)),
                beq(Rcpp::NumericVector(0)), vol(std::numeric_limits<double>::signaling_NaN()), dimension(_A.ncol()), type(1) {}
    Hpolytope(Rcpp::NumericMatrix _A, Rcpp::NumericVector _b, Rcpp::NumericMatrix _Aeq, Rcpp::NumericVector _beq) : 
                A(_A), b(_b), Aeq(_Aeq), beq(_beq), vol(std::numeric_limits<double>::signaling_NaN()), dimension(_A.ncol()), type(1) {}
    Hpolytope(Rcpp::NumericMatrix _A, Rcpp::NumericVector _b, double volume) : A(_A), b(_b),
                Aeq(Rcpp::NumericMatrix(0,0)), beq(Rcpp::NumericVector(0)), vol(volume), dimension(_A.ncol()), type(1) {}
    Hpolytope(Rcpp::NumericMatrix _A, Rcpp::NumericVector _b, Rcpp::NumericMatrix _Aeq, Rcpp::NumericVector _beq, 
                double volume) : A(_A), b(_b), Aeq(_Aeq), beq(_beq), vol(volume), dimension(_A.ncol()), type(1) {}
    Rcpp::NumericMatrix A;
    Rcpp::NumericVector b;
    Rcpp::NumericMatrix Aeq;
    Rcpp::NumericVector beq;
    double vol;
    unsigned int dimension;
    int type;
};

class Vpolytope {
public:
    Vpolytope() {}
    Vpolytope(Rcpp::NumericMatrix _V) : V(_V), vol(std::numeric_limits<double>::signaling_NaN()), dimension(_V.ncol()), type(2) {}
    Vpolytope(Rcpp::NumericMatrix _V, double volume) : V(_V), vol(volume), dimension(_V.ncol()), type(2) {}
    Rcpp::NumericMatrix V;
    double vol;
    unsigned int dimension;
    int type;
};

class Zonotope {
public:
    Zonotope() {}
    Zonotope(Rcpp::NumericMatrix _G) : G(_G), vol(std::numeric_limits<double>::signaling_NaN()), dimension(_G.ncol()), type(3) {}
    Zonotope(Rcpp::NumericMatrix _G, double volume) : G(_G), vol(volume), dimension(_G.ncol()), type(3) {}
    Rcpp::NumericMatrix G;
    double vol;
    unsigned int dimension;
    int type;
};

class VPinterVP {
public:
    VPinterVP() {}
    VPinterVP(Rcpp::NumericMatrix _V1, Rcpp::NumericMatrix _V2) : V1(_V1), V2(_V2), vol(std::numeric_limits<double>::signaling_NaN()), dimension(_V1.ncol()), type(4) {}
    VPinterVP(Rcpp::NumericMatrix _V1, Rcpp::NumericMatrix _V2, double volume) : V1(_V1), V2(_V2), vol(volume), dimension(_V1.ncol()), type(4) {}
    Rcpp::NumericMatrix V1;
    Rcpp::NumericMatrix V2;
    double vol;
    unsigned int dimension;
    int type;
};
#include <RcppEigen.h>
using SpMat=Eigen::SparseMatrix<double> ;
class sparse_constraint_problem {
public:
    sparse_constraint_problem() {}
    sparse_constraint_problem(SpMat _Aineq, Rcpp::NumericVector _bineq, SpMat _Aeq, Rcpp::NumericVector _beq) : Aineq(_Aineq), bineq(_bineq), Aeq(_Aeq), beq(_beq),
                        lb(Rcpp::NumericVector(_Aineq.cols(),-1e9)), ub(Rcpp::NumericMatrix(_Aineq.cols(),1e9)), dimension(_Aineq.cols()), type(5) {}
    sparse_constraint_problem(SpMat _Aineq, Rcpp::NumericVector _bineq, SpMat _Aeq, Rcpp::NumericVector _beq,
                        Rcpp::NumericVector _lb, Rcpp::NumericVector _ub) : Aineq(_Aineq), bineq(_bineq), Aeq(_Aeq), beq(_beq),
                        lb(_lb), ub(_ub), dimension(_Aineq.cols()), type(5) {}
    sparse_constraint_problem(Rcpp::NumericVector _lb, Rcpp::NumericVector _ub) : Aineq(SpMat(0, _lb.length())), bineq(Rcpp::NumericVector(0)), Aeq(SpMat(0, _lb.length())), beq(Rcpp::NumericVector(0)),
                        lb(_lb), ub(_ub), dimension(_lb.length()), type(5) {}
    SpMat Aineq;
    Rcpp::NumericVector bineq;
    SpMat Aeq;
    Rcpp::NumericVector beq;
    Rcpp::NumericVector lb;
    Rcpp::NumericVector ub;
    unsigned int dimension;
    int type;
};
RCPP_MODULE(polytopes){
    using namespace Rcpp ;

    //' An exposed class to represent a H-polytope
    //'
    //' @description A H-polytope is a convex polytope defined by a set of linear inequalities or equivalently a \eqn{d}-dimensional H-polytope with \eqn{m} facets is defined by a \eqn{m\times d} matrix A and a \eqn{m}-dimensional vector b, s.t.: \eqn{Ax\leq b}.
    //'
    //' @field A \eqn{m\times d} numerical matrix A
    //' @field b \eqn{m}-dimensional vector b
    //' @field Aeq \eqn{q\times d} numerical matrix Aeq
    //' @field beq \eqn{q}-dimensional vector beq
    //' @field volume The volume of the polytope.
    //' @field dimension An integer that declares the dimension of the polytope. It has not be given to the constructor.
    //' @field type An integer that declares the representation of the polytope. For H-representation the default value is 1. It has not be given to the constructor.
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
    .constructor<Rcpp::NumericMatrix, Rcpp::NumericVector, Rcpp::NumericMatrix, Rcpp::NumericVector>()
    .constructor<Rcpp::NumericMatrix, Rcpp::NumericVector, double>()
    .constructor<Rcpp::NumericMatrix, Rcpp::NumericVector, Rcpp::NumericMatrix, Rcpp::NumericVector, double>()

    .field( "A", &Hpolytope::A )
    .field( "b", &Hpolytope::b )
    .field( "Aeq", &Hpolytope::Aeq )
    .field( "beq", &Hpolytope::beq )
    .field( "volume", &Hpolytope::vol )
    .field( "dimension", &Hpolytope::dimension )
    .field( "type", &Hpolytope::type );

    //' An exposed C++ class to represent a V-polytope
    //'
    //' @description A V-polytope is a convex polytope defined by the set of its vertices.
    //'
    //' @field V \eqn{m\times d} numerical matrix that contains the vertices row-wise
    //' @field volume The volume of the polytope.
    //' @field dimension An integer that declares the dimension of the polytope. It has not be given to the constructor.
    //' @field type An integer that declares the representation of the polytope. For H-representation the default value is 2. It has not be given to the constructor.
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
    .constructor<Rcpp::NumericMatrix, double>()

    .field( "V", &Vpolytope::V )
    .field( "volume", &Vpolytope::vol )
    .field( "dimension", &Vpolytope::dimension )
    .field( "type", &Vpolytope::type );

    //' An exposed C++ class to represent a zonotope
    //'
    //' @description A zonotope is a convex polytope defined by the Minkowski sum of \eqn{m} \eqn{d}-dimensional segments.
    //'
    //' @field G \eqn{m\times d} numerical matrix that contains the segments (or generators) row-wise
    //' @field volume The volume of the polytope.
    //' @field dimension An integer that declares the dimension of the polytope. It has not be given to the constructor.
    //' @field type An integer that declares the representation of the polytope. For H-representation the default value is 3. It has not be given to the constructor.
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
    .constructor<Rcpp::NumericMatrix, double>()

    .field( "G", &Zonotope::G )
    .field( "volume", &Zonotope::vol )
    .field( "dimension", &Zonotope::dimension )
    .field( "type", &Zonotope::type );

    //' An exposed C++ class to represent an intersection of two V-polytopes
    //'
    //' @description An intersection of two V-polytopes is defined by the intersection of the two coresponding convex hulls.
    //'
    //' @field V1 \eqn{m\times d} numerical matrix that contains the vertices of the first V-polytope (row-wise)
    //' @field V2 \eqn{q\times d} numerical matrix that contains the vertices of the second V-polytope (row-wise)
    //' @field volume The volume of the polytope.
    //' @field dimension An integer that declares the dimension of the polytope. It has not be given to the constructor.
    //' @field type An integer that declares the representation of the polytope. For H-representation the default value is 4. It has not be given to the constructor.
    //'
    //' @example
    //' # define the intwrsection of a 2-d simplex with a 2-d cross polytope
    //' P1 = gen_simplex(2,'V')
    //' P2 = gen_cross(2,'V')
    //' P = VpolytopeIntersection$new(P1$V, P2$V)
    //' @export
    class_<VPinterVP>("VpolytopeIntersection")
    // expose the default constructor
    .constructor()
    .constructor<Rcpp::NumericMatrix, Rcpp::NumericMatrix>()
    .constructor<Rcpp::NumericMatrix, Rcpp::NumericMatrix, double>()

    .field( "V1", &VPinterVP::V1 )
    .field( "V2", &VPinterVP::V2 )
    .field( "volume", &VPinterVP::vol )
    .field( "dimension", &VPinterVP::dimension )
    .field( "type", &VPinterVP::type );

    //' An exposed class to represent a sparse constraint problem
    //'
    //' @description A constraint problem is defined by a set of linear inequalities and equalities or equivalently a \eqn{d}-dimensional constraint problem is defined by a \eqn{mineq\times d} matrix Aineq and a \eqn{mineq}-dimensional vector bineq, s.t.: \eqn{Aineqx\leq bineq}, a \eqn{meq\times d} matrix Aeq and a \eqn{meq}-dimensional vector beq, s.t.: \eqn{Aeqx\eq beq} and two \eqn{d} vectors lb, ub such that \eqn{lb\leq x \leq ub}.
    //'
    //' @field Aineq \eqn{mineq\times d} sparse matrix Aineq
    //' @field bineq \eqn{mineq}-dimensional vector bineq
    //' @field Aeq \eqn{meq\times d} sparse matrix Aeq
    //' @field beq \eqn{meq}-dimensional vector beq
    //' @field lb \eqn{d}-dimensional vector lb
    //' @field ub \eqn{d}-dimensional vector ub
    //' @field dimension An integer that declares the dimension of the polytope. It has not be given to the constructor.
    //' @field type An integer that declares the representation of the polytope. For sparse_constraint_problem the default value is 5. It has not be given to the constructor.
    //'
    //' @example
    //' # create a 2-d unit simplex
    //' Aineq = matrix(, nrow=0, ncol=2, byrow = TRUE)
    //' Aineq = as( Aineq, 'dgCMatrix' )
    //' bineq= matrix(,nrow=0, ncol=1, byrow=TRUE)
    //' Aeq = matrix(c(1,1), ncol=2, nrow=1, byrow=TRUE)
    //' Aeq = as( Aeq, 'dgCMatrix' )
    //' beq = c(1)
    //' lb = c(0,0)
    //' ub = c(1,1)
    //' P = sparse_constraint_problem$new(Aineq, bineq, Aeq, beq, lb, ub)
    //' @export
    class_<sparse_constraint_problem>("sparse_constraint_problem")
    // expose the default constructor
    .constructor()
    .constructor<SpMat, Rcpp::NumericVector, SpMat, Rcpp::NumericVector>()
    .constructor<SpMat, Rcpp::NumericVector, SpMat, Rcpp::NumericVector, Rcpp::NumericVector, Rcpp::NumericVector>()
    .constructor<Rcpp::NumericVector, Rcpp::NumericVector>()

    .field( "Aineq", &sparse_constraint_problem::Aineq )
    .field( "bineq", &sparse_constraint_problem::bineq )
    .field( "Aeq", &sparse_constraint_problem::Aeq )
    .field( "beq", &sparse_constraint_problem::beq )
    .field( "lb", &sparse_constraint_problem::lb )
    .field( "ub", &sparse_constraint_problem::ub )
    .field( "dimension", &sparse_constraint_problem::dimension )
    .field( "type", &sparse_constraint_problem::type );
}

extern SEXP _rcpp_module_boot_polytopes(void);
