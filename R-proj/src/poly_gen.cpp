// [[Rcpp::depends(BH)]]

// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

#include <Rcpp.h>
#include <RcppEigen.h>
#include <chrono>
#include "cartesian_geom/cartesian_kernel.h"
#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include "polytopes.h"
#include "polytope_generators.h"
#include "extractMatPoly.h"

//' An internal Rccp function as a polytope generator
//'
//' @param kind_gen An integer to declare the type of the polytope.
//' @param Vpoly_gen A boolean parameter to declare if the requested polytope has to be in V-representation.
//' @param dim_gen An integer to declare the dimension of the requested polytope.
//' @param m_gen An integer to declare the number of generators for the requested random zonotope
//'
//' @return A numerical matrix describing the requested polytope
// [[Rcpp::export]]
Rcpp::NumericMatrix poly_gen (int kind_gen, bool Vpoly_gen, int dim_gen, int m_gen){

    typedef double NT;
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef boost::mt19937    RNGType;
    typedef HPolytope<Point> Hpolytope;
    typedef VPolytope<Point, RNGType > Vpolytope;
    typedef Zonotope<Point> zonotope;
    //typedef copula_ellipsoid<Point> CopEll;

    Hpolytope HP;
    Vpolytope VP;
    zonotope ZP;

    Rcpp::NumericMatrix Mat;
    if (kind_gen == 0) {
        zonotope ZP = gen_zonotope<zonotope, RNGType>(dim_gen, m_gen);
        Mat = extractMatPoly(ZP);
    } else if (Vpoly_gen) {
        if (kind_gen == 1) {
            VP = gen_cube<Vpolytope>(dim_gen, true);
            Mat = extractMatPoly(VP);
        } else if (kind_gen == 2) {
            VP = gen_cross<Vpolytope>(dim_gen, true);
            Mat = extractMatPoly(VP);
        } else if (kind_gen == 3) {
            VP = gen_simplex<Vpolytope>(dim_gen, true);
            Mat = extractMatPoly(VP);
        } else {
            return Mat;
        }
    } else {
        if (kind_gen == 1) {
            HP = gen_cube<Hpolytope>(dim_gen, false);
            Mat = extractMatPoly(HP);
        } else if (kind_gen == 2) {
            HP = gen_cross<Hpolytope>(dim_gen, false);
            Mat = extractMatPoly(HP);
        } else if (kind_gen == 3) {
            HP = gen_simplex<Hpolytope>(dim_gen, false);
            Mat = extractMatPoly(HP);
        } else if (kind_gen == 4) {
            HP = gen_prod_simplex<Hpolytope>(dim_gen);
            Mat = extractMatPoly(HP);
        } else if (kind_gen == 5) {
            HP = gen_skinny_cube<Hpolytope>(dim_gen);
            Mat = extractMatPoly(HP);
        } else {
            return Mat;
        }
    }
    return Mat;

}
