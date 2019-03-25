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
//' @param m_gen An integer to declare the number of generators for the requested random zonotope.
//'
//' @section warning:
//' Do not use this function.
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

    Hpolytope HP;
    Vpolytope VP;
    zonotope ZP;

    Rcpp::NumericMatrix Mat;
    if (kind_gen == 0) {
        zonotope ZP = gen_zonotope<zonotope, RNGType>(dim_gen, m_gen);
        Mat = extractMatPoly(ZP);
    } else if (Vpoly_gen) {
        switch (kind_gen) {

            case 1: {
                VP = gen_cube<Vpolytope>(dim_gen, true);
                break;
            }
            case 2: {
                VP = gen_cross<Vpolytope>(dim_gen, true);
                break;
            }
            case 3: {
                VP = gen_simplex<Vpolytope>(dim_gen, true);
                break;
            }
            case 4: {
                VP = random_vpoly<Vpolytope, RNGType>(dim_gen, m_gen);
                break;
            }
        }
        Mat = extractMatPoly(VP);
    } else {
        switch (kind_gen) {

            case 1:
                HP = gen_cube<Hpolytope>(dim_gen, false);
                break;

            case 2: {
                HP = gen_cross<Hpolytope>(dim_gen, false);
                break;
            }
            case 3: {
                HP = gen_simplex<Hpolytope>(dim_gen, false);
                break;
            }
            case 4: {
                HP = gen_prod_simplex<Hpolytope>(dim_gen);
                break;
            }
            case 5: {
                HP = gen_skinny_cube<Hpolytope>(dim_gen);
                break;
            }
            case 6: {
                HP = random_hpoly<Hpolytope, RNGType>(dim_gen, m_gen);
                break;
            }
        }
        Mat = extractMatPoly(HP);
    }

    return Mat;

}
