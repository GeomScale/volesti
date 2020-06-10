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
#include "convex_bodies/hpolytope.h"
#include "convex_bodies/vpolytope.h"
#include "convex_bodies/zpolytope.h"
#include "generators/known_polytope_generators.h"
#include "generators/h_polytopes_generator.h"
#include "generators/v_polytopes_generators.h"
#include "generators/z_polytopes_generators.h"
#include "extractMatPoly.h"

//' An internal Rccp function as a polytope generator
//'
//' @param kind_gen An integer to declare the type of the polytope.
//' @param Vpoly_gen A boolean parameter to declare if the requested polytope has to be in V-representation.
//' @param Zono_gen A boolean parameter to declare if the requested polytope has to be a zonotope.
//' @param dim_gen An integer to declare the dimension of the requested polytope.
//' @param m_gen An integer to declare the number of generators for the requested random zonotope or the number of vertices for a V-polytope.
//' @param seed Optional. A fixed seed for the random polytope generator.
//'
//' @keywords internal
//'
//' @return A numerical matrix describing the requested polytope
// [[Rcpp::export]]
Rcpp::NumericMatrix poly_gen (int kind_gen, bool Vpoly_gen, bool Zono_gen, int dim_gen, int m_gen,
             Rcpp::Nullable<double> seed = R_NilValue) {

    typedef double NT;
    typedef Cartesian <NT> Kernel;
    typedef typename Kernel::Point Point;
    typedef boost::mt19937 RNGType;
    typedef HPolytope <Point> Hpolytope;
    typedef VPolytope <Point> Vpolytope;
    typedef Zonotope <Point> zonotope;

    double seed2 = (!seed.isNotNull()) ? std::numeric_limits<double>::signaling_NaN() : Rcpp::as<double>(seed);

    if (Zono_gen) {
        switch (kind_gen) {

            case 1:
                return extractMatPoly(gen_zonotope_uniform<zonotope, RNGType>(dim_gen, m_gen, seed2));
            case 2:
                return extractMatPoly(gen_zonotope_gaussian<zonotope, RNGType>(dim_gen, m_gen, seed2));
            case 3:
                return extractMatPoly(gen_zonotope_exponential<zonotope, RNGType>(dim_gen, m_gen, seed2));

        }

    } else if (Vpoly_gen) {
        switch (kind_gen) {

            case 1:
                return extractMatPoly(gen_cube<Vpolytope>(dim_gen, true));

            case 2:
                return extractMatPoly(gen_cross<Vpolytope>(dim_gen, true));

            case 3:
                return extractMatPoly(gen_simplex<Vpolytope>(dim_gen, true));

            case 4:
                return extractMatPoly(random_vpoly<Vpolytope, RNGType>(dim_gen, m_gen, seed2));

            case 5:
                return extractMatPoly(random_vpoly_incube<Vpolytope, RNGType>(dim_gen, m_gen, seed2));

        }
    } else {
        switch (kind_gen) {

            case 1:
                return extractMatPoly(gen_cube<Hpolytope>(dim_gen, false));

            case 2:
                return extractMatPoly(gen_cross<Hpolytope>(dim_gen, false));

            case 3:
                return extractMatPoly(gen_simplex<Hpolytope>(dim_gen, false));

            case 4:
                return extractMatPoly(gen_prod_simplex<Hpolytope>(dim_gen));

            case 5:
                return extractMatPoly(gen_skinny_cube<Hpolytope>(dim_gen));

            case 6:
                return extractMatPoly(random_hpoly<Hpolytope, RNGType>(dim_gen, m_gen, seed2));

        }
    }

    throw Rcpp::exception("Wrong inputs!");

}
