// VolEsti (volume computation and sampling library)

// Copyright (c) 2018 Vissarion Fisikopoulos, Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.
//Contributed and/or modified by Panagiotis Repouskos, as part of Google Summer of Code 2019 program.

// VolEsti is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or (at
// your option) any later version.
//
// VolEsti is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
//
// See the file COPYING.LESSER for the text of the GNU Lesser General
// Public License.  If you did not receive this file along with HeaDDaCHe,
// see <http://www.gnu.org/licenses/>.

#ifndef VOLESTI_LP_GENERATOR_H
#define VOLESTI_LP_GENERATOR_H

#include "lp_problem.h"
#include "polytope_generators.h"
#include "polytopes.h"
#include "Eigen"
#include <cstdlib>

typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VT;

void generate_objective_function(const int& dim, VT& objFunction) {
    objFunction.resize(dim);

    for (int i = 0; i < dim; i++) {
        objFunction(i) = -5 + ((double)rand() / RAND_MAX)*30;
    }
}

template<class Point, typename NT>
optimization::lp_problem<Point, NT> generate_lp_cube(int dim) {
    VT objFunction;
    generate_objective_function(dim, objFunction);

    HPolytope<Point> HP = gen_cube<HPolytope<Point> >(dim, false);

    return optimization::lp_problem<Point, NT>(HP, objFunction, optimization::minimize);
}

template<class Point, typename NT>
optimization::lp_problem<Point, NT> generate_lp_cross(int dim) {
    VT objFunction;
    generate_objective_function(dim, objFunction);

    HPolytope<Point> HP = gen_cross<HPolytope<Point> >(dim, false);

    return optimization::lp_problem<Point, NT>(HP, objFunction, optimization::minimize);
}

template<class Point, typename NT>
optimization::lp_problem<Point, NT> generate_lp_simplex(int dim) {
    VT objFunction;
    generate_objective_function(dim, objFunction);

    HPolytope<Point> HP = gen_simplex<HPolytope<Point> >(dim, false);

    return optimization::lp_problem<Point, NT>(HP, objFunction, optimization::minimize);
}

template<class Point, typename NT>
optimization::lp_problem<Point, NT> generate_lp_prod_simplex(int dim) {
    VT objFunction;
    generate_objective_function(dim, objFunction);

    HPolytope<Point> HP = gen_prod_simplex<HPolytope<Point> >(dim, false);

    return optimization::lp_problem<Point, NT>(HP, objFunction, optimization::minimize);
}

template<class Point, typename NT>
optimization::lp_problem<Point, NT> generate_lp_skinny_cube(int dim) {
    VT objFunction;
    generate_objective_function(dim, objFunction);

    HPolytope<Point> HP = gen_skinny_cube<HPolytope<Point> >(dim, false);

    return optimization::lp_problem<Point, NT>(HP, objFunction, optimization::minimize);
}

template<class Point, typename NT, class RNGType>
optimization::lp_problem<Point, NT> generate_lp_zonotope(int dim, int m) {
    VT objFunction;
    generate_objective_function(dim, objFunction);

    HPolytope<Point> HP = gen_zonotope<HPolytope<Point>, RNGType>(dim, m);

    return optimization::lp_problem<Point, NT>(HP, objFunction, optimization::minimize);
}

template<class Point, typename NT, class RNGType>
optimization::lp_problem<Point, NT> generate_lp(int dim, int m) {
    VT objFunction;
    generate_objective_function(dim, objFunction);

    HPolytope<Point> HP = random_hpoly<HPolytope<Point>, RNGType>(dim, m);

    return optimization::lp_problem<Point, NT>(HP, objFunction, optimization::minimize);
}

#endif //VOLESTI_LP_GENERATOR_H
