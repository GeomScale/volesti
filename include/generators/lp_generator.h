//
// Created by panagiotis on 27/5/2019.
//

#ifndef VOLESTI_LP_GENERATOR_H
#define VOLESTI_LP_GENERATOR_H

#include "lp_problem.h"
#include "polytope_generators.h"
#include "polytopes.h"
#include "Eigen"
#include <cstdlib>

template<class Point, typename NT>
optimization::lp_problem<Point, NT> generate_lp_cube(int dim) {
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
    VT obj;
    obj.resize(dim);

    for (int i = 0; i < dim; i++) {
        obj(i) = rand();
    }

    HPolytope<Point> HP = gen_cube<HPolytope<Point> >(dim, false);

    return optimization::lp_problem<Point, NT>(HP, obj, optimization::minimize);
}

template<class Point, typename NT>
optimization::lp_problem<Point, NT> generate_lp_cross(int dim) {
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
    VT obj;
    obj.resize(dim);

    for (int i = 0; i < dim; i++) {
        obj(i) = rand();
    }

    HPolytope<Point> HP = gen_cross<HPolytope<Point> >(dim, false);

    return optimization::lp_problem<Point, NT>(HP, obj, optimization::minimize);
}

template<class Point, typename NT>
optimization::lp_problem<Point, NT> generate_lp_simplex(int dim) {
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
    VT obj;
    obj.resize(dim);

    for (int i = 0; i < dim; i++) {
        obj(i) = rand();
    }

    HPolytope<Point> HP = gen_simplex<HPolytope<Point> >(dim, false);

    return optimization::lp_problem<Point, NT>(HP, obj, optimization::minimize);
}

template<class Point, typename NT>
optimization::lp_problem<Point, NT> generate_lp_prod_simplex(int dim) {
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
    VT obj;
    obj.resize(dim);

    for (int i = 0; i < dim; i++) {
        obj(i) = rand();
    }

    HPolytope<Point> HP = gen_prod_simplex<HPolytope<Point> >(dim, false);

    return optimization::lp_problem<Point, NT>(HP, obj, optimization::minimize);
}

template<class Point, typename NT>
optimization::lp_problem<Point, NT> generate_lp_skinny_cube(int dim) {
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
    VT obj;
    obj.resize(dim);

    for (int i = 0; i < dim; i++) {
        obj(i) = rand();
    }

    HPolytope<Point> HP = gen_skinny_cube<HPolytope<Point> >(dim, false);

    return optimization::lp_problem<Point, NT>(HP, obj, optimization::minimize);
}

template<class Point, typename NT, class RNGType>
optimization::lp_problem<Point, NT> generate_lp_zonotope(int dim, int m) {
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
    VT obj;
    obj.resize(dim);

    for (int i = 0; i < dim; i++) {
        obj(i) = rand();
    }

    HPolytope<Point> HP = gen_zonotope<HPolytope<Point>, RNGType>(dim, m);

    return optimization::lp_problem<Point, NT>(HP, obj, optimization::minimize);
}

template<class Point, typename NT, class RNGType>
optimization::lp_problem<Point, NT> generate_lp(int dim, int m) {
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
    VT obj;
    obj.resize(dim);

    for (int i = 0; i < dim; i++) {
        obj(i) = rand();
    }

    HPolytope<Point> HP = random_hpoly<HPolytope<Point>, RNGType>(dim, m);

    return optimization::lp_problem<Point, NT>(HP, obj, optimization::minimize);
}

#endif //VOLESTI_LP_GENERATOR_H
