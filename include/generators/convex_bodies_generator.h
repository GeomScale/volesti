// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef CONVEX_BODIES_GEN_H
#define CONVEX_BODIES_GEN_H

#include <exception>

#include "convex_bodies/convex_body.h"

#ifndef isnan
  using std::isnan;
#endif

template <class ConvexBody>
ConvexBody generate_unit_ball(unsigned int dim) {

    typedef typename ConvexBody::MT    MT;
    typedef typename ConvexBody::VT    VT;
    typedef typename ConvexBody::NT    NT;
    typedef typename ConvexBody::PointType Point;
    typedef std::function<NT(const Point&)> func;
    typedef std::function<Point(const Point&)> grad;

    func unit_ball_func = [](const Point &x) {
        return x.dot(x) - 1;
    };

    grad unit_ball_grad = [](const Point &x) {
        return 2 * x;
    };

    std::vector<func> unit_ball_funcs{unit_ball_func};
    std::vector<grad> unit_ball_grads{unit_ball_grad};

    return ConvexBody(unit_ball_funcs, unit_ball_grads, dim);
}

#endif
