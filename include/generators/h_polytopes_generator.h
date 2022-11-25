// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef H_POLYTOPES_GEN_H
#define H_POLYTOPES_GEN_H

#include <exception>


#ifndef isnan
  using std::isnan;
#endif

/// This function generates a random H-polytope of given dimension and number of hyperplanes $m$
/// @tparam Polytope Type of returned polytope
/// @tparam RNGType RNGType Type
template <class Polytope, class RNGType>
Polytope random_hpoly(unsigned int dim, unsigned int m, double seed = std::numeric_limits<double>::signaling_NaN()) {

    typedef typename Polytope::MT    MT;
    typedef typename Polytope::VT    VT;
    typedef typename Polytope::NT    NT;
    typedef typename Polytope::PointType Point;

    unsigned rng_seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(rng_seed);
    if (!isnan(seed)) {
        unsigned rng_seed = seed;
        rng.seed(rng_seed);
    }

    MT A(m, dim);
    VT b(m);
    Point p(dim);

    for (int i = 0; i < m; ++i) {
        boost::normal_distribution<> rdist(0, 1);
        NT normal = NT(0);
        NT *data = p.pointerToData();

        //RNGType rng2 = var.rng;
        for (unsigned int i = 0; i < dim; ++i) {
            *data = rdist(rng);
            normal += *data * *data;
            data++;
        }

        normal = 1.0 / std::sqrt(normal);
        p *= normal;
        A.row(i) = p.getCoefficients();
        b(i) = 10.0;
    }

    return Polytope(dim, A, b);
}

#endif
