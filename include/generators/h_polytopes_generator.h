// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef H_POLYTOPES_GEN_H
#define H_POLYTOPES_GEN_H

#include <exception>
#include <chrono>

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
        b(i) = 0.1;
    }
    /*
        // Print each column of A
    for (unsigned int j = 0; j < dim; ++j) {
        std::cout << "Column " << j << ": ";
        for (unsigned int i = 0; i < m; ++i) {
            std::cout << A(i, j) << " ";
        }
        std::cout << std::endl;
    }

    // Print elements of b
    std::cout << "Vector b: ";
    for (unsigned int i = 0; i < m; ++i) {
        std::cout << b(i) << " ";
    }
    std::cout << std::endl;
    */
    return Polytope(dim, A, b);

    return Polytope(dim, A, b);
}

#endif
