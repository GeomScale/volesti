// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef H_POLYTOPES_GEN_H
#define H_POLYTOPES_GEN_H

#include <exception>
#include <chrono>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>

#include "preprocess/inscribed_ellipsoid_rounding.hpp"

#ifndef isnan
  using std::isnan;
#endif

/// This function generates a random H-polytope of given dimension and number of hyperplanes $m$
/// @tparam Polytope Type of returned polytope
/// @tparam RNGType RNGType Type
template <class Polytope, class RNGType>
Polytope random_hpoly(unsigned int dim, unsigned int m, int seed = std::numeric_limits<int>::signaling_NaN()) {

    typedef typename Polytope::MT    MT;
    typedef typename Polytope::VT    VT;
    typedef typename Polytope::NT    NT;
    typedef typename Polytope::PointType Point;

    int rng_seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(rng_seed);
    if (!isnan(seed)) {
        int rng_seed = seed;
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

/// This function generates a transformation that maps the unit ball to a skinny ellipsoid
/// with given ratio between the lengths of its maximum and minimum axis
template <class MT, class VT, class RNGType, typename NT>
MT get_skinny_transformation(const int d, NT const eig_ratio, int const seed)
{
    boost::normal_distribution<> gdist(0, 1);
    RNGType rng(seed);
    
    MT W(d, d);
    for (int i = 0; i < d; i++) {
        for (int j = 0; j < d; j++) {
            W(i, j) = gdist(rng);
        }
    }
    
    Eigen::HouseholderQR<MT> qr(W);
    MT Q = qr.householderQ();
    
    VT diag(d);
    const NT eig_min = NT(1), eig_max = eig_ratio;
    diag(0) = eig_min;
    diag(d-1) = eig_max;
    boost::random::uniform_real_distribution<NT> udist(NT(0), NT(1));
    NT rand;
    for (int i = 1; i < d-1; i++) {
        rand = udist(rng);
        diag(i) = rand * eig_max + (NT(1)-rand) * eig_min;
    }
    std::sort(diag.begin(), diag.end());
    MT cov = Q * diag.asDiagonal() * Q.transpose();

    return cov;
}

/// This function generates a skinny random H-polytope of given dimension and number of hyperplanes $m$
/// @tparam Polytope Type of returned polytope
/// @tparam NT Number type
/// @tparam RNGType RNGType Type
template <class Polytope, typename NT, class RNGType>
Polytope skinny_random_hpoly(unsigned int dim, unsigned int m, const bool pre_rounding = false,
                             const NT eig_ratio = NT(1000.0), int seed = std::numeric_limits<int>::signaling_NaN()) {

    typedef typename Polytope::MT    MT;
    typedef typename Polytope::VT    VT;
    typedef typename Polytope::PointType Point;

    int rng_seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(rng_seed);
    if (!isnan(seed)) {
        int rng_seed = seed;
        rng.seed(rng_seed);
    }

    Polytope P = random_hpoly<Polytope, RNGType>(dim, m, seed);

    // rounding the polytope before applying the skinny transformation
    if (pre_rounding) {
        Point x0(dim);
        // run only one iteration
        inscribed_ellipsoid_rounding<MT, VT, NT>(P, x0, 1);
    }

    MT cov = get_skinny_transformation<MT, VT, RNGType, NT>(dim, eig_ratio, seed);
    Eigen::LLT<MT> lltOfA(cov);
    MT L = lltOfA.matrixL();
    P.linear_transformIt(L.inverse());

    return P;
}



#endif
