// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef Z_POLYTOPES_GEN_H
#define Z_POLYTOPES_GEN_H

#include <exception>

#ifndef isnan
  using std::isnan;
#endif

/// Generates a random Zonotope with generators draw from Gaussian distribution
/// @tparam Polytope polytope type
/// @tparam RNGType RNGType type
template <class Polytope, class RNGType>
Polytope gen_zonotope_gaussian(int dim, int m, double seed = std::numeric_limits<double>::signaling_NaN()) {

    typedef typename Polytope::MT    MT;
    typedef typename Polytope::VT    VT;
    typedef typename Polytope::NT    NT;

    unsigned rng_seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(rng_seed);
    if (!isnan(seed)) {
        unsigned rng_seed = seed;
        rng.seed(rng_seed);
    }
    boost::normal_distribution<> rdist(0, 1);
    boost::normal_distribution<> rdist2(50, 33.3);

    MT A;
    VT b;
    A.resize(m, dim);
    b.resize(m);
    NT rand_gaus;

    for (unsigned int i = 0; i < m; ++i) {
        b(i) = 1.0;
        for (unsigned int j = 0; j < dim; ++j) {
            A(i,j) = rdist(rng);
        }
        A.row(i)=A.row(i)/A.row(i).norm();
        while(true){
            rand_gaus = rdist2(rng);
            if (rand_gaus > 0.0 && rand_gaus<100.0){
                A.row(i) = A.row(i) * rand_gaus;
                break;
            }
        }
    }

    Polytope P(dim, A, b);
    return P;
}


/// Generates a random Zonotope with generators draw from uniform distribution
/// @tparam Polytope polytope type
/// @tparam RNGType RNGType type
template <class Polytope, class RNGType>
Polytope gen_zonotope_uniform(int dim, int m, double seed = std::numeric_limits<double>::signaling_NaN()) {

    typedef typename Polytope::MT    MT;
    typedef typename Polytope::VT    VT;
    typedef typename Polytope::NT    NT;

    unsigned rng_seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(rng_seed);
    if (!isnan(seed)) {
        unsigned rng_seed = seed;
        rng.seed(rng_seed);
    }
    boost::normal_distribution<> rdist(0, 1);
    boost::random::uniform_real_distribution<> urdist1(0, 100);

    MT A;
    VT b;
    A.resize(m, dim);
    b.resize(m);

    for (unsigned int i = 0; i < m; ++i) {
        b(i) = 1.0;
        for (unsigned int j = 0; j < dim; ++j) {
            A(i,j) = rdist(rng);
        }
        A.row(i)=A.row(i)/A.row(i).norm();
        A.row(i) = A.row(i) * urdist1(rng);
    }

    Polytope P(dim, A, b);
    return P;

}


/// Generates a random Zonotope with generators draw from exponential distribution
/// @tparam Polytope polytope type
/// @tparam RNGType RNGType type
template <class Polytope, class RNGType>
Polytope gen_zonotope_exponential(int dim, int m, double seed = std::numeric_limits<double>::signaling_NaN()) {

    typedef typename Polytope::MT    MT;
    typedef typename Polytope::VT    VT;
    typedef typename Polytope::NT    NT;

    unsigned rng_seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(rng_seed);
    if (!isnan(seed)) {
        unsigned rng_seed = seed;
        rng.seed(rng_seed);
    }
    boost::normal_distribution<> rdist(0, 1);
    boost::normal_distribution<> expdist(1.0/30.0);

    MT A;
    VT b;
    A.resize(m, dim);
    b.resize(m);
    NT rand_exp;

    for (unsigned int i = 0; i < m; ++i) {
        b(i) = 1.0;
        for (unsigned int j = 0; j < dim; ++j) {
            A(i,j) = rdist(rng);
        }
        A.row(i)=A.row(i)/A.row(i).norm();
        while(true){
            rand_exp = expdist(rng);
            if (rand_exp > 0.0 && rand_exp<100.0){
                A.row(i) = A.row(i) * rand_exp;
                break;
            }
        }
    }

    Polytope P(dim, A, b);
    return P;
}

#endif
