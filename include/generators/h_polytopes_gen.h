// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef H_POLYTOPES_GEN_H
#define H_POLYTOPES_GEN_H

#include <exception>
#include "samplers.h"

template <class Polytope, class RNGType>
Polytope random_hpoly(unsigned int dim, unsigned int m) {

    typedef typename Polytope::MT    MT;
    typedef typename Polytope::VT    VT;
    typedef typename Polytope::NT    NT;
    typedef typename Polytope::PolytopePoint Point;

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(seed);
    boost::random::uniform_real_distribution<> urdist1(-10, 10);
    Point p(dim);
    typename std::vector<NT>::iterator pit;
    MT A(m, dim);
    VT b(m);
    unsigned int j;

    for(unsigned int i=0; i<m; ++i){
        p = get_direction<RNGType, Point, NT>(dim);
        pit = p.iter_begin();
        j = 0;
        for ( ;  pit!=p.iter_end(); ++pit, ++j) {
            A(i,j) = *pit;
        }
        b(i) = 10.0;

    }
    Polytope HP;
    HP.init(dim, A, b);

    return HP;
}

#endif
