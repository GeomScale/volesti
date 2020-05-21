// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

// Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef RANDOM_WALKS_UNIFORM_RDHR_WALK_HPP
#define RANDOM_WALKS_UNIFORM_RDHR_WALK_HPP

#include "convex_bodies/ball.h"
#include "convex_bodies/ballintersectconvex.h"
#include "convex_bodies/hpolytope.h"
#include "convex_bodies/vpolytope.h"
#include "convex_bodies/vpolyintersectvpoly.h"
#include "convex_bodies/zpolytope.h"
#include "convex_bodies/zonoIntersecthpoly.h"
#include "generators/boost_random_number_generator.hpp"

// Random directions hit-and-run walk with uniform target distribution

struct RDHRWalk
{
    struct parameters {};
    parameters param;

template
<
    typename Polytope,
    typename RandomNumberGenerator
>
struct Walk
{
    typedef typename Polytope::PointType Point;
    typedef typename Point::FT NT;

    template <typename GenericPolytope>
    Walk(GenericPolytope const& P, Point const& p, RandomNumberGenerator& rng)
    {
        initialize(P, p, rng);
    }

    template <typename GenericPolytope>
    Walk(GenericPolytope const& P, Point const& p,
         RandomNumberGenerator& rng, parameters const& params)
    {
        initialize(P, p, rng);
    }

    template
    <
        typename BallPolytope
    >
    inline void apply(BallPolytope const& P,
                      Point& p,   // a point to start
                      unsigned int const& walk_length,
                      RandomNumberGenerator& rng)
    {
        for (auto j=0u; j<walk_length; ++j)
        {
            Point v = GetDirection<Point>::apply(p.dimension(), rng);
            std::pair<NT, NT> bpair = P.line_intersect(_p, v, _lamdas, _Av,
                                                       _lambda);
            _lambda = rng.sample_urdist() * (bpair.first - bpair.second)
                    + bpair.second;
            _p += (_lambda * v);
        }
        p = _p;
    }

private :

    template <typename BallPolytope>
    inline void initialize(BallPolytope const& P,
                           Point const& p,
                           RandomNumberGenerator &rng)
    {
        _lamdas.setZero(P.num_of_hyperplanes());
        _Av.setZero(P.num_of_hyperplanes());

        Point v = GetDirection<Point>::apply(p.dimension(), rng);
        std::pair<NT, NT> bpair = P.line_intersect(p, v, _lamdas, _Av);
        _lambda = rng.sample_urdist() * (bpair.first - bpair.second) + bpair.second;
        _p = (_lambda * v) + p;
    }

    Point _p;
    NT _lambda;
    typename Point::Coeff _lamdas;
    typename Point::Coeff _Av;
};

};


#endif // RANDOM_WALKS_UNIFORM_RDHR_WALK_HPP
