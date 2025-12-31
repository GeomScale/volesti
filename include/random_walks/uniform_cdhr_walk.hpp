// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

// Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef RANDOM_WALKS_UNIFORM_CDHR_WALK_HPP
#define RANDOM_WALKS_UNIFORM_CDHR_WALK_HPP
#include <cassert>
#include "sampling/sphere.hpp"

// coordinate directions hit-and-run walk with uniform target distribution

struct CDHRWalk
{
    struct parameters {};
    parameters param;

template
<
    typename Polytope,
    typename RandomNumberGenerator
>

// NOTE:
// This implementation of CDHR (Coordinate Direction Hit-and-Run)
// is stateful. The walk depends on the previously selected coordinate
// direction and previous point (_p_prev) when computing line intersections.
// This behavior should be preserved when restarting chains.
struct Walk
{
    typedef typename Polytope::PointType Point;
    typedef typename Point::FT NT;

    template <typename GenericPolytope>
    Walk(GenericPolytope& P, Point const& p, RandomNumberGenerator& rng)
    {
        initialize(P, p, rng);
    }

    template <typename GenericPolytope>
    Walk(GenericPolytope& P, Point const& p,
         RandomNumberGenerator& rng, parameters const& params)
    {
        initialize(P, p, rng);
    }

    template
    <
        typename BallPolytope
    >
    inline void apply(BallPolytope const& P,
                      Point &p,   // a point to start
                      unsigned int const& walk_length,
                      RandomNumberGenerator &rng)
    {
        for (auto j = 0u; j < walk_length; ++j)
    {
        auto rand_coord_prev = _rand_coord;
        _rand_coord = rng.sample_uidist();

        // Ensure sampled coordinate is valid
        assert(_rand_coord < _p.dimension());

        NT kapa = rng.sample_urdist();

        std::pair<NT, NT> bpair =
            P.line_intersect_coord(_p,
                                _p_prev,
                                _rand_coord,
                                rand_coord_prev,
                                _lamdas);

        // Intersection interval must be valid
        assert(bpair.first <= bpair.second);

        _p_prev = _p;
        _p.set_coord(
            _rand_coord,
            _p[_rand_coord] + bpair.first + kapa * (bpair.second - bpair.first)
        );
    }
        p = _p;
    }

private :

    template <typename BallPolytope>
    inline void initialize(BallPolytope const& P,
                        Point const& p,
                        RandomNumberGenerator &rng)
    {
        // CDHR requires a non-empty point dimension
        assert(p.dimension() > 0);

        _lamdas.setZero(P.num_of_hyperplanes());

        _rand_coord = rng.sample_uidist();
        // Ensure sampled coordinate is valid
        assert(_rand_coord < p.dimension());

        NT kapa = rng.sample_urdist();
        _p = p;

        std::pair<NT, NT> bpair =
            P.line_intersect_coord(_p, _rand_coord, _lamdas);

        // Intersection interval must be valid
        assert(bpair.first <= bpair.second);

        _p_prev = _p;
        _p.set_coord(
            _rand_coord,
            _p[_rand_coord] + bpair.first + kapa * (bpair.second - bpair.first)
        );
    }


    unsigned int _rand_coord;
    Point _p;
    Point _p_prev;
    typename Point::Coeff _lamdas;
};

};


#endif // RANDOM_WALKS_UNIFORM_CDHR_WALK_HPP
