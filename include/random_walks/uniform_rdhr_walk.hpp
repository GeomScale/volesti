// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

// Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef RANDOM_WALKS_UNIFORM_RDHR_WALK_HPP
#define RANDOM_WALKS_UNIFORM_RDHR_WALK_HPP

#include "sampling/sphere.hpp"

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

private:
    static constexpr NT eps = NT(1e-12);
    static constexpr unsigned int max_direction_tries = 10;

public:
    template <typename GenericPolytope>
    Walk(GenericPolytope& P, Point const& p, RandomNumberGenerator& rng)
    {
        initialize(P, p, rng);
    }

    template <typename GenericPolytope>
    Walk(GenericPolytope& P, Point const& p,
         RandomNumberGenerator& rng, parameters const&)
    {
        initialize(P, p, rng);
    }

    template <typename BallPolytope>
    inline void apply(BallPolytope& P,
                      Point& p,
                      unsigned int const& walk_length,
                      RandomNumberGenerator& rng)
    {
        for (unsigned int j = 0; j < walk_length; ++j)
        {
            bool moved = false;

            for (unsigned int tries = 0; tries < max_direction_tries; ++tries)
            {
                Point v = GetDirection<Point>::apply(p.dimension(), rng);
                auto bpair = P.line_intersect(_p, v, _lamdas, _Av, _lambda);

                if (bpair.first > bpair.second + eps)
                {
                    _lambda = rng.sample_urdist() * (bpair.first - bpair.second)
                            + bpair.second;
                    _p += (_lambda * v);
                    moved = true;
                    break;
                }
            }

            // If no valid direction was found, keep the current point
            if (!moved)
            {
                // no-op
            }
        }
        p = _p;
    }

private:
    template <typename BallPolytope>
    inline void initialize(BallPolytope& P,
                           Point const& p,
                           RandomNumberGenerator& rng)
    {
        _lamdas.setZero(P.num_of_hyperplanes());
        _Av.setZero(P.num_of_hyperplanes());

        bool initialized = false;

        for (unsigned int tries = 0; tries < max_direction_tries; ++tries)
        {
            Point v = GetDirection<Point>::apply(p.dimension(), rng);
            auto bpair = P.line_intersect(p, v, _lamdas, _Av);

            if (bpair.first > bpair.second + eps)
            {
                _lambda = rng.sample_urdist() * (bpair.first - bpair.second)
                        + bpair.second;
                _p = (_lambda * v) + p;
                initialized = true;
                break;
            }
        }

        if (!initialized)
        {
            _p = p; // safe fallback
        }
    }

    Point _p;
    NT _lambda;
    typename Point::Coeff _lamdas;
    typename Point::Coeff _Av;
};

};

#endif // RANDOM_WALKS_UNIFORM_RDHR_WALK_HPP

