// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

// Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef RANDOM_WALKS_GAUSSIAN_RDHR_WALK_HPP
#define RANDOM_WALKS_GAUSSIAN_RDHR_WALK_HPP

#include "random_walks/gaussian_helpers.hpp"
#include "generators/boost_random_number_generator.hpp"


// Pick a point from the distribution exp(-a_i||x||^2) on the chord
template
<
    typename Point,
    typename NT,
    typename RandomNumberGenerator
>
bool chord_random_point_generator_exp(Point &lower,
                                      Point & upper,
                                      const NT &a_i,
                                      Point &p,
                                      RandomNumberGenerator& rng,
                                      unsigned int max_sampling_tries)
{
    NT r, r_val, fn;
    Point bef = upper - lower;
    // pick from 1-dimensional gaussian if enough weight is inside polytope P
    if (a_i > EXP_CHORD_TOLERENCE && std::sqrt(bef.squared_length()) >= (2.0 / std::sqrt(2.0 * a_i)))
    {
        Point a = -1.0 * lower;
        Point b = (1.0 / std::sqrt(bef.squared_length())) * bef;
        Point z = (a.dot(b) * b) + lower;
        NT low_bd = (lower[0] - z[0]) / b[0];
        NT up_bd = (upper[0] - z[0]) / b[0];
        for (unsigned int tries = 0; tries < max_sampling_tries; ++tries)
        {
            r = rng.sample_ndist();
            r = r / std::sqrt(2.0 * a_i);

            if (r >= low_bd && r <= up_bd)
            {
                p = (r * b) + z;
                return true;
            }
        }
        return false;

    // select using rejection sampling from a bounding rectangle
    } else {
        NT M = get_max(lower, upper, a_i);
        for (unsigned int tries = 0; tries < max_sampling_tries; ++tries)
        {
            r = rng.sample_urdist();
            Point pef = r * upper;
            p = ((1.0 - r) * lower) + pef;

            r_val = M * rng.sample_urdist();
            fn = eval_exp(p, a_i);

            if (r_val < fn)
            {
                return true;
            }
        }
        return false;
    }
}

// Random directions hit-and-run walk with spherical Gaussian target distribution

struct GaussianRDHRWalk
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

    static constexpr NT eps = NT(1e-12);
    static constexpr unsigned int max_direction_tries = 10;
    static constexpr unsigned int max_sampling_tries = 100;

    Walk(Polytope&, Point const&, NT const&, RandomNumberGenerator&)
    {}

    Walk(Polytope&, Point const&, NT const&, RandomNumberGenerator&,
         parameters&)
    {}

    template <typename BallPolytope>
    inline void apply(BallPolytope const& P,
                    Point& p,
                    NT const& a_i,
                    unsigned int const& walk_length,
                    RandomNumberGenerator& rng)
    {
        for (auto j = 0u; j < walk_length; ++j)
        {
            bool moved = false;

            for (unsigned int tries = 0; tries < max_direction_tries; ++tries)
            {
                Point v = GetDirection<Point>::apply(p.dimension(), rng);
                auto dbpair = P.line_intersect(p, v);

                if (dbpair.first > dbpair.second + eps)
                {
                    Point upper = (dbpair.first * v) + p;
                    Point lower = (dbpair.second * v) + p;

                    bool ok = chord_random_point_generator_exp(
                        lower, upper, a_i, p, rng, max_sampling_tries);

                    if (ok)
                    {
                        moved = true;
                        break;
                    }
                }
            }

            if (!moved)
            {
                // no-op: keep p unchanged
            }
        }
    }

};

};


#endif // RANDOM_WALKS_GAUSSIAN_RDHR_WALK_HPP
