// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

// Contributed and/or modified by Konstantinos Pallikaris, as part of Google Summer of Code 2021 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef RANDOM_WALKS_GAUSSIAN_RDHR_WALK_MULTITHREAD_HPP
#define RANDOM_WALKS_GAUSSIAN_RDHR_WALK_MULTITHREAD_HPP

#include "random_walks/gaussian_helpers.hpp"
#include "random_walks/gaussian_rdhr_walk.hpp"
#include "generators/boost_random_number_generator.hpp"


// Random directions hit-and-run walk with spherical Gaussian target distribution
// multithread version

struct GaussianRDHRWalk_multithread
{

    template<typename NT, typename Point>
    struct thread_parameters
    {
        thread_parameters(unsigned int d, unsigned int m)
        {
            p = Point(d);
            v = Point(d);
        }

        Point p;
        Point v;
    };

template
<
    typename Polytope,
    typename RandomNumberGenerator
>
struct Walk
{
    typedef typename Polytope::PointType Point;
    typedef typename Point::FT NT;

    Walk(Polytope const&, Point const&, NT const&, RandomNumberGenerator&)
    {}

    Walk(Polytope const&, Point const&, NT const&, RandomNumberGenerator&,
         parameters&)
    {}

    template
    <
        typename BallPolytope
    >
    inline void apply(BallPolytope const& P,
                      thread_params &params, // parameters 
                      NT const& a_i,
                      unsigned int const& walk_length,
                      RandomNumberGenerator &rng)
    {
        for (auto j = 0u; j < walk_length; ++j)
        {
            params.v = GetDirection<Point>::apply(P.dimension(), rng);
            std::pair <NT, NT> dbpair = P.line_intersect(params.p, params.v);

            NT min_plus = dbpair.first;
            NT max_minus = dbpair.second;
            Point upper = (min_plus * params.v) + params.p;
            Point lower = (max_minus * params.v) + params.p;

            chord_random_point_generator_exp(lower, upper, a_i, params.p, rng);
        }
    }
};

};


#endif // RANDOM_WALKS_GAUSSIAN_RDHR_WALK_MULTITHREAD_HPP
