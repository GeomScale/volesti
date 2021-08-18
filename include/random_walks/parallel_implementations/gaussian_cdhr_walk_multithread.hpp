// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

// Contributed and/or modified by Konstantinos Pallikaris, as part of Google Summer of Code 2021 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef RANDOM_WALKS_GAUSSIAN_CDHR_WALK_MULTITHREAD_HPP
#define RANDOM_WALKS_GAUSSIAN_CDHR_WALK_MULTITHREAD_HPP

#include "sampling/sphere.hpp"
#include "generators/boost_random_number_generator.hpp"
#include "random_walks/gaussian_helpers.hpp"
#include "random_walks/gaussian_cdhr_walk.hpp"


// Coordinate directions hit-and-run walk with spherical Gaussian target distribution
// Multithred version

struct GaussianCDHRWalk_multithread
{

    template<typename NT, typename Point>
    struct thread_parameters
    {
        thread_parameters(unsigned int d, unsigned int m)
        {
            p = Point(d);
            p_prev = Point(d);
            lambdas.setZero(m);
        }

        Point p;
        Point p_prev;
        unsigned int rand_coord_prev;
        unsigned int rand_coord;
        typename Point::Coeff lambdas;
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

    Walk(Polytope const& P,
         thread_params &params,
         NT const& a_i,
         RandomNumberGenerator &rng)
    {
        initialize(P, params, a_i, rng);
    }
    

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
            params.rand_coord_prev = params.rand_coord;
            params.rand_coord = rng.sample_uidist();
            std::pair <NT, NT> bpair =
                    P.line_intersect_coord(params.p, params.p_prev, params.rand_coord,
                                           params.rand_coord_prev, params.lamdas);
            NT dis = chord_random_point_generator_exp_coord
                        (params.p[params.rand_coord] + bpair.second,
                         params.p[params.rand_coord] + bpair.first,
                         a_i,
                         rng);
            params.p_prev = params.p;
            params.p.set_coord(params.rand_coord, dis);
        }
    }

private :

    template <typename BallPolytope>
    inline void initialize(BallPolytope const& P,
                           thread_params &params, // parameters 
                           NT const& a_i,
                           RandomNumberGenerator &rng)
    {
        params.lamdas.setZero(P.num_of_hyperplanes());
        params.rand_coord = rng.sample_uidist();

        std::pair <NT, NT> bpair = P.line_intersect_coord(params.p, params.rand_coord, params.lamdas);
        NT dis = chord_random_point_generator_exp_coord
                    (params.p[params.rand_coord] + bpair.second,
                     params.p[params.rand_coord] + bpair.first,
                     a_i, rng);
        params.p_prev = params.p;
        params.p.set_coord(params.rand_coord, dis);
    }

};

};

#endif // RANDOM_WALKS_GAUSSIAN_CDHR_WALK_MULTITHREAD_HPP
