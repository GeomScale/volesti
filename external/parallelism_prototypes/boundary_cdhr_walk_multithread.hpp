// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

// Contributed and/or modified by Konstantinos Pallikaris, as part of Google Summer of Code 2021 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef RANDOM_WALKS_BOUNDARY_CDHR_WALK_MULTITHREAD_HPP
#define RANDOM_WALKS_BOUNDARY_CDHR_WALK_MULTITHREAD_HPP

#include "sampling/sphere.hpp"

// random directions hit-and-run walk with uniform target distribution
// from boundary

struct BCDHRWalk_multithread
{
    template<typename NT, typename Point>
    struct thread_parameters
    {
        thread_parameters(unsigned int d, unsigned int m)
        {
            p = Point(d);
            p1 = Point(d);
            p2 = Point(d);
            p_prev = Point(d);
            lambdas.setZero(m);
        }

        Point p;
        Point p1;
        Point p2;
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

        template <typename GenericPolytope>
        Walk(GenericPolytope const& P, thread_params &parameters, RandomNumberGenerator& rng)
        {
            initialize(P, parameters, rng);
        }

        template
        <
            typename BallPolytope
        >
        inline void apply(BallPolytope const& P,
                          thread_params &params, // parameters 
                          unsigned int const& walk_length,
                          RandomNumberGenerator& rng)
        {
            std::pair<NT, NT> bpair;
            for (auto j = 0u; j < walk_length; ++j)
            {
                params.rand_coord_prev = params.rand_coord;
                params.rand_coord = rng.sample_uidist();
                NT kapa = rng.sample_urdist();
                bpair = P.line_intersect_coord(params.p,
                                               params.p_prev,
                                               params.rand_coord,
                                               params.rand_coord_prev,
                                               params.lamdas);
                params.p_prev = params.p;
                params.p.set_coord(params.rand_coord, params.p[params.rand_coord] + bpair.first + kapa
                                          * (bpair.second - bpair.first));
            }
            params.p1 = params.p_prev;
            params.p2 = params.p_prev;
            params.p1.set_coord(params.rand_coord, params.p_prev[params.rand_coord] + bpair.first);
            params.p2.set_coord(params.rand_coord, params.p_prev[params.rand_coord] + bpair.second);
        }

    private :

        template <typename GenericBody>
        inline void initialize(GenericBody const& P,
                               thread_params &params, // parameters 
                               RandomNumberGenerator& rng)
        {
            params.lamdas.setZero(P.num_of_hyperplanes());
            params.rand_coord = rng.sample_uidist();
            NT kapa = rng.sample_urdist();
            
            std::pair<NT, NT> bpair = P.line_intersect_coord(params.p, params.rand_coord,
                                                             params.lamdas);
            params.p_prev = params.p;
            params.p.set_coord(params.rand_coord, params.p[params.rand_coord] + bpair.first + kapa
                                      * (bpair.second - bpair.first));
        }

    };

};

#endif // RANDOM_WALKS_BOUNDARY_CDHR_WALK_MULTITHREAD_HPP
