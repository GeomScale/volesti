// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

// Contributed and/or modified by Konstantinos Pallikaris, as part of Google Summer of Code 2021 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef RANDOM_WALKS_BOUNDARY_RDHR_WALK_MULTITHREAD_HPP
#define RANDOM_WALKS_BOUNDARY_RDHR_WALK_MULTITHREAD_HPP

#include "sampling/sphere.hpp"

// Random directions hit-and-run walk with uniform target distribution
// from the boundary, multithread version

struct BRDHRWalk_multithread
{

    template<typename NT, typename Point>
    struct thread_parameters
    {
        thread_parameters(unsigned int d, unsigned int m)
        {
            p = Point(d);
            p1 = Point(d);
            p2 = Point(d);
            v = Point(d);
            lambdas.setZero(m);
            Av.setZero(m);
            lambda_prev = NT(0);
        }

        Point p;
        Point p1;
        Point p2;
        Point v;
        NT lambda_prev;
        typename Point::Coeff lambdas;
        typename Point::Coeff Av;
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
            for (auto j=0u; j<walk_length; ++j)
            {
                params.v = GetDirection<Point>::apply(P.dimension(), rng);
                std::pair<NT, NT> bpair = P.line_intersect(params.p, params.v, params.lamdas, params.Av,
                                                           params.lambda_prev);
                params.lambda_prev = rng.sample_urdist() * (bpair.first - bpair.second)
                          + bpair.second;
                params.p1 = (bpair.first * params.v);
                params.p1 += params.p;
                params.p2 = (bpair.second * params.v);
                params.p2 += params.p;
                params.p += (params.lambda_prev * params.v);
            }
        }

    private :

        template <typename GenericBody>
        inline void initialize(GenericBody const& P,
                               thread_params &params, // parameters 
                               RandomNumberGenerator& rng)
        {
            params.lamdas.setZero(P.num_of_hyperplanes());
            params.Av.setZero(P.num_of_hyperplanes());

            params.v = GetDirection<Point>::apply(P.dimension(), rng);
            std::pair<NT, NT> bpair = P.line_intersect(params.p, params.v, params.lamdas, params.Av);
            params.lambda_prev = rng.sample_urdist() * (bpair.first - bpair.second) + bpair.second;
            params.p = (params.lambda_prev * params.v) + params.p;
        }

    };

};


#endif // RANDOM_WALKS_BOUNDARY_RDHR_WALK_MULTITHREAD_HPP
