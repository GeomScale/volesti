// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2021 Vissarion Fisikopoulos
// Copyright (c) 2018-2021 Apostolos Chalkis

// Contributed and/or modified by Konstantinos Pallikaris, as part of Google Summer of Code 2021 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef RANDOM_WALKS_UNIFORM_RDHR_WALK_MULTITHREAD_HPP
#define RANDOM_WALKS_UNIFORM_RDHR_WALK_MULTITHREAD_HPP


#include "sampling/sphere.hpp"

// Random directions hit-and-run walk with uniform target distribution
// Parallel version

struct RDHRWalk_multithread
{

    template<typename NT, typename Point>
    struct thread_parameters
    {
        thread_parameters(unsigned int d, unsigned int m)
        {
            p = Point(d);
            v = Point(d);
            lambdas.setZero(m);
            Av.setZero(m);
            lambda_prev = NT(0);
        }

        Point p;
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
            params.p += (params.lambda_prev * params.v);
        }
    }

private :

    template <typename BallPolytope>
    inline void initialize(BallPolytope const& P,
                           thread_params &params, // parameters 
                           RandomNumberGenerator &rng)
    {
        params.lamdas.setZero(P.num_of_hyperplanes());
        params.Av.setZero(P.num_of_hyperplanes());

        params.v = GetDirection<Point>::apply(P.dimension(), rng);
        std::pair<NT, NT> bpair = P.line_intersect(params.p, params.v, params.lamdas, params.Av);
        params.lambda_prev = rng.sample_urdist() * (bpair.first - bpair.second) + bpair.second;
        params.p += (params.lambda_prev * params.v);
    }

};

};


#endif // RANDOM_WALKS_UNIFORM_RDHR_WALK_HPP
