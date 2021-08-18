// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

// Contributed and/or modified by Konstantinos Pallikaris, as part of Google Summer of Code 2019 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef RANDOM_WALKS_UNIFORM_BILLIARD_WALK_MULTITHREAD_HPP
#define RANDOM_WALKS_UNIFORM_BILLIARD_WALK_MULTITHREAD_HPP

#include "sampling/sphere.hpp"



// Billiard walk for uniform distribution

struct BilliardWalk_multithread
{

    template<typename NT, typename Point>
    struct thread_parameters
    {
        thread_parameters(unsigned int d, unsigned int m)
        {
            p = Point(d);
            p0 = Point(d);
            v = Point(d);
            lambdas.setZero(m);
            Av.setZero(m);
            lambda_prev = NT(0);
        }

        Point p;
        Point p0;
        Point v;
        NT lambda_prev;
        typename Point::Coeff lambdas;
        typename Point::Coeff Av;
    };

    BilliardWalk(double L)
            :   param(L, true)
    {}

    BilliardWalk()
            :   param(0, false)
    {}

    struct parameters
    {
        parameters(double L, bool set)
                :   m_L(L), set_L(set)
        {}
        double m_L;
        bool set_L;
    };

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
    Walk(GenericPolytope const& P, thread_params &parameters, RandomNumberGenerator &rng)
    {
        _Len = compute_diameter<GenericPolytope>
                ::template compute<NT>(P);
        initialize(P, parameters, rng);
    }

    template <typename GenericPolytope>
    Walk(GenericPolytope const& P, thread_params &parameters, RandomNumberGenerator &rng,
         parameters const& params)
    {
        _Len = params.set_L ? params.m_L
                          : compute_diameter<GenericPolytope>
                            ::template compute<NT>(P);
        initialize(P, parameters, rng);
    }

    template
    <
        typename GenericPolytope
    >
    inline void apply(GenericPolytope const& P,
                      thread_params &parameters,
                      unsigned int const& walk_length,
                      RandomNumberGenerator &rng)
    {
        unsigned int n = P.dimension();
        NT T = rng.sample_urdist() * _Len;
        const NT dl = 0.995;

        for (auto j=0u; j<walk_length; ++j)
        {
            T = rng.sample_urdist() * _Len;
            parameters.v = GetDirection<Point>::apply(n, rng);
            parameters.p0 = parameters.p;
            int it = 0;
            while (it < 50*n)
            {
                auto pbpair = P.line_positive_intersect(parameters.p, parameters.v, parameters.lambdas,
                                                        parameters.Av, parameters.lambda_prev);
                if (T <= pbpair.first) {
                    parameters.p += (T * parameters.v);
                    parameters.lambda_prev = T;
                    break;
                }
                parameters.lambda_prev = dl * pbpair.first;
                parameters.p += (parameters.lambda_prev * parameters.v);
                T -= parameters.lambda_prev;
                P.compute_reflection(parameters.v, parameters.p, pbpair.second);
                it++;
            }
            if (it == 50*n) 
            {
                parameters.p = parameters.p0;
            }
        }
    }

    inline void update_delta(NT L)
    {
        _Len = L;
    }

private :

    template
    <
        typename GenericPolytope
    >
    inline void initialize(GenericPolytope const& P,
                           thread_params &parameters,
                           RandomNumberGenerator &rng)
    {
        unsigned int n = P.dimension();
        const NT dl = 0.995;
        parameters.lambdas.setZero(P.num_of_hyperplanes());
        parameters.Av.setZero(P.num_of_hyperplanes());

        parameters.v = GetDirection<Point>::apply(n, rng);
        NT T = rng.sample_urdist() * _Len;
        int it = 0;

        std::pair<NT, int> pbpair
                = P.line_positive_intersect(parameters.p, parameters.v, parameters.lambdas, parameters.Av);
        if (T <= pbpair.first) 
        {
            parameters.p += (T * parameters.v);
            parameters.lambda_prev = T;
            return;
        }
        parameters.lambda_prev = dl * pbpair.first;
        parameters.p += (parameters.lambda_prev * parameters.v);
        T -= parameters.lambda_prev;
        P.compute_reflection(parameters.v, parameters.p, pbpair.second);

        while (it <= 50*n)
        {
            std::pair<NT, int> pbpair
                    = P.line_positive_intersect(parameters.p, parameters.v, parameters.lambdas, 
                                                parameters.Av, parameters.lambda_prev);
            if (T <= pbpair.first) 
            {
                parameters.p += (T * parameters.v);
                parameters.lambda_prev = T;
                break;
            }else if (it == 50*n) {
                parameters.lambda_prev = rng.sample_urdist() * pbpair.first;
                parameters.p += (parameters.lambda_prev * parameters.v);
                break;
            }
            parameters.lambda_prev = dl * pbpair.first;
            parameters.p += (parameters.lambda_prev * parameters.v);
            T -= parameters.lambda_prev;
            P.compute_reflection(parameters.v, parameters.p, pbpair.second);
            it++;
        }
    }

    NT _Len;
};

};


#endif // RANDOM_WALKS_UNIFORM_BILLIARD_WALK_MULTITHREAD_HPP
