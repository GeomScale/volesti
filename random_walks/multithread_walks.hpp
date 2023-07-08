// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

// Contributed and/or modified by Konstantinos Pallikaris, as part of Google Summer of Code 2021 program.

// Licensed under GNU LGPL.3, see LICENCE file

// THIS IS A TEMPORAL FILE. IT CONTAINS A FEW PROTOTYPES

#ifndef RANDOM_WALKS_MULTITHREAD_WALKS_HPP
#define RANDOM_WALKS_MULTITHREAD_WALKS_HPP

#include "sampling/sphere.hpp"
#include "generators/boost_random_number_generator.hpp"
#include "random_walks/gaussian_helpers.hpp"
#include "random_walks/gaussian_cdhr_walk.hpp"
#include "random_walks/gaussian_rdhr_walk.hpp"


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
        typedef thread_parameters<NT, Point> thread_parameters_;
        //typedef thread_params thread_parameters_;

        template <typename GenericPolytope>
        Walk(GenericPolytope& P, thread_parameters_ &parameters, RandomNumberGenerator& rng)
        {
            initialize(P, parameters, rng);
        }

        template
        <
            typename BallPolytope
        >
        inline void apply(BallPolytope const& P,
                          thread_parameters_ &params, // parameters
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
                                               params.lambdas);
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
                               thread_parameters_ &params, // parameters
                               RandomNumberGenerator& rng)
        {
            params.lambdas.setZero(P.num_of_hyperplanes());
            params.rand_coord = rng.sample_uidist();
            NT kapa = rng.sample_urdist();

            std::pair<NT, NT> bpair = P.line_intersect_coord(params.p, params.rand_coord,
                                                             params.lambdas);
            params.p_prev = params.p;
            params.p.set_coord(params.rand_coord, params.p[params.rand_coord] + bpair.first + kapa
                                      * (bpair.second - bpair.first));
        }

    };

};


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
        typedef thread_parameters<NT, Point> thread_parameters_;
        //typedef thread_params thread_parameters_;

        template <typename GenericPolytope>
        Walk(GenericPolytope& P, thread_parameters_ &parameters, RandomNumberGenerator& rng)
        {
            initialize(P, parameters, rng);
        }

        template
        <
                typename BallPolytope
        >
        inline void apply(BallPolytope const& P,
                          thread_parameters_ &params, // parameters
                          unsigned int const& walk_length,
                          RandomNumberGenerator& rng)
        {
            for (auto j=0u; j<walk_length; ++j)
            {
                params.v = GetDirection<Point>::apply(P.dimension(), rng);
                std::pair<NT, NT> bpair = P.line_intersect(params.p, params.v, params.lambdas, params.Av,
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
                               thread_parameters_ &params, // parameters
                               RandomNumberGenerator& rng)
        {
            params.lambdas.setZero(P.num_of_hyperplanes());
            params.Av.setZero(P.num_of_hyperplanes());

            params.v = GetDirection<Point>::apply(P.dimension(), rng);
            std::pair<NT, NT> bpair = P.line_intersect(params.p, params.v, params.lambdas, params.Av);
            params.lambda_prev = rng.sample_urdist() * (bpair.first - bpair.second) + bpair.second;
            params.p = (params.lambda_prev * params.v) + params.p;
        }

    };

};


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
    typedef thread_parameters<NT, Point> thread_parameters_;
    //typedef thread_params thread_parameters_;

    Walk(Polytope& P,
         thread_parameters_ &params,
         NT const& a_i,
         RandomNumberGenerator &rng)
    {
        initialize(P, params, a_i, rng);
    }

    template <typename parameters>
    Walk(Polytope& P,
         thread_parameters_ &params,
         NT const& a_i,
         RandomNumberGenerator &rng,
         parameters&)
    {
        initialize(P, params, a_i, rng);
    }


    template
    <
        typename BallPolytope
    >
    inline void apply(BallPolytope const& P,
                      thread_parameters_ &params, // parameters
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
                                           params.rand_coord_prev, params.lambdas);
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
                           thread_parameters_ &params, // parameters
                           NT const& a_i,
                           RandomNumberGenerator &rng)
    {
        params.lambdas.setZero(P.num_of_hyperplanes());
        params.rand_coord = rng.sample_uidist();

        std::pair <NT, NT> bpair = P.line_intersect_coord(params.p, params.rand_coord, params.lambdas);
        NT dis = chord_random_point_generator_exp_coord
                    (params.p[params.rand_coord] + bpair.second,
                     params.p[params.rand_coord] + bpair.first,
                     a_i, rng);
        params.p_prev = params.p;
        params.p.set_coord(params.rand_coord, dis);
    }

};

};


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
    typedef thread_parameters<NT, Point> thread_parameters_;
    //typedef thread_params thread_parameters_;

    Walk(Polytope&, thread_parameters_ &, NT const&, RandomNumberGenerator&)
    {}

    template <typename parameters>
    Walk(Polytope&, thread_parameters_ &, NT const&, RandomNumberGenerator&,
         parameters&)
    {}

    template
    <
        typename BallPolytope
    >
    inline void apply(BallPolytope const& P,
                      thread_parameters_ &params, // parameters
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

    BilliardWalk_multithread(double L)
            :   param(L, true)
    {}

    BilliardWalk_multithread()
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
    typedef thread_parameters<NT, Point> thread_parameters_;
    //typedef thread_params thread_parameters_;

    template <typename GenericPolytope>
    Walk(GenericPolytope& P, thread_parameters_ &parameters, RandomNumberGenerator &rng)
    {
        _Len = compute_diameter<GenericPolytope>
                ::template compute<NT>(P);
        initialize(P, parameters, rng);
    }

    template <typename GenericPolytope, typename parameters_>
    Walk(GenericPolytope& P, thread_parameters_ &parameters, RandomNumberGenerator &rng,
         parameters_ const& params)
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
    inline void apply(GenericPolytope& P,
                      thread_parameters_ &parameters,
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
    inline void initialize(GenericPolytope& P,
                           thread_parameters_ &parameters,
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


// coordinate directions hit-and-run walk with uniform target distribution
// Parallel version

struct CDHRWalk_multithread
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
    typedef thread_parameters<NT, Point> thread_parameters_;
    //typedef thread_params thread_parameters_;

    template <typename GenericPolytope>
    Walk(GenericPolytope& P, thread_parameters_ &parameters, RandomNumberGenerator& rng)
    {
        initialize(P, parameters, rng);
    }


    template
    <
        typename BallPolytope
    >
    inline void apply(BallPolytope const& P,
                      thread_parameters_ &params, // parameters
                      unsigned int const& walk_length,
                      RandomNumberGenerator &rng)
    {
        for (auto j=0u; j<walk_length; ++j)
        {
            params.rand_coord_prev = params.rand_coord;
            params.rand_coord = rng.sample_uidist();
            NT kapa = rng.sample_urdist();
            std::pair<NT, NT> bpair = P.line_intersect_coord(params.p,
                                                             params.p_prev,
                                                             params.rand_coord,
                                                             params.rand_coord_prev,
                                                             params.lambdas);
            params.p_prev = params.p;
            params.p.set_coord(params.rand_coord, params.p[params.rand_coord] + bpair.first + kapa
                         * (bpair.second - bpair.first));
        }
    }

private :

    template <typename BallPolytope>
    inline void initialize(BallPolytope const& P,
                           thread_parameters_ &params, // parameters
                           RandomNumberGenerator &rng)
    {
        params.lambdas.setZero(P.num_of_hyperplanes());
        params.rand_coord = rng.sample_uidist();
        NT kapa = rng.sample_urdist();

        std::pair<NT, NT> bpair = P.line_intersect_coord(params.p, params.rand_coord, params.lambdas);
        params.p_prev = params.p;
        params.p.set_coord(params.rand_coord, params.p[params.rand_coord] + bpair.first + kapa
                    * (bpair.second - bpair.first));
    }

};

};


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
    typedef thread_parameters<NT, Point> thread_parameters_;
    //typedef thread_params thread_parameters_;

    template <typename GenericPolytope>
    Walk(GenericPolytope& P, thread_parameters_ &parameters, RandomNumberGenerator& rng)
    {
        initialize(P, parameters, rng);
    }

    template
    <
        typename BallPolytope
    >
    inline void apply(BallPolytope const& P,
                      thread_parameters_ &params, // parameters
                      unsigned int const& walk_length,
                      RandomNumberGenerator& rng)
    {
        for (auto j=0u; j<walk_length; ++j)
        {
            params.v = GetDirection<Point>::apply(P.dimension(), rng);
            std::pair<NT, NT> bpair = P.line_intersect(params.p, params.v, params.lambdas, params.Av,
                                                       params.lambda_prev);
            params.lambda_prev = rng.sample_urdist() * (bpair.first - bpair.second)
                    + bpair.second;
            params.p += (params.lambda_prev * params.v);
        }
    }

private :

    template <typename BallPolytope>
    inline void initialize(BallPolytope const& P,
                           thread_parameters_ &params, // parameters
                           RandomNumberGenerator &rng)
    {
        params.lambdas.setZero(P.num_of_hyperplanes());
        params.Av.setZero(P.num_of_hyperplanes());

        params.v = GetDirection<Point>::apply(P.dimension(), rng);
        std::pair<NT, NT> bpair = P.line_intersect(params.p, params.v, params.lambdas, params.Av);
        params.lambda_prev = rng.sample_urdist() * (bpair.first - bpair.second) + bpair.second;
        params.p += (params.lambda_prev * params.v);
    }

};

};


#endif // RANDOM_WALKS_MULTITHREAD_WALKS_HPP
