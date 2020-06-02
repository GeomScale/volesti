// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

// Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef RANDOM_WALKS_UNIFORM_BALL_WALK_HPP
#define RANDOM_WALKS_UNIFORM_BALL_WALK_HPP

#include "generators/boost_random_number_generator.hpp"

// Ball walk with uniform target distribution

struct BallWalk
{
    BallWalk(double L)
        :   param(L, true)
    {}

    BallWalk()
            :   param(0, false)
    {}

    struct parameters
    {
        parameters(double L, bool set)
            :   m_L(L), set_delta(set)
        {}
        double m_L;
        bool set_delta;
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
        Walk(GenericPolytope const& P, Point const& /*p*/,
             RandomNumberGenerator& /*rng*/)
        {
            _delta = compute_delta(P);
        }

        template <typename GenericPolytope>
        Walk(GenericPolytope const& P, Point const& /*p*/,
             RandomNumberGenerator& /*rng*/, parameters const& params)
        {
            _delta = params.set_delta ? params.m_L
                                      : compute_delta(P);
        }

        template <typename GenericPolytope>
        static inline NT compute_delta(GenericPolytope const& P)
        {
            //return ((P.InnerBall()).second * NT(4)) / NT(P.dimension());
            return (NT(4) * (P.InnerBall()).second) / std::sqrt(NT(P.dimension()));
        }

        template<typename BallPolytope>
        inline void apply(BallPolytope const& P,
                          Point &p,   // a point to start
                          unsigned int const& walk_length,
                          RandomNumberGenerator &rng)
        {
            for (auto j = 0u; j < walk_length; ++j)
            {
                Point y = GetPointInDsphere<Point>::apply(P.dimension(),
                                                          _delta,
                                                          rng);
                y += p;
                if (P.is_in(y) == -1) p = y;
            }
        }

        inline void update_delta(NT delta)
        {
            _delta = delta;
        }

    private:
        double _delta;
    };
};

#endif // RANDOM_WALKS_UNIFORM_BALL_WALK_HPP
