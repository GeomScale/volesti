// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

// Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef RANDOM_WALKS_UNIFORM_BALL_WALK_HPP
#define RANDOM_WALKS_UNIFORM_BALL_WALK_HPP

#include "convex_bodies/ball.h"
#include "convex_bodies/ballintersectconvex.h"
#include "convex_bodies/hpolytope.h"
#include "convex_bodies/vpolytope.h"
#include "convex_bodies/vpolyintersectvpoly.h"
#include "convex_bodies/zpolytope.h"
#include "convex_bodies/zonoIntersecthpoly.h"

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
        typedef Ball<Point> BallType;
        typedef BallIntersectPolytope<Polytope,BallType> BallPolytope;
        typedef HPolytope<Point> Hpolytope;
        typedef Zonotope<Point> zonotope;
        typedef ZonoIntersectHPoly <zonotope, Hpolytope> ZonoHPoly;

        Walk (Polytope const& P, Point&, RandomNumberGenerator&) {
            _delta = ((P.InnerBall()).second * NT(4)) / NT(P.dimension());
        }

        Walk (BallPolytope const& P, Point &, RandomNumberGenerator &) {
            _delta = ((P.InnerBall()).second * NT(4)) / NT(P.dimension());
        }

        Walk (BallType const&, Point &, RandomNumberGenerator &) {}

        Walk(ZonoHPoly const& P, Point & p, RandomNumberGenerator &) {
            _delta = ((P.InnerBall()).second * NT(4)) / NT(P.dimension());
        }

        Walk (Polytope const& P, Point&, RandomNumberGenerator&, parameters const& params) {
            if (params.set_delta) {
                _delta = params.m_L;
            } else {
                _delta = ((P.InnerBall()).second * NT(4)) / NT(P.dimension());
            }
        }
        Walk (BallPolytope const& P, Point &, RandomNumberGenerator &, parameters const& params) {
            if (params.set_delta) {
                _delta = params.m_L;
            } else {
                _delta = ((P.InnerBall()).second * NT(4)) / NT(P.dimension());
            }
        }
        Walk (BallType const&, Point &, RandomNumberGenerator &, parameters &) {}

        Walk(ZonoHPoly const& P, Point & p, RandomNumberGenerator &, parameters const& params) {
            if (params.set_delta) {
                _delta = params.m_L;
            } else {
                _delta = ((P.InnerBall()).second * NT(4)) / NT(P.dimension());
            }
        }

        template<typename BallPolytope>
        inline void apply(BallPolytope const& P,
                          Point &p,   // a point to start
                          unsigned int const& walk_length,
                          RandomNumberGenerator &rng)
        {
            //const NT delta = ((P.InnerBall()).second * NT(4)) / NT(P.dimension());

            for (auto j = 0u; j < walk_length; ++j)
            {
                Point y = GetPointInDsphere<Point>::apply(P.dimension(),
                                                          _delta,
                                                          rng);
                y += p;
                if (P.is_in(y) == -1) p = y;
            }
            //std::cout << "use" << parameters.m_L << std::endl;
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
