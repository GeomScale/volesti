// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

// Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef RANDOM_WALKS_UNIFORM_IMPROVED_BILLIARD_WALK_HPP
#define RANDOM_WALKS_UNIFORM_IMPROVED_BILLIARD_WALK_HPP

#include "convex_bodies/ball.h"
#include "convex_bodies/ballintersectconvex.h"
#include "convex_bodies/hpolytope.h"
#include "convex_bodies/vpolytope.h"
#include "convex_bodies/vpolyintersectvpoly.h"
#include "convex_bodies/zpolytope.h"
#include "convex_bodies/zonoIntersecthpoly.h"


// Billiard walk for uniform distribution

struct ImprovedBilliardWalk
{
    ImprovedBilliardWalk(double L)
            :   param(L, true)
    {}

    ImprovedBilliardWalk()
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

    struct update_parameters
    {
        update_parameters()
                :   facet_prev(0), hit_ball(false), inner_vi_ak(0.0), ball_inner_norm(0.0)
        {}
        int facet_prev;
        bool hit_ball;
        double inner_vi_ak;
        double ball_inner_norm;
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
        typedef typename Polytope::MT MT;
        typedef typename Point::FT NT;
        typedef HPolytope<Point> Hpolytope;
        typedef Zonotope<Point> zonotope;
        typedef ZonoIntersectHPoly <zonotope, Hpolytope> ZonoHPoly;
        typedef Ball<Point> BallType;
        typedef BallIntersectPolytope<Polytope,BallType> BallPolytope;

        Walk(Polytope const& P, Point &p, RandomNumberGenerator &rng) : _update_params()
        {
            _L = compute_diameter<Polytope>::template compute<NT>(P);
            _AA = P.get_AA();
            initialize(P, p, rng);
        }

        Walk(BallPolytope const& P, Point &p, RandomNumberGenerator &rng) : _update_params()
        {
            _L = compute_diameter<BallPolytope>::template compute<NT>(P);
            _AA = P.get_AA();
            initialize(P, p, rng);
        }

        Walk(ZonoHPoly const& P, Point & p, RandomNumberGenerator &rng) : _update_params()
        {
            _L = compute_diameter<ZonoHPoly>::template compute<NT>(P);
            _AA = P.get_AA();
            initialize(P, p, rng);
        }

        Walk(Polytope const& P, Point & p, RandomNumberGenerator &rng, parameters const& params) : _update_params()
        {
            if(params.set_L)
            {
                _L = params.m_L;
            }
            else
            {
                _L = compute_diameter<Polytope>::template compute<NT>(P);
            }
            _AA = P.get_AA();
            initialize(P, p, rng);
        }

        Walk(BallPolytope const& P, Point & p, RandomNumberGenerator &rng, parameters const& params) : _update_params()
        {
            if(params.set_L)
            {
                _L = params.m_L;
            }
            else
            {
                _L = compute_diameter<BallPolytope>::template compute<NT>(P);
            }
            _AA = P.get_AA();
            initialize(P, p, rng);
        }

        Walk(ZonoHPoly const& P, Point & p, RandomNumberGenerator &rng, parameters const& params) : _update_params()
        {
            if(params.set_L)
            {
                _L = params.m_L;
            }
            else
            {
                _L = compute_diameter<ZonoHPoly>::template compute<NT>(P);
            }
            _AA = P.get_AA();
            initialize(P, p, rng);
        }

        Walk (BallType const&, Point &, RandomNumberGenerator &,  parameters &) {}

        Walk (BallType const&, Point &, RandomNumberGenerator &) {}

        template
                <
                        typename GenericPolytope
                >
        inline void apply(GenericPolytope const& P,
                          Point &p,   // a point to start
                          unsigned int const& walk_length,
                          RandomNumberGenerator &rng)
        {
            unsigned int n = P.dimension();
            //NT diameter = P.get_diameter();
            NT T = -std::log(rng.sample_urdist()) * _L;
            const NT dl = 0.995;
            int it;

            for (auto j=0u; j<walk_length; ++j)
            {
                T = rng.sample_urdist() * _L;
                _v = GetDirection<Point>::apply(n, rng);
                Point p0 = _p;

                it = 0;
                std::pair<NT, int> pbpair = P.line_positive_intersect(_p, _v, _lambdas, _Av, _lambda_prev, _update_params);
                if (T <= pbpair.first) {
                    _p += (T * _v);
                    _lambda_prev = T;
                    return;
                }

                _lambda_prev = dl * pbpair.first;
                _p += (_lambda_prev * _v);
                T -= _lambda_prev;
                P.compute_reflection(_v, _p, _update_params);
                it++;

                while (it < 100*n)
                {
                    std::pair<NT, int> pbpair
                            = P.line_positive_intersect(_p, _v, _lambdas, _Av, _lambda_prev, _AA, _update_params);
                    if (T <= pbpair.first) {
                        _p += (T * _v);
                        _lambda_prev = T;
                        break;
                    }
                    _lambda_prev = dl * pbpair.first;
                    _p += (_lambda_prev * _v);
                    T -= _lambda_prev;
                    P.compute_reflection(_v, _p, _update_params);
                    it++;
                }
                if (it == 100*n) _p = p0;
            }
            p = _p;
        }

        inline void update_delta(NT L)
        {
            _L = L;
        }

    private :

        template
                <
                        typename GenericPolytope
                >
        inline void initialize(GenericPolytope const& P,
                               Point &p,
                               RandomNumberGenerator &rng)
        {
            //std::cout<<
            unsigned int n = P.dimension();
            const NT dl = 0.995;
            _lambdas.setZero(P.num_of_hyperplanes());
            _Av.setZero(P.num_of_hyperplanes());
            _p = p;
            _v = GetDirection<Point>::apply(n, rng);

            NT T = -std::log(rng.sample_urdist()) * _L;
            Point p0 = _p;
            int it = 0;

            std::pair<NT, int> pbpair
                    = P.line_first_positive_intersect(_p, _v, _lambdas, _Av, _update_params);
            if (T <= pbpair.first) {
                _p += (T * _v);
                _lambda_prev = T;
                return;
            }
            _lambda_prev = dl * pbpair.first;
            _p += (_lambda_prev * _v);
            T -= _lambda_prev;
            P.compute_reflection(_v, _p, _update_params);

            while (it < 100*n)
            {
                std::pair<NT, int> pbpair
                        = P.line_positive_intersect(_p, _v, _lambdas, _Av, _lambda_prev, _AA, _update_params);
                if (T <= pbpair.first) {
                    _p += (T * _v);
                    _lambda_prev = T;
                    break;
                }
                _lambda_prev = dl * pbpair.first;
                _p += (_lambda_prev * _v);
                T -= _lambda_prev;
                P.compute_reflection(_v, _p, _update_params);
                it++;
            }
            if (it == 100*n) _p = p0;
        }

        double _L;
        Point _p;
        Point _v;
        NT _lambda_prev;
        MT _AA;
        update_parameters _update_params;
        typename Point::Coeff _lambdas;
        typename Point::Coeff _Av;
    };

};


#endif


