// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

// Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2019 program.
// Contributed and modified by Huu Phuoc Le as part of Google Summer of Code 2022 program

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef RANDOM_WALKS_UNIFORM_BILLIARD_WALK_HPP
#define RANDOM_WALKS_UNIFORM_BILLIARD_WALK_HPP

#include <Eigen/Eigen>

#include "convex_bodies/ball.h"
#include "convex_bodies/ballintersectconvex.h"
#include "convex_bodies/hpolytope.h"
#include "convex_bodies/spectrahedra/spectrahedron.h"
#include "convex_bodies/correlation_matrices/correlation_spectrahedron.hpp"
#include "convex_bodies/correlation_matrices/correlation_spectrahedron_MT.hpp"
#ifndef DISABLE_LPSOLVE
    #include "convex_bodies/vpolytope.h"
    #include "convex_bodies/vpolyintersectvpoly.h"
    #include "convex_bodies/zpolytope.h"
    #include "convex_bodies/zonoIntersecthpoly.h"
#endif
#include "sampling/sphere.hpp"
#include "random_walks/boundary_cdhr_walk.hpp"
#include "generators/boost_random_number_generator.hpp"
#include "sampling/random_point_generators.hpp"
#include "volume/sampling_policies.hpp"
#include "random_walks/compute_diameter.hpp"

// Billiard walk for uniform distribution

struct BilliardWalk
{
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
    Walk(GenericPolytope &P, Point const& p, RandomNumberGenerator &rng)
    {
        _Len = compute_diameter<GenericPolytope>
                ::template compute<NT>(P);
        initialize(P, p, rng);
    }

    template <typename GenericPolytope>
    Walk(GenericPolytope &P, Point const& p, RandomNumberGenerator &rng,
         parameters const& params)
    {
        _Len = params.set_L ? params.m_L
                          : compute_diameter<GenericPolytope>
                            ::template compute<NT>(P);
        initialize(P, p, rng);
    }

    template
    <
        typename GenericPolytope
    >
    inline void apply(GenericPolytope &P,
                      Point& p,   // a point to start
                      unsigned int const& walk_length,
                      RandomNumberGenerator &rng)
    {
        unsigned int n = P.dimension();
        NT T = rng.sample_urdist() * _Len;
        const NT dl = 0.995;

        for (auto j=0u; j<walk_length; ++j)
        {
            T = rng.sample_urdist() * _Len;
            _v = GetDirection<Point>::apply(n, rng);

            Point p0 = _p;
            int it = 0;
            while (it < 50*n)
            {
                auto pbpair = P.line_positive_intersect(_p, _v, _lambdas,
                                                        _Av, _lambda_prev);

                if (T <= pbpair.first) {
                    _p += (T * _v);
                    _lambda_prev = T;
                    break;
                }

                _lambda_prev = dl * pbpair.first;
                _p += (_lambda_prev * _v);
                T -= _lambda_prev;

                P.compute_reflection(_v, _p, pbpair.second);

                it++;
            }
            if (it == 50*n){
                _p = p0;
            }
        }
        p = _p;
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
    inline void initialize(GenericPolytope &P,
                           Point const& p,
                           RandomNumberGenerator &rng)
    {
        unsigned int n = P.dimension();
        const NT dl = 0.995;
        _lambdas.setZero(P.num_of_hyperplanes());
        _Av.setZero(P.num_of_hyperplanes());
        _p = p;
        _v = GetDirection<Point>::apply(n, rng);

        NT T = rng.sample_urdist() * _Len;
        Point p0 = _p;
        int it = 0;
        std::pair<NT, int> pbpair
                = P.line_positive_intersect(_p, _v, _lambdas, _Av);
        if (T <= pbpair.first) {
            _p += (T * _v);
            _lambda_prev = T;
            return;
        }
        _lambda_prev = dl * pbpair.first;
        _p += (_lambda_prev * _v);
        T -= _lambda_prev;
        P.compute_reflection(_v, _p, pbpair.second);
        while (it <= 50*n)
        {
            std::pair<NT, int> pbpair
                    = P.line_positive_intersect(_p, _v, _lambdas, _Av, _lambda_prev);
            if (T <= pbpair.first) {
                _p += (T * _v);
                _lambda_prev = T;
                break;
            }else if (it == 50*n) {
                _lambda_prev = rng.sample_urdist() * pbpair.first;
                _p += (_lambda_prev * _v);
                break;
            }
            _lambda_prev = dl * pbpair.first;
            _p += (_lambda_prev * _v);
            T -= _lambda_prev;
            P.compute_reflection(_v, _p, pbpair.second);
            it++;
        }
    }

    NT _Len;
    Point _p;
    Point _v;
    NT _lambda_prev;
    typename Point::Coeff _lambdas;
    typename Point::Coeff _Av;
};

};

#endif // RANDOM_WALKS_UNIFORM_BILLIARD_WALK_HPP
