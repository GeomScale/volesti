// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

// Gaussian Ball Walk assumes:
// - polytope dimension > 0
// - existence of a non-zero inner ball
// - positive Gaussian scale parameter a
// Violating these assumptions leads to undefined sampling behavior.


#ifndef RANDOM_WALKS_GAUSSIAN_BALL_WALK_HPP
#define RANDOM_WALKS_GAUSSIAN_BALL_WALK_HPP
#include <stdexcept>
#include <cmath>
#include <algorithm>

#include "sampling/sphere.hpp"
#include "random_walks/gaussian_helpers.hpp"

// Ball walk with spherical Gaussian target distribution

struct GaussianBallWalk
{

    GaussianBallWalk(double L)
            :   param(L, true)
    {}

    GaussianBallWalk()
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
    static inline NT compute_delta(GenericPolytope& P, NT const& a)
    {
        const auto dim = P.dimension();
        if (dim <= 0)
        {
            throw std::runtime_error(
                "GaussianBallWalk requires a polytope of positive dimension");
        }

        if (a <= NT(0))
        {
            throw std::runtime_error(
                "GaussianBallWalk requires a strictly positive Gaussian scale parameter");
        }

        const NT radius = (P.InnerBall()).second;
        if (radius <= NT(0))
        {
            throw std::runtime_error(
                "GaussianBallWalk requires a polytope with a non-zero inner ball");
        }

        return (NT(4) * radius) / std::sqrt(a * NT(dim));
    }


    Walk (Polytope& P, Point const& p, NT const& a,
          RandomNumberGenerator &rng)
    {
        _delta = compute_delta(P, a);
    }

    Walk (Polytope& P,
          Point const& p,
          NT const& a,
          RandomNumberGenerator &rng,
          parameters const& params)
    {
        _delta = params.set_delta ? NT(params.m_L) : compute_delta(P, a);
        if (_delta <= NT(0))
        {
            throw std::runtime_error(
                "GaussianBallWalk requires a strictly positive step size");
        }
    }

    template<typename BallPolytope>
    inline void apply(BallPolytope const& P,
                      Point &p,   // a point to start
                      NT const& a_i,
                      unsigned int const& walk_length,
                      RandomNumberGenerator& rng)
    {
        for (auto j = 0u; j < walk_length; ++j)
        {
            Point y = GetPointInDsphere<Point>::apply(P.dimension(),
                                                      _delta,
                                                      rng);
            y += p;
            if (P.is_in(y) == 1)
            {
                NT f_x = eval_exp(p, a_i);
                NT f_y = eval_exp(y, a_i);
                NT rnd = rng.sample_urdist();
                if (rnd <= f_y / f_x) {
                    p = y;
                }
            }
        }
    }

    inline void update_delta(NT delta)
    {
        if (delta <= NT(0))
        {
            throw std::runtime_error(
                "GaussianBallWalk step size must be positive");
        }
        _delta = delta;
    }

private :
    NT _delta;
};

};

#endif // RANDOM_WALKS_GAUSSIAN_BALL_WALK_HPP
