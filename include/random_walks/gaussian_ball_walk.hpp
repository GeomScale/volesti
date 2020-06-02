// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef RANDOM_WALKS_GAUSSIAN_BALL_WALK_HPP
#define RANDOM_WALKS_GAUSSIAN_BALL_WALK_HPP

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
    static inline NT compute_delta(GenericPolytope const& P, NT const& a)
    {
        //return ((P.InnerBall()).second * NT(4)) / NT(P.dimension());
        return NT(4) * (P.InnerBall()).second / std::sqrt(std::max(NT(1), a) * NT(P.dimension()));
    }

    Walk (Polytope const& P, Point const& p, NT const& a,
          RandomNumberGenerator &rng)
    {
        _delta = compute_delta(P, a);
    }

    Walk (Polytope const& P,
          Point const& p,
          NT const& a,
          RandomNumberGenerator &rng,
          parameters const& params)
    {
        _delta = params.set_delta ? params.m_L
                                  : compute_delta(P, a);
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
            if (P.is_in(y) == -1)
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
        _delta = delta;
    }

private :
    NT _delta;
};

};

#endif // RANDOM_WALKS_GAUSSIAN_BALL_WALK_HPP
