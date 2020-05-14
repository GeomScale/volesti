// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef RANDOM_WALKS_GAUSSIAN_CDHR_WALK_HPP
#define RANDOM_WALKS_GAUSSIAN_CDHR_WALK_HPP

#include "volume/new_basic_sampling_features.hpp"

// Pick a point from the distribution exp(-a_i||x||^2) on the coordinate chord
template
<
    typename NT,
    typename RandomNumberGenerator
>
NT chord_random_point_generator_exp_coord(const NT &l,
                                          const NT &u,
                                          const NT &a_i,
                                          RandomNumberGenerator& rng)
{
    NT r, r_val, fn, dis;
    const NT tol = 0.00000001;
    // pick from 1-dimensional gaussian if enough weight is inside polytope P
    if (a_i > tol && u - l >= 2.0 / std::sqrt(2.0 * a_i)) {
        while (true) {
            r = rng.sample_ndist();//rdist(rng2);
            r = r / std::sqrt(2.0 * a_i);
            if (r >= l && r <= u) {
                break;
            }
        }
        dis = r;

    // select using rejection sampling from a bounding rectangle
    } else {
        NT M = get_max_coord(l, u, a_i);
        while (true) {
            r = rng.sample_urdist();
            dis = (1.0 - r) * l + r * u;
            r_val = M * rng.sample_urdist();
            fn = std::exp(-a_i * dis * dis);
            if (r_val < fn) {
                break;
            }
        }
    }
    return dis;
}


// Coordinate directions hit-and-run walk with spherical Gaussian target distribution

struct GaussianCDHRWalk
{

    struct parameters {};
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

    Walk(Polytope const& P, Point & p, NT const& a_i, RandomNumberGenerator &rng)
    {
        initialize(P, p, a_i, rng);
    }

    Walk(Polytope const& P, Point & p, NT const& a_i, RandomNumberGenerator &rng, parameters &)
    {
        initialize(P, p, a_i, rng);
    }

    template
    <
        typename BallPolytope
    >
    inline void apply(BallPolytope const& P,
                      Point &p,   // a point to start
                      NT const& a_i,
                      unsigned int const& walk_length,
                      RandomNumberGenerator &rng)
    {
        for (auto j=0u; j<walk_length; ++j)
        {
            auto rand_coord_prev = _rand_coord;
            _rand_coord = rng.sample_uidist();
            std::pair <NT, NT> bpair =
                    P.line_intersect_coord(_p, _p_prev, _rand_coord,
                                           rand_coord_prev, _lamdas);
            NT dis = chord_random_point_generator_exp_coord
                        (_p[_rand_coord] + bpair.second,
                         _p[_rand_coord] + bpair.first,
                         a_i,
                         rng);
            _p_prev = _p;
            _p.set_coord(_rand_coord, dis);
            /*
            auto rand_coord_prev = _rand_coord;
            _rand_coord = rng.sample_uidist();
            NT kapa = rng.sample_urdist();
            std::pair<NT, NT> bpair = P.line_intersect_coord(_p,
                                                             _p_prev,
                                                             _rand_coord,
                                                             rand_coord_prev,
                                                             _lamdas);
            _p_prev = _p;
            _p.set_coord(_rand_coord, _p[_rand_coord] + bpair.first + kapa
                         * (bpair.second - bpair.first));
            */
        }
        p = _p;
    }

private :

    template <typename BallPolytope>
    inline void initialize(BallPolytope const& P,
                           Point &p,
                           NT const& a_i,
                           RandomNumberGenerator &rng)
    {
        _lamdas.setZero(P.num_of_hyperplanes());
        _rand_coord = rng.sample_uidist();//(rng2);
        _p = p;
        std::pair <NT, NT> bpair = P.line_intersect_coord(_p, _rand_coord, _lamdas);
        NT dis = chord_random_point_generator_exp_coord
                    (_p[_rand_coord] + bpair.second,
                     _p[_rand_coord] + bpair.first,
                     a_i, rng);
        _p_prev = p;
        _p.set_coord(_rand_coord, dis);
        /*
        _lamdas.setZero(P.num_of_hyperplanes());
        _rand_coord = rng.sample_uidist();
        NT kapa = rng.sample_urdist();
        _p = p;
        std::pair<NT, NT> bpair = P.line_intersect_coord(_p, _rand_coord, _lamdas);
        _p_prev = _p;
        _p.set_coord(_rand_coord, _p[_rand_coord] + bpair.first + kapa
                    * (bpair.second - bpair.first));
        */
    }

    unsigned int _rand_coord;
    Point _p;
    Point _p_prev;
    typename Point::Coeff _lamdas;
};

};

#endif // RANDOM_WALKS_GAUSSIAN_CDHR_WALK_HPP
