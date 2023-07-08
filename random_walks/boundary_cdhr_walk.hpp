// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef RANDOM_WALKS_BOUNDARY_CDHR_WALK_HPP
#define RANDOM_WALKS_BOUNDARY_CDHR_WALK_HPP

#include "sampling/sphere.hpp"

// random directions hit-and-run walk with uniform target distribution
// from boundary

struct BCDHRWalk
{
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
        Walk(GenericPolytope& P, Point const& p, RandomNumberGenerator& rng)
        {
            initialize(P, p, rng);
        }

        template
        <
            typename BallPolytope
        >
        inline void apply(BallPolytope const& P,
                          Point& p1,   // a point to start
                          Point& p2,
                          unsigned int const& walk_length,
                          RandomNumberGenerator& rng)
        {
            std::pair<NT, NT> bpair;
            for (auto j = 0u; j < walk_length; ++j)
            {
                auto rand_coord_prev = _rand_coord;
                _rand_coord = rng.sample_uidist();
                NT kapa = rng.sample_urdist();
                bpair = P.line_intersect_coord(_p,
                                               _p_prev,
                                               _rand_coord,
                                               rand_coord_prev,
                                               _lamdas);
                _p_prev = _p;
                _p.set_coord(_rand_coord, _p[_rand_coord] + bpair.first + kapa
                                          * (bpair.second - bpair.first));
            }
            p1 = _p_prev;
            p2 = _p_prev;
            p1.set_coord(_rand_coord, _p_prev[_rand_coord] + bpair.first);
            p2.set_coord(_rand_coord, _p_prev[_rand_coord] + bpair.second);
        }

    private :

        template <typename GenericBody>
        inline void initialize(GenericBody const& P,
                               Point const& p,
                               RandomNumberGenerator& rng)
        {
            _lamdas.setZero(P.num_of_hyperplanes());
            _rand_coord = rng.sample_uidist();
            NT kapa = rng.sample_urdist();
            _p = p;
            std::pair<NT, NT> bpair = P.line_intersect_coord(_p, _rand_coord,
                                                             _lamdas);
            _p_prev = _p;
            _p.set_coord(_rand_coord, _p[_rand_coord] + bpair.first + kapa
                                      * (bpair.second - bpair.first));
        }

        unsigned int _rand_coord;
        Point _p;
        Point _p_prev;
        typename Point::Coeff _lamdas;
    };

};

#endif // RANDOM_WALKS_BOUNDARY_CDHR_WALK_HPP
