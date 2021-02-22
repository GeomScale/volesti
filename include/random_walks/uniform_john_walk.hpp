// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

// Contributed and/or modified by Alexandros Manochis, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef RANDOM_WALKS_JOHN_WALK_HPP
#define RANDOM_WALKS_JOHN_WALK_HPP


#include "ellipsoid_walks/john_walker.h"


// John walk for uniform distribution

struct JohnWalk
{
    JohnWalk(double L)
            :   param(L, true)
    {}

    JohnWalk()
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
        typedef typename Polytope::MT MT;
        typedef typename Polytope::VT VT;
        typedef typename Point::FT NT;

        Walk(Polytope &P, Point &p, RandomNumberGenerator &)
        {
            MT A = P.get_mat();
            VT b = P.get_vec(), _vec_point = VT::Zero(P.dimension()), p0 = p.getCoefficients();
            NT r = P.ComputeInnerBall().second;
            johnw = JohnWalker<NT>(p0, A, b, r);
        }

        Walk(Polytope &P, Point & p, RandomNumberGenerator &, parameters const& params)
        {
            MT A = P.get_mat();
            VT b = P.get_vec(), _vec_point = VT::Zero(P.dimension()), p0 = p.getCoefficients();
            NT r = params.set_L ? params.m_L
                          : P.ComputeInnerBall().second;
            johnw = JohnWalker<NT>(p0, A, b, r);
        }

        template
                <
                        typename GenericPolytope
                >
        inline void apply(GenericPolytope &,
                          Point &p,   // a point to start
                          unsigned int const& walk_length,
                          RandomNumberGenerator &)
        {
            for (auto j=0u; j<walk_length; ++j)
            {
                johnw.doSample(_vec_point);
            }
            p.set_coeffs(_vec_point);
        }

    private :

        JohnWalker<NT> johnw;
        VT _vec_point;
    };

};


#endif
