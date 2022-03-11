// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

//Contributed and/or modified by Alexandros Manochis, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file


#ifndef ESTIMATE_L_SMOOTH_PARAMETER_HPP
#define ESTIMATE_L_SMOOTH_PARAMETER_HPP

#include "random_walks/random_walks.hpp"

template
<
    typename WalkTypePolicy = AcceleratedBilliardWalk,
    typename Polytope,
    typename Point,
    typename NegativeGradientFunctor,
    typename RandomNumberGenerator
>
double estimate_L_smooth(Polytope &P, Point &p, unsigned int const& walk_length, 
                         NegativeGradientFunctor F, RandomNumberGenerator &rng)
{
    typedef typename Point::FT NT;
    typedef typename WalkTypePolicy::template Walk
            <
                Polytope,
                RandomNumberGenerator
            > RandomWalk;

    P.ComputeInnerBall();

    unsigned int d = P.dimension();
    unsigned int rnum = 5 * d;
    std::vector<Point> randPoints(1), vecPoint1, vecPoint2;
    std::vector< std::vector<Point> > listOfPoints;

    RandomWalk walk(P, p, rng);
    for (unsigned int i=0; i<rnum; ++i)
    {
        walk.template apply(P, p, walk_length, rng);
        randPoints[0] = p;
        listOfPoints.push_back(randPoints);
    }
    
    NT L = std::numeric_limits<NT>::lowest(), Ltemp;

    for (auto pit=listOfPoints.begin(); pit!=(listOfPoints.end()-1); ++pit)
    {
        vecPoint1 = *pit;

        for (auto qit=(pit+1); qit!=listOfPoints.end(); ++qit)
        {
            vecPoint2 = *qit;
            Ltemp = (F(1, vecPoint1, 0) - F(1, vecPoint2, 0)).length() / (vecPoint1[0] - vecPoint2[0]).length();

            if (Ltemp > L)
            {
                L = Ltemp;
            }
        }
    }
    return L;
}


#endif
