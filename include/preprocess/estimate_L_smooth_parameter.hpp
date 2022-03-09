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
double estimate_L_smooth(Polytope &P, Point &p, unsigned int const& walk_length, NegativeGradientFunctor F, RandomNumberGenerator &rng)
{
    typedef typename Point::FT NT;
    typedef typename WalkTypePolicy::template Walk
            <
                    Polytope,
                    RandomNumberGenerator
            > RandomWalk;

    P.ComputeInnerBall();

    //WalkTypePolicy::parameters params(NT(2) * std::sqrt(NT(P.dimension())) * P.InnerBall().second, true);

    unsigned int d = P.dimension();
    unsigned int rnum = 20 * d;
    std::vector<Point> randPoints(1), vecPoint1, vecPoint2;
    std::vector< std::vector<Point> > listOfPoints;

    RandomWalk walk(P, p, rng);
    for (unsigned int i=0; i<rnum; ++i)
    {
        walk.template apply(P, p, walk_length, rng);
        randPoints[0] = p;
        
        listOfPoints.push_back(randPoints);
        //std::cout<<(listOfPoints[i])[0].getCoefficients().transpose()<<std::endl;
    }
    //std::cout<<"length = "<<listOfPoints.size()<<std::endl;
    NT L = std::numeric_limits<NT>::lowest(), Ltemp;

    for (int i=0; i<rnum-1; i++)
    {
        vecPoint1 = listOfPoints[i];
        //std::cout<<"vecPoint1 = "<<vecPoint1[0].getCoefficients().transpose()<<std::endl;
        for (int j=i+1; j<rnum; j++)
        {
            vecPoint2 = listOfPoints[j];
            //std::cout<< "Fi = " <<F(1, vecPoint1, 0).getCoefficients().transpose()<<std::endl;
            //std::cout<< "Fj = " <<F(1, vecPoint2, 0).getCoefficients().transpose()<<std::endl;
            //std::cout <<" (vecPoint1[0] - vecPoint2[0]).length() = "<<(vecPoint1[0] - vecPoint2[0]).length()<<std::endl;
            Ltemp = (F(1, vecPoint1, 0) - F(1, vecPoint2, 0)).length() / (vecPoint1[0] - vecPoint2[0]).length();
            //std::cout<<"Ltemp = "<<Ltemp<<std::endl;
            if (Ltemp > L)
            {
                L = Ltemp;
            }
        }
    }
    return L;
}


#endif
