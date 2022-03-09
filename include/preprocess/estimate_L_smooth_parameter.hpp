// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

//Contributed and/or modified by Alexandros Manochis, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file


#ifndef ESTIMATE_L_SMOOTH_PARAMETER_HPP
#define ESTIMATE_L_SMOOTH_PARAMETER_HPP


template
<
    typename WalkTypePolicy,
    typename Polytope,
    typename Point,
    typename MT,
    typename VT,
    typename RandomNumberGenerator
>
double svd_on_sample(Polytope &P, Point &p, unsigned int const& num_rounding_steps, MT &V, VT &s, VT &Means,
                   unsigned int const& walk_length, RandomNumberGenerator &rng)
{
    typedef typename WalkTypePolicy::template Walk
            <
                    Polytope,
                    RandomNumberGenerator
            > walk;

    typedef RandomPointGenerator <walk> RandomPointGenerator;
    PushBackWalkPolicy push_back_policy;

    unsigned int N = num_rounding_steps;

    std::list<Point> randPoints;
    MT RetMat(N, P.dimension());
    RandomPointGenerator::apply(P, p, N, walk_length, randPoints,
                                push_back_policy, rng);

    
}


#endif
