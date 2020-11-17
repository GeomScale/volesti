#ifndef SAMPLING_MULTI_H
#define SAMPLING_MULTI_H

#include "random_point_gen_multi.hpp"


template <typename WalkTypePolicy,
        typename Polytope,
        typename RandomNumberGenerator,
        typename MT,
        typename VT
>
void uniform_sampling_speedup(Polytope &P,
                   RandomNumberGenerator &rng,
                   const unsigned int &walk_len,
                   const unsigned int &Neff,
                   unsigned int &window,
                   MT &EssRandPoints,
                   unsigned int &Neff_sampled,
                   MT &TotalRandPoints,
                   bool rounding_requested,
                   bool &complete,
                   const VT &starting_point,
                   unsigned int const& nburns)
{
    typedef typename WalkTypePolicy::template Walk
            <
                    Polytope,
                    RandomNumberGenerator
            > walk;

    //RandomNumberGenerator rng(P.dimension());
    //PushBackWalkPolicy push_back_policy;
    typedef RandomPointEfficientGenerator<walk> RandomPointGenerator;

    VT p = starting_point;

    //RandomPointGenerator::apply(P, p, nburns, walk_len, randPoints,
    //                            push_back_policy, rng, WalkType.param);
    EssRandPoints.setZero(P.dimension(), Neff);
    if (rounding_requested) TotalRandPoints.setZero(20*P.dimension(), P.dimension());
    unsigned int min_skip = 2;
    RandomPointGenerator::apply(P, p, Neff, walk_len, min_skip, window, EssRandPoints, Neff_sampled, TotalRandPoints,
                                nburns, rounding_requested, complete, rng);
}


template <
        typename Polytope,
        typename RandomNumberGenerator,
        typename MT,
        typename VT,
        typename WalkTypePolicy
>
void uniform_sampling_speedup(Polytope &P,
                   RandomNumberGenerator &rng,
                   const unsigned int &walk_len,
                   const unsigned int &Neff,
                   unsigned int &window,
                   MT &EssRandPoints,
                   unsigned int &Neff_sampled,
                   MT &TotalRandPoints,
                   bool rounding_requested,
                   bool &complete,
                   const VT &starting_point,
                   unsigned int const& nburns,
                   WalkTypePolicy &WalkType)
{
        
    typedef typename WalkTypePolicy::template Walk
            <
                    Polytope,
                    RandomNumberGenerator
            > walk;

    //RandomNumberGenerator rng(P.dimension());
    //PushBackWalkPolicy push_back_policy;
    typedef RandomPointEfficientGenerator<walk> RandomPointGenerator;

    VT p = starting_point;

    //RandomPointGenerator::apply(P, p, nburns, walk_len, randPoints,
    //                            push_back_policy, rng, WalkType.param);
    EssRandPoints.setZero(P.dimension(), Neff);
    if (rounding_requested) TotalRandPoints.setZero(20*P.dimension(), P.dimension());
    unsigned int min_skip = 2;
    RandomPointGenerator::apply(P, p, Neff, walk_len, min_skip, window, EssRandPoints, Neff_sampled, TotalRandPoints,
                                nburns, rounding_requested, complete, rng, WalkType.param);

}

#endif
