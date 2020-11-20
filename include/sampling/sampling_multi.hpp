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
                   unsigned int const& max_num_samples,
                   unsigned int &window,
                   unsigned int &Neff_sampled,
                   unsigned int & total_samples,
                   unsigned int const& num_rounding_steps,
                   MT &TotalRandPoints,
                   bool &complete,
                   const VT &starting_point,
                   unsigned int const& nburns,
                   bool request_rounding)
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
    std::cout<<"[2] num_rounding_steps = "<<num_rounding_steps<<std::endl;
    if (request_rounding) {
        TotalRandPoints.setZero(num_rounding_steps, P.dimension());
    } else {
        TotalRandPoints.setZero(max_num_samples, P.dimension());
    }
    //unsigned int min_skip = 2;
    RandomPointGenerator::apply(P, p, Neff, max_num_samples, walk_len, window, 
                                Neff_sampled, total_samples, TotalRandPoints, nburns, complete, 
                                request_rounding, rng);
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
                   unsigned int const& max_num_samples,
                   unsigned int &window,
                   unsigned int &Neff_sampled,
                   unsigned int & total_samples,
                   unsigned int const& num_rounding_steps,
                   MT &TotalRandPoints,
                   bool &complete,
                   const VT &starting_point,
                   unsigned int const& nburns,
                   bool request_rounding,
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
    std::cout<<"[2] num_rounding_steps = "<<num_rounding_steps<<std::endl;
    if (request_rounding) {
        TotalRandPoints.setZero(num_rounding_steps, P.dimension());
    } else {
        TotalRandPoints.setZero(max_num_samples, P.dimension());
    }
    //unsigned int min_skip = 2;
    RandomPointGenerator::apply(P, p, Neff, max_num_samples, walk_len, window,
                                Neff_sampled, total_samples, TotalRandPoints, nburns, complete,
                                rng, request_rounding, WalkType.param);

                                /*Polytope& P,
                      VT &p,   // a point to start
                      unsigned int const& rnum,
                      unsigned int const& max_num_samples,
                      unsigned int const& walk_length,
                      unsigned int &window,
                      unsigned int &Neff_sampled,
                      MT &TotalRandPoints,
                      unsigned int const& nburns,
                      bool &complete,
                      RandomNumberGenerator &rng,
                      bool request_rounding,
                      Parameters const& parameters*/

}

#endif
