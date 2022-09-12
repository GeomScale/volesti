// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2020 Apostolos Chalkis

// Contributed by Huu Phuoc Le as part of Google Summer of Code 2022 program

// Licensed under GNU LGPL.3, see LICENCE file

/// Functions to sample correlation matrices w.r.t. a truncated density

#ifndef VOLESTI_SAMPLING_SAMPLE_CORRELATION_MATRICES_HPP
#define VOLESTI_SAMPLING_SAMPLE_CORRELATION_MATRICES_HPP

#include <sampling/sampling.hpp>

// New implementations for sampling correlation matrices with CorrelationSpectrahedron and CorrelationSpectrahedron_MT

template
<
    typename WalkTypePolicy,
    typename PointType,
    typename RNGType,
    typename PointList
>
void uniform_correlation_sampling(const unsigned int &n,
                                    PointList &randPoints,
                                    const unsigned int &walkL,
                                    const unsigned int &num_points,
                                    unsigned int const& nburns){
    CorrelationSpectrahedron<PointType> P(n);
    const unsigned int d = P.dimension();
    PointType startingPoint(d);
    RNGType rng(d);

    uniform_sampling<WalkTypePolicy>(randPoints, P, rng, walkL, num_points, startingPoint, nburns);
}

template
<
    typename WalkTypePolicy,
    typename PointType,
    typename RNGType,
    typename PointList
>
void uniform_correlation_sampling_MT(   const unsigned int &n,
                                        PointList &randPoints,
                                        const unsigned int &walkL,
                                        const unsigned int &num_points,
                                        unsigned int const& nburns){
    CorrelationSpectrahedron_MT<PointType> P(n);
    const unsigned int d = P.dimension();
    PointType startingPoint(n);
    RNGType rng(d);

    uniform_sampling<WalkTypePolicy>(randPoints, P, rng, walkL, num_points, startingPoint, nburns);
}

template
<
    typename WalkTypePolicy,
    typename PointType,
    typename RNGType,
    typename PointList,
    typename NT
>
void gaussian_correlation_sampling( const unsigned int &n,
                                    PointList &randPoints,
                                    const unsigned int &walkL,
                                    const unsigned int &num_points,
                                    const NT &a,
                                    unsigned int const& nburns = 0){

    typedef CorrelationSpectrahedron<PointType> SpectrahedronType;
    typedef typename WalkTypePolicy::template Walk<SpectrahedronType, RNGType> Walk;

    PushBackWalkPolicy push_back_policy;

    SpectrahedronType P(n);
    const unsigned int d = P.dimension();
    PointType p(d);
    RNGType rng(d);

    Walk walk(P, p, a, rng);

    for (unsigned int i = 0; i < nburns; ++i){
        walk.template apply(P, p, a, walkL, rng);
    }

    for (unsigned int i = 0; i < num_points; ++i){
        walk.template apply(P, p, a, walkL, rng);
        push_back_policy.apply(randPoints, p);
    }
}

template
<
    typename WalkTypePolicy,
    typename PointType,
    typename RNGType,
    typename PointList,
    typename NT
>
void gaussian_correlation_sampling_MT(  const unsigned int &n,
                                        PointList &randPoints,
                                        const unsigned int &walkL,
                                        const unsigned int &num_points,
                                        const NT &a,
                                        unsigned int const& nburns = 0){

    typedef CorrelationSpectrahedron_MT<PointType> SpectrahedronType;
    typedef typename WalkTypePolicy::template Walk<SpectrahedronType, RNGType> Walk;

    PushBackWalkPolicy push_back_policy;

    SpectrahedronType P(n);
    const unsigned int d = P.dimension();
    PointType p(n);
    RNGType rng(d);

    Walk walk(P, p, a, rng);

    for (unsigned int i = 0; i < nburns; ++i){
        walk.template apply(P, p, a, walkL, rng);
    }

    for (unsigned int i = 0; i < num_points; ++i){
        walk.template apply(P, p, a, walkL, rng);
        push_back_policy.apply(randPoints, p);
    }
}

template
<
        typename WalkTypePolicy,
        typename PointType,
        typename RNGType,
        typename PointList,
        typename NT,
        typename VT
>
void exponential_correlation_sampling(  const unsigned int &n,
                                        PointList &randPoints,
                                        const unsigned int &walkL,
                                        const unsigned int &num_points,
                                        const VT &c,
                                        const NT &T,
                                        unsigned int const& nburns = 0){

    typedef CorrelationSpectrahedron<PointType> SpectrahedronType;
    typedef typename WalkTypePolicy::template Walk<SpectrahedronType, RNGType> Walk;

    PushBackWalkPolicy push_back_policy;

    SpectrahedronType P(n);
    const unsigned int d = P.dimension();
    PointType p(d);
    RNGType rng(d);
    PointType _c(c);

    Walk walk(P, p, _c, T, rng);

    for (unsigned int i = 0; i < nburns; ++i){
        walk.template apply(P, p, walkL, rng);
        push_back_policy.apply(randPoints, p);
    }

    for (unsigned int i = 0; i < num_points; ++i){
        walk.template apply(P, p, walkL, rng);
        push_back_policy.apply(randPoints, p);
    }
}

template
<
        typename WalkTypePolicy,
        typename PointType,
        typename RNGType,
        typename PointList,
        typename NT,
        typename VT
>
void exponential_correlation_sampling_MT(   const unsigned int &n,
                                            PointList &randPoints,
                                            const unsigned int &walkL,
                                            const unsigned int &num_points,
                                            const VT &c,
                                            const NT &T,
                                            unsigned int const& nburns = 0){


    typedef CorrelationSpectrahedron<PointType> SpectrahedronType;
    typedef typename WalkTypePolicy::template Walk<SpectrahedronType, RNGType> Walk;

    PushBackWalkPolicy push_back_policy;

    SpectrahedronType P(n);
    const unsigned int d = P.dimension();
    PointType p(n);
    RNGType rng(d);
    PointType _c(c);

    Walk walk(P, p, _c, T, rng);

    for (unsigned int i = 0; i < nburns; ++i){
        walk.template apply(P, p, walkL, rng);
        push_back_policy.apply(randPoints, p);
    }

    for (unsigned int i = 0; i < num_points; ++i){
        walk.template apply(P, p, walkL, rng);
        push_back_policy.apply(randPoints, p);
    }
}

#endif //VOLESTI_SAMPLING_SAMPLE_CORRELATION_MATRICES_HPP
