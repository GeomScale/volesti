// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2020 Apostolos Chalkis

// Contributed by Huu Phuoc Le as part of Google Summer of Code 2022 program

// Licensed under GNU LGPL.3, see LICENCE file

/// Functions to sample correlation matrices w.r.t. a truncated density

#ifndef VOLESTI_SAMPLING_SAMPLE_CORRELATION_MATRICES_HPP
#define VOLESTI_SAMPLING_SAMPLE_CORRELATION_MATRICES_HPP

#include <sampling/sampling.hpp>

template<typename NT, typename MT>
MT getCoefficientsFromMatrix(const MT& mat) {
    int n = mat.rows();
    int d = n * (n - 1) / 2;
    Eigen::Matrix<NT, Eigen::Dynamic, 1> coeffs(d);
    int k = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < i; ++j) {
            coeffs(k) = mat(i, j);
            ++k;
    	}
    }
    return coeffs;
}

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
    typename MT
>
void uniform_correlation_sampling_MT(   const unsigned int &n,
                                        std::list<MT> &randCorMatrices,
                                        const unsigned int &walkL,
                                        const unsigned int &num_points,
                                        unsigned int const& nburns){
    using PointList = std::list<PointType>;
    PointList randPoints;
    
    CorrelationSpectrahedron_MT<PointType> P(n);
    const unsigned int d = P.dimension();
    PointType startingPoint(n);
    RNGType rng(d);

    uniform_sampling<WalkTypePolicy>(randPoints, P, rng, walkL, num_points, startingPoint, nburns);

    for(const auto&p : randPoints){
        MT final_cor_mat = p.mat + p.mat.transpose() - MT::Identity(n, n);
    	randCorMatrices.push_back(final_cor_mat);
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
void gaussian_correlation_sampling( const unsigned int &n,
                                    PointList &randPoints,
                                    const unsigned int &walkL,
                                    const unsigned int &num_points,
                                    const NT &a,
                                    unsigned int const& nburns = 0){
    CorrelationSpectrahedron<PointType> P(n);
    const unsigned int d = P.dimension();
    PointType startingPoint(d);
    RNGType rng(d);

    gaussian_sampling<WalkTypePolicy>(randPoints, P, rng, walkL, num_points, a, startingPoint, nburns);
}

template
<
    typename WalkTypePolicy,
    typename PointType,
    typename RNGType,
    typename MT,
    typename NT
>
void gaussian_correlation_sampling_MT(  const unsigned int &n,
                                        std::list<MT> &randCorMatrices,
                                        const unsigned int &walkL,
                                        const unsigned int &num_points,
                                        const NT &a,
                                        unsigned int const& nburns = 0){
    using PointList = std::list<PointType>;
    PointList randPoints;
    
    CorrelationSpectrahedron_MT<PointType> P(n);
    const unsigned int d = P.dimension();
    PointType startingPoint(n);
    RNGType rng(d);

    gaussian_sampling<WalkTypePolicy>(randPoints, P, rng, walkL, num_points, a, startingPoint, nburns);

    for(const auto&p : randPoints){
        MT final_cor_mat = p.mat + p.mat.transpose() - MT::Identity(n, n);
        randCorMatrices.push_back(final_cor_mat);
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
    CorrelationSpectrahedron<PointType> P(n);
    const unsigned int d = P.dimension();
    PointType startingPoint(d);
    RNGType rng(d);
    PointType _c(c);

    exponential_sampling<WalkTypePolicy>(randPoints, P, rng, walkL, num_points, _c, T, startingPoint, nburns);
}

template
<
        typename WalkTypePolicy,
        typename PointType,
        typename RNGType,
        typename MT,
        typename NT,
        typename VT
>
void exponential_correlation_sampling_MT(   const unsigned int &n,
                                            std::list<MT> &randCorMatrices,
                                            const unsigned int &walkL,
                                            const unsigned int &num_points,
                                            const VT &c,
                                            const NT &T,
                                            unsigned int const& nburns = 0){
    using PointList = std::list<PointType>;
    PointList randPoints;
    
    CorrelationSpectrahedron_MT<PointType> P(n);
    const unsigned int d = P.dimension();
    PointType startingPoint(n);
    RNGType rng(d);
    PointType _c(c);

    exponential_sampling<WalkTypePolicy>(randPoints, P, rng, walkL, num_points, _c, T, startingPoint, nburns);

    for(const auto&p : randPoints){
        MT final_cor_mat = p.mat + p.mat.transpose() - MT::Identity(n, n);
        randCorMatrices.push_back(final_cor_mat);
    }
}

#endif //VOLESTI_SAMPLING_SAMPLE_CORRELATION_MATRICES_HPP
