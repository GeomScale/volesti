// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis
// Copyright (c) 2021 Vaibhav Thakkar

//Contributed and/or modified by Vaibhav Thakkar, as part of Google Summer of Code 2021 program.


#ifndef SAMPLERS_ELLIPSOID_HPP
#define SAMPLERS_ELLIPSOID_HPP

#include "sphere.hpp"
#include "Eigen/Eigen"


template <typename Point>
struct GetPointInDellipsoid
{
    typedef typename Point::FT NT;

    template <typename VT, typename MT, typename RandomNumberGenerator>
    inline static Point apply(unsigned int const& dim,
                              VT const& eigenvalues_inv_sqrt, // sqrt of inverse of eigenvalues of matrix A in (x'Ax <= 1)
                              MT const& EigenVectors,         // eigenvectors of matrix A in (x'Ax <= 1)
                              RandomNumberGenerator& rng)
    {
        // Generate a point inside a sphere of radius 1.0
        Point p = GetPointInDsphere<Point>::apply(dim, NT(1.0), rng);

        // scale points to the ellipsoid using the eigenvalues
        VT scaled_vec = p.getCoefficients().cwiseProduct(eigenvalues_inv_sqrt);

        // rotate with the eigenvectors
        return Point(EigenVectors * scaled_vec);
    }
};


template <typename Point>
struct GetGaussianDirection
{
    typedef typename Point::FT NT;
    typedef typename Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;

    template <typename MT, typename RandomNumberGenerator>
    inline static Point apply(unsigned int const& dim,
                              MT const& L, // cholesky matrix L of the covariance matrix, LL' = Sigma
                              RandomNumberGenerator &rng,
                              bool normalize=true)
    {
        // Generate a point inside a sphere of radius 1.0
        Point p = GetDirection<Point>::apply(dim ,rng);

        // Multiply with cholesky matrix
        VT gaussian_vec = L.template triangularView<Eigen::Lower>() * p.getCoefficients();
        if (normalize) {
            gaussian_vec.normalize();
        }

        // convert to point
        return Point(gaussian_vec);
    }
};


template <typename Point>
struct GetPointOnDellipsoid
{
    typedef typename Point::FT NT;
    typedef typename Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;

    template <typename MT, typename RandomNumberGenerator>
    inline static Point apply(unsigned int const& dim,
                              MT const& L, // cholesky matrix L of the covariance matrix, LL' = Sigma
                              RandomNumberGenerator &rng)
    {
        return GetGaussianDirection<Point>::apply(dim, L, rng, false);
    }
};


#endif // ELLIPSOID_HPP
