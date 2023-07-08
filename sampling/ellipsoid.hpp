// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2021 Vissarion Fisikopoulos
// Copyright (c) 2018-2021 Apostolos Chalkis
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

    template <typename Ellipsoid, typename RandomNumberGenerator>
    inline static Point apply(unsigned int const& dim,
                              Ellipsoid const& E,
                              RandomNumberGenerator& rng)
    {
        // Generate a point inside a sphere of radius 1.0
        Point p = GetPointInDsphere<Point>::apply(dim, NT(1.0), rng);

        // transform it to a point inside an ellipsoid
        return Point(E.mult_Lcov(p.getCoefficients()));
    }
};


template <typename Point>
struct GetGaussianDirection
{
    typedef typename Point::FT NT;
    typedef typename Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;

    template <typename Ellipsoid, typename RandomNumberGenerator>
    inline static Point apply(unsigned int const& dim,
                              Ellipsoid const& E,   // ellipsoid representing the Gaussian distribution
                              RandomNumberGenerator &rng)
    {
        // Generate a point inside a sphere of radius 1.0
        Point p = GetDirection<Point>::apply(dim, rng, false);

        // Multiply with cholesky matrix
        VT gaussian_vec = E.mult_Lcov(p.getCoefficients());

        // convert to point
        return Point(gaussian_vec);
    }
};


template <typename Point>
struct GetPointOnDellipsoid
{
    typedef typename Point::FT NT;
    typedef typename Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;

    template <typename Ellipsoid, typename RandomNumberGenerator>
    inline static Point apply(unsigned int const& dim,
                              Ellipsoid const& E,   // ellipsoid representing the Gaussian distribution
                              RandomNumberGenerator &rng)
    {
        // Generate a point inside a sphere of radius 1.0
        Point p = GetDirection<Point>::apply(dim, rng, true);

        // Multiply with cholesky matrix
        VT gaussian_vec = E.mult_Lcov(p.getCoefficients());

        // convert to point
        return Point(gaussian_vec);
    }
};


#endif // ELLIPSOID_HPP