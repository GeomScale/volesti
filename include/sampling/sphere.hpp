// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

//Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.
// Contributed and modified by Huu Phuoc Le as part of Google Summer of Code 2022 program

// Licensed under GNU LGPL.3, see LICENCE file


#ifndef SAMPLERS_SPHERE_HPP
#define SAMPLERS_SPHERE_HPP

#include "convex_bodies/correlation_matrices/corre_matrix.hpp"

template <typename Point>
struct GetDirection
{
    typedef typename Point::FT NT;

    template <typename RandomNumberGenerator>
    inline static Point apply(unsigned int const& dim,
                              RandomNumberGenerator &rng,
                              bool normalize=true)
    {
        NT normal = NT(0);
        Point p(dim);
        NT* data = p.pointerToData();

        if(normalize)
        {
            for (unsigned int i=0; i<dim; ++i)
            {
                *data = rng.sample_ndist();
                normal += *data * *data;
                data++;
            }

            normal = NT(1)/std::sqrt(normal);
            p *= normal;
        }else
        {
            for (unsigned int i=0; i<dim; ++i)
            {
                *data = rng.sample_ndist();
                data++;
            }
        }
        return p;
    }
};

/// Return a random direction for sampling correlation matrices with matrix PointType
template <typename NT>
struct GetDirection<CorreMatrix<NT>>
{
    template <typename RandomNumberGenerator>
    inline static CorreMatrix<NT> apply(unsigned int const& dim,
                              RandomNumberGenerator &rng,
                              bool normalize=true)
    {
        typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>   MT;
        typedef Eigen::Matrix<NT, Eigen::Dynamic, 1>                VT;

        unsigned int n = std::ceil(std::sqrt(2*dim));
        MT mat = MT::Zero(n,n);
        NT normal = NT(0), coeff;

        int i, j;

        if(normalize)
        {
            for(i = 0; i < n ; ++i)
            {
                for(j = 0; j < i; ++j)
                {
                    coeff = rng.sample_ndist();
                    mat(i,j) = coeff;
                    normal +=  coeff * coeff;
                }
            }
            normal = NT(1)/std::sqrt(normal);
            mat *= normal;
        }else
        {
            for(i = 0; i < n ; ++i)
            {
                for(j = 0; j < i; ++j)
                {
                    coeff = rng.sample_ndist();
                    mat(i,j) = coeff;
                }
            }
        }
        return CorreMatrix<NT>(mat);
    }
};

template <typename Point>
struct GetPointInDsphere
{
    template <typename NT, typename RandomNumberGenerator>
    inline static Point apply(unsigned int const& dim,
                              NT const& radius,
                              RandomNumberGenerator &rng)
    {
        Point p = GetDirection<Point>::apply(dim, rng);
        NT U = rng.sample_urdist();
        U = std::pow(U, NT(1)/(NT(dim)));
        p *= radius * U;
        return p;
    }
};

template <typename Point>
struct GetPointOnDsphere
{
    template <typename NT, typename RandomNumberGenerator>
    inline static Point apply(unsigned int const& dim,
                              NT const& radius,
                              RandomNumberGenerator &rng)
    {
        Point p = GetDirection<Point>::apply(dim, rng);
        if (radius != 0) p *= radius;
        return p;
    }
};



#endif // SPHERE_HPP
