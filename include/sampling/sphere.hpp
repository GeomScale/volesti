// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

//Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.

// Licensed under GNU LGPL.3, see LICENCE file


#ifndef SAMPLERS_SPHERE_HPP
#define SAMPLERS_SPHERE_HPP

#include <Eigen/Eigen>

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

        for (unsigned int i=0; i<dim; ++i)
        {
            *data = rng.sample_ndist();
            normal += *data * *data;
            data++;
        }

        normal = NT(1)/std::sqrt(normal);
        if (normalize) p *= normal;

        return p;
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

template <typename VT>
struct GetDirectionTangentPlane
{
    typedef double NT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;

    template <typename RandomNumberGenerator>
    inline static void apply(VT const& p, VT &v,
                              RandomNumberGenerator &rng)
    {
        unsigned int dim = p.rows();
        //NT normal = NT(0);
        //VT ut(dim);
        NT* data = v.data();

        for (unsigned int i=0; i<dim; ++i)
        {
            *data = rng.sample_ndist();
            //normal += *data * *data;
            data++;
        }

        //normal = NT(1)/std::sqrt(normal);
        //v *= normal;
        //std::cout<<"I-v = "<<((MT::Identity(dim, dim) - p * p.transpose()) * v).transpose()<<std::endl;

        v = ((MT::Identity(dim, dim) - p * p.transpose()) * v).eval(); //optimize it
        //std::cout<<"v = "<<v.transpose()<<"\n"<<std::endl;
        v *= (NT(1) / v.norm());
        //Point q(u);
        //return u;
    }
};


template <typename VT>
struct GetGaussianDirectionTangentPlane
{
    typedef double NT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;

    template <typename RandomNumberGenerator>
    inline static void apply(VT const& p, VT &v, MT const& L_chol, MT const& sigma,
                              RandomNumberGenerator &rng)
    {
        unsigned int dim = p.rows();
        //NT normal = NT(0);
        //VT ut(dim);
        NT* data = v.data();

        for (unsigned int i=0; i<dim; ++i)
        {
            *data = rng.sample_ndist();
            //normal += *data * *data;
            data++;
        }
        v = (L_chol.template triangularView<Eigen::Lower>() * v).eval();
        VT sigma_p = sigma*p;
        NT a = (-p.dot(v)) / (p.dot(sigma_p));
        //normal = NT(1)/std::sqrt(normal);
        v += (a*sigma_p);
        //std::cout<<"I-v = "<<((MT::Identity(dim, dim) - p * p.transpose()) * v).transpose()<<std::endl;

        //v = ((MT::Identity(dim, dim) - p * p.transpose()) * v).eval(); //optimize it
        //std::cout<<"v = "<<v.transpose()<<"\n"<<std::endl;
        v *= (NT(1) / v.norm());
        //Point q(u);
        //return u;
    }
};

#endif // SPHERE_HPP
