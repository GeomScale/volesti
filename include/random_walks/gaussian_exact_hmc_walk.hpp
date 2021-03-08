// VolEsti (volume computation and sampling library)
// Copyright (c) 2021 Vissarion Fisikopoulos
// Copyright (c) 2021 Apostolos Chalkis


// Licensed under GNU LGPL.3, see LICENCE file

#ifndef RANDOM_WALKS_GAUSSIAN_EXACT_HMC_WALK_HPP
#define RANDOM_WALKS_GAUSSIAN_EXACT_HMC_WALK_HPP

#include "sampling/sphere.hpp"



// Exact HMC for sampling from the spherical Gaussian distribution

struct ExactHMCGaussianWalk
{
    ExactHMCGaussianWalk(double L)
            :   param(L, true)
    {}

    ExactHMCGaussianWalk()
            :   param(0, false)
    {}

    
    struct parameters
    {
        parameters(double L, bool set)
                :   m_L(L), set_L(set)
        {}
        double m_L;
        bool set_L;
    };

    parameters param;


template
<
    typename Polytope,
    typename RandomNumberGenerator
>
struct Walk
{
    typedef typename Polytope::PointType Point;
    typedef typename Point::FT NT;
    typedef typename Polytope::VT VT;

    template <typename GenericPolytope>
    Walk(GenericPolytope const& P, Point const& p, NT const& a_i, RandomNumberGenerator &rng)
    {
        _Len = compute_diameter<GenericPolytope>
                ::template compute<NT>(P);
        _omega = std::sqrt(NT(2) * a_i);
        initialize(P, p, a_i, rng);
    }

    template <typename GenericPolytope>
    Walk(GenericPolytope const& P, Point const& p, NT const& a_i, RandomNumberGenerator &rng,
         parameters const& params)
    {
        _Len = params.set_L ? params.m_L
                          : compute_diameter<GenericPolytope>
                            ::template compute<NT>(P);
        _omega = std::sqrt(NT(2) * a_i);
        initialize(P, p, a_i, rng);
    }

    template
    <
        typename GenericPolytope
    >
    inline void apply(GenericPolytope const& P,
                      Point& p,
                      NT const& a_i,
                      unsigned int const& walk_length,
                      RandomNumberGenerator &rng)
    {
        unsigned int n = P.dimension();
        NT T;
        const NT dl = 0.995;

        for (auto j=0u; j<walk_length; ++j)
        {
            T = rng.sample_urdist() * _Len;
            _v = GetDirection<Point>::apply(n, rng, false);
            Point p0 = _p;
            int it = 0;
            while (it < 100*n)
            {
                auto pbpair = P.trigonometric_positive_intersect(_p, _v, _omega);
                if (T <= pbpair.first) {
                    update_position(_p, _v, T, _omega);
                    break;
                }
                _lambda_prev = dl * pbpair.first;
                T -= _lambda_prev;
                update_position(_p, _v, _lambda_prev, _omega);
                P.compute_reflection(_v, _p, pbpair.second);
                it++;
            }
            if (it == 100*n){
                _p = p0;
            }
        }
        p = _p;
    }

    inline void update_delta(NT L)
    {
        _Len = L;
    }

private :

    template
    <
        typename GenericPolytope
    >
    inline void initialize(GenericPolytope const& P,
                           Point const& p,
                           NT const& a_i,
                           RandomNumberGenerator &rng)
    {
        unsigned int n = P.dimension();
        const NT dl = 0.995;
        _p = p;
        _v = GetDirection<Point>::apply(n, rng, false);

        NT T = rng.sample_urdist() * _Len;
        Point p0 = _p;
        int it = 0;

        while (it <= 100*n)
        {
            auto pbpair
                    = P.trigonometric_positive_intersect(_p, _v, _omega);
            if (T <= pbpair.first) {
                update_position(_p, _v, T, _omega);
                break;
            }else if (it == 100*n) {
                _lambda_prev = rng.sample_urdist() * pbpair.first;
                update_position(_p, _v, _lambda_prev, _omega);
                break;
            }
            _lambda_prev = dl * pbpair.first;
            update_position(_p, _v, _lambda_prev, _omega);
            T -= _lambda_prev;
            P.compute_reflection(_v, _p, pbpair.second);
            it++;
        }
    }

    inline void update_position(Point &p, Point &v, NT const& T, NT const& omega)
    {
        NT C, Phi;
        for (size_t i = 0; i < p.dimension(); i++)
        {
            C = std::sqrt(p[i] * p[i] + (v[i] * v[i]) / (omega * omega));
            Phi = std::atan((-v[i]) / (p[i] * omega));
            if (v[i] < 0.0 && Phi < 0.0) {
                Phi += M_PI;
            } else if (v[i] > 0.0 && Phi > 0.0) {
                Phi -= M_PI;
            }
            p.set_coord(i, C * std::cos(omega * T + Phi));
            v.set_coord(i, -C * omega * std::sin(omega * T + Phi));
        }
        
    }

    NT _Len;
    Point _p;
    Point _v;
    NT _omega;
    NT _lambda_prev;
};

};








#endif // RANDOM_WALKS_GAUSSIAN_HMC_WALK_HPP

