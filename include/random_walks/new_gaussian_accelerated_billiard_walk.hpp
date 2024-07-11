// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

// Contributed and/or modified by Luca Perju, as part of Google Summer of Code 2024 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef RANDOM_WALKS_NEW_GAUSSIAN_ACCELERATED_BILLIARD_WALK_HPP
#define RANDOM_WALKS_NEW_GAUSSIAN_ACCELERATED_BILLIARD_WALK_HPP

#include <Eigen/Eigen>
#include <cmath>
#include <iomanip>
#include "convex_bodies/orderpolytope.h"
#include "convex_bodies/ellipsoid.h"
#include "convex_bodies/ballintersectconvex.h"
#include "generators/boost_random_number_generator.hpp"
#include "sampling/ellipsoid.hpp"
#include "random_walks/uniform_billiard_walk.hpp"

#include "random_walks/compute_diameter.hpp"


struct GABW
{
    GABW(double L)
            :   param(L, true)
    {}

    GABW()
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

    struct update_parameters
    {
        update_parameters()
                :   facet_prev(0), hit_ball(false), inner_vi_ak(0.0), ball_inner_norm(0.0)
        {}
        int facet_prev;
        bool hit_ball;
        double inner_vi_ak;
        double ball_inner_norm;
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
        typedef typename Polytope::MT MT;
        typedef typename Polytope::VT VT;
        typedef typename Point::FT NT;

        template <typename GenericPolytope>
        Walk(GenericPolytope& P,
             Point const& p,
             MT const& E,   // covariance matrix representing the Gaussian distribution
             RandomNumberGenerator &rng)
        {
            _update_parameters = update_parameters();
            _Len = compute_diameter<GenericPolytope>
                ::template compute<NT>(P);

            _AA.noalias() = P.get_mat() * P.get_mat().transpose();
            _rho = 1000 * P.dimension(); // upper bound for the number of reflections (experimental)
            initialize(P, p, E, rng);
        }

        template <typename GenericPolytope>
        Walk(GenericPolytope& P,
             Point const& p,
             MT const& E,   // covariance matrix representing the Gaussian distribution
             RandomNumberGenerator &rng,
             parameters const& params)
        {
            _update_parameters = update_parameters();
            _Len = params.set_L ? params.m_L
                              : compute_diameter<GenericPolytope>
                                ::template compute<NT>(P);
            
            _AA.noalias() = P.get_mat() * P.get_mat().transpose();
            _rho = 1000 * P.dimension(); // upper bound for the number of reflections (experimental)
            initialize(P, p, E, rng);
        }

        template <typename GenericPolytope>
        inline void apply(GenericPolytope& P,
                          Point &p,       // a point to return the result
                          MT const& E,   // covariance matrix representing the Gaussian distribution
                          unsigned int const& walk_length,
                          RandomNumberGenerator &rng)
        {
            unsigned int n = P.dimension();
            NT T;
            const NT dl = 0.995;
            int it;
            NT coef = 1.0;
            NT vEv;

            for (auto j=0u; j<walk_length; ++j)
            {
                T = -std::log(rng.sample_urdist()) * _Len;


                Eigen::LLT<MT> lltOfE(E.inverse()); // compute the Cholesky decomposition of inv(E)
                if (lltOfE.info() != Eigen::Success) {
                    throw std::runtime_error("Cholesky decomposition failed!");
                }
                _L_cov = lltOfE.matrixL();
                _v = GetDirection<Point>::apply(n, rng, false);
                _v = Point(_L_cov.template triangularView<Eigen::Lower>() * _v.getCoefficients());
                coef = 1.0;

                vEv = (_v.getCoefficients().transpose() * E.template selfadjointView<Eigen::Upper>()).dot(_v.getCoefficients());

                Point p0 = _p;

                it = 0;
                std::pair<NT, int> pbpair;
                if(!was_reset)
                    pbpair = P.line_positive_intersect(_p, _v, _lambdas, _Av,
                                                                      _lambda_prev, _update_parameters);
                else {
                    pbpair = P.line_first_positive_intersect(_p, _v, _lambdas,
                                                                     _Av, _update_parameters);
                    was_reset = false;
                }                                                    

                if (T <= pbpair.first) {
                    _p += (T * _v);
                    _lambda_prev = T;
                    continue;
                }

                _lambda_prev = dl * pbpair.first;
                _p += (_lambda_prev * _v);
                T -= _lambda_prev;
                coef = P.compute_reflection(_v, _p, _AE, _AEA, vEv, _update_parameters);
                it++;


                while (it < _rho)
                {
                    std::pair<NT, int> pbpair
                        = P.line_positive_intersect(_p, _v, _lambdas, _Av, _lambda_prev, _AA, _update_parameters);
                    _Av *= coef;
                    _update_parameters.inner_vi_ak *= coef;
                    pbpair.first /= coef;

                    if (T <= pbpair.first) {
                        _p += (T * _v);
                        _lambda_prev = T;
                        break;
                    }
                    _lambda_prev = dl * pbpair.first;
                    _p += (_lambda_prev * _v);
                    T -= _lambda_prev;
                    coef = P.compute_reflection(_v, _p, _AE, _AEA, vEv, _update_parameters);
                    it++;
                }
                if (it == _rho){
                    _p = p0;
                    was_reset = true;
                }
            }

            p = _p;
        }

    private :

        template <typename GenericPolytope>
        inline void initialize(GenericPolytope& P,
                               Point const& p,  // a point to start
                               MT const& E,   // covariance matrix representing the Gaussian distribution
                               RandomNumberGenerator &rng)
        {
            unsigned int n = P.dimension();
            const NT dl = 0.995;
            was_reset = false;
            _lambdas.setZero(P.num_of_hyperplanes());
            _Av.setZero(P.num_of_hyperplanes());
            _p = p;
            _AE.noalias() = P.get_mat() * E;
            _AEA.resize(P.num_of_hyperplanes());
            for(int i = 0; i < P.num_of_hyperplanes(); ++i)
            {
                _AEA(i) = _AE.row(i).dot(P.get_mat().row(i));
            }

            Eigen::LLT<MT> lltOfE(E.inverse()); // compute the Cholesky decomposition of inv(E)
            if (lltOfE.info() != Eigen::Success) {
                throw std::runtime_error("Cholesky decomposition failed!");
            }
            _L_cov = lltOfE.matrixL();
            _v = GetDirection<Point>::apply(n, rng, false);
            _v = Point(_L_cov.template triangularView<Eigen::Lower>() * _v.getCoefficients());

            NT T = -std::log(rng.sample_urdist()) * _Len;
            Point p0 = _p;
            int it = 0;
            NT coef = 1.0;
            NT vEv = (_v.getCoefficients().transpose() * E.template selfadjointView<Eigen::Upper>()).dot(_v.getCoefficients());

            std::pair<NT, int> pbpair
                    = P.line_first_positive_intersect(_p, _v, _lambdas, _Av, _update_parameters);
            if (T <= pbpair.first) {
                _p += (T * _v);
                _lambda_prev = T;
                return;
            }
            _lambda_prev = dl * pbpair.first;
            _p += (_lambda_prev * _v);
            T -= _lambda_prev;
            coef = P.compute_reflection(_v, _p, _AE, _AEA, vEv, _update_parameters);

            while (it <= _rho)
            {
                std::pair<NT, int> pbpair
                        = P.line_positive_intersect(_p, _v, _lambdas, _Av, _lambda_prev, _AA, _update_parameters);
                _Av *= coef;
                _update_parameters.inner_vi_ak *= coef;
                pbpair.first /= coef;

                if (T <= pbpair.first) {
                    _p += (T * _v);
                    _lambda_prev = T;
                    break;
                } else if (it == _rho) {
                    _lambda_prev = rng.sample_urdist() * pbpair.first;
                    _p += (_lambda_prev * _v);
                    break;
                }
                _lambda_prev = dl * pbpair.first;
                _p += (_lambda_prev * _v);
                T -= _lambda_prev;
                coef = P.compute_reflection(_v, _p, _AE, _AEA, vEv, _update_parameters);
                it++;
            }
        }

        NT _Len;
        Point _p;
        Point _v;
        NT _lambda_prev;
        MT _AA;
        MT _L_cov;   // LL' = inv(E)
        MT _AE;
        VT _AEA;
        unsigned int _rho;
        update_parameters _update_parameters;
        typename Point::Coeff _lambdas;
        typename Point::Coeff _Av;
        bool was_reset;
    };

};


#endif
