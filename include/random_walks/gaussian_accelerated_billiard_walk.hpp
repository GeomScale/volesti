// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis
// Copyright (c) 2021 Vaibhav Thakkar
// Copyright (c) 2024 Luca Perju

// Contributed and/or modified by Vaibhav Thakkar, as part of Google Summer of Code 2021 program.
// Contributed and/or modified by Luca Perju, as part of Google Summer of Code 2024 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef RANDOM_WALKS_GAUSSIAN_ACCELERATED_BILLIARD_WALK_HPP
#define RANDOM_WALKS_GAUSSIAN_ACCELERATED_BILLIARD_WALK_HPP

#include <Eigen/Eigen>
#include <cmath>
#include "convex_bodies/orderpolytope.h"
#include "convex_bodies/ellipsoid.h"
#include "convex_bodies/ballintersectconvex.h"
#include "generators/boost_random_number_generator.hpp"
#include "sampling/ellipsoid.hpp"
#include "random_walks/uniform_billiard_walk.hpp"

#include "random_walks/compute_diameter.hpp"


struct GaussianAcceleratedBilliardWalk
{
    GaussianAcceleratedBilliardWalk(double L)
            :   param(L, true)
    {}

    GaussianAcceleratedBilliardWalk()
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
                :   facet_prev(0), hit_ball(false), inner_vi_ak(0.0), ball_inner_norm(0.0), moved_dist(0.0)
        {}
        int facet_prev;
        bool hit_ball;
        double inner_vi_ak;
        double ball_inner_norm;
        double moved_dist;
    };

    parameters param;


    template
    <
        typename Polytope,
        typename RandomNumberGenerator,
        typename E_type = typename Polytope::DenseMT
    >
    struct Walk
    {
        typedef typename Polytope::PointType Point;
        typedef typename Polytope::MT MT;
        typedef typename Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> DenseMT;
        typedef typename Polytope::VT VT;
        typedef typename Point::FT NT;
        using AA_type = std::conditional_t< std::is_same_v<MT, typename Eigen::SparseMatrix<NT, Eigen::RowMajor>>, typename Eigen::SparseMatrix<NT>, DenseMT >; 
        // AA is sparse colMajor if MT is sparse rowMajor, and Dense otherwise
        using AE_type = std::conditional_t< std::is_same_v<MT, typename Eigen::SparseMatrix<NT, Eigen::RowMajor>> && std::is_base_of<Eigen::SparseMatrixBase<E_type>, E_type>::value, typename Eigen::SparseMatrix<NT, Eigen::RowMajor>, DenseMT >; 
      //                                  (               (                                (                   ))                   (                       (      )        )                                     (                   )          )
        // AE is sparse rowMajor if (MT is sparse rowMajor and E is sparse), and Dense otherwise

        void computeLcov(const E_type E)
        {
            if constexpr (std::is_base_of<Eigen::SparseMatrixBase<E_type>, E_type >::value) {
                Eigen::SimplicialLLT<E_type> lltofE;
                lltofE.compute(E);
                if (lltofE.info() != Eigen::Success) {
                    throw std::runtime_error("First Cholesky decomposition failed for sparse matrix!");
                }
                Eigen::SparseMatrix<NT> I(E.cols(), E.cols());
                I.setIdentity();
                Eigen::SparseMatrix<NT> E_inv = lltofE.solve(I);
                Eigen::SimplicialLLT<E_type> lltofEinv;
                lltofEinv.compute(E_inv);
                if (lltofE.info() != Eigen::Success) {
                    throw std::runtime_error("Second Cholesky decomposition failed for sparse matrix!");
                }
                _L_cov = lltofEinv.matrixL();
            } else {
                Eigen::LLT<E_type> lltOfE(E.llt().solve(E_type::Identity(E.cols(), E.cols()))); // compute the Cholesky decomposition of inv(E)
                if (lltOfE.info() != Eigen::Success) {
                    throw std::runtime_error("Cholesky decomposition failed for dense matrix!");
                }
                _L_cov = lltOfE.matrixL();
            }
        }

        template <typename GenericPolytope>
        Walk(GenericPolytope& P,
             Point const& p,
             E_type const& E,   // covariance matrix representing the Gaussian distribution
             RandomNumberGenerator &rng)
        {
            if(!P.is_normalized()) {
                P.normalize();
            }
            _update_parameters = update_parameters();
            _L = compute_diameter<GenericPolytope>::template compute<NT>(P);
            computeLcov(E);
            _E = E;
            if constexpr (std::is_same<AA_type, Eigen::SparseMatrix<NT>>::value) {
                _AA = (P.get_mat() * P.get_mat().transpose());
            } else {
                _AA.noalias() = (DenseMT)(P.get_mat() * P.get_mat().transpose());
            }
            _rho = 1000 * P.dimension(); // upper bound for the number of reflections (experimental)
            initialize(P, p, rng);
        }

        template <typename GenericPolytope>
        Walk(GenericPolytope& P,
             Point const& p,
             E_type const& E,   // covariance matrix representing the Gaussian distribution
             RandomNumberGenerator &rng,
             parameters const& params)
        {
            if(!P.is_normalized()) {
                P.normalize();
            }
            _update_parameters = update_parameters();
            _L = params.set_L ? params.m_L
                              : compute_diameter<GenericPolytope>
                                ::template compute<NT>(P);
            computeLcov(E);
            _E = E;
            if constexpr (std::is_same<AA_type, Eigen::SparseMatrix<NT>>::value) {
                _AA = (P.get_mat() * P.get_mat().transpose());
            } else {
                _AA.noalias() = (DenseMT)(P.get_mat() * P.get_mat().transpose());
            }
            _rho = 1000 * P.dimension(); // upper bound for the number of reflections (experimental)
            initialize(P, p, rng);
        }

        template <typename GenericPolytope>
        inline void apply(GenericPolytope& P,
                          Point &p,       // a point to return the result
                          unsigned int const& walk_length,
                          RandomNumberGenerator &rng)
        {
            unsigned int n = P.dimension();
            NT T;
            const NT dl = 0.995;
            int it;
            NT coef = 1.0;
            NT vEv;
            typename Point::Coeff b;
            NT* b_data;
            if constexpr (std::is_same<MT, Eigen::SparseMatrix<NT, Eigen::RowMajor>>::value) {
                b = P.get_vec();
                b_data = b.data();
            }

            for (auto j=0u; j<walk_length; ++j)
            {
                T = -std::log(rng.sample_urdist()) * _L;

                _v = GetDirection<Point>::apply(n, rng, false);
                _v = Point(_L_cov.template triangularView<Eigen::Lower>() * _v.getCoefficients());
                coef = 1.0;

                vEv = (_v.getCoefficients().transpose() * _E.template selfadjointView<Eigen::Upper>()).dot(_v.getCoefficients());

                Point p0 = _p;

                it = 0;
                std::pair<NT, int> pbpair = P.line_first_positive_intersect(_p, _v, _lambdas, _Av, _update_parameters);

                if (T <= pbpair.first) {
                    _p += (T * _v);
                    _lambda_prev = T;
                    continue;
                }

                _lambda_prev = dl * pbpair.first;
                if constexpr (std::is_same<MT, Eigen::SparseMatrix<NT, Eigen::RowMajor>>::value) {
                    _update_parameters.moved_dist = _lambda_prev;
                    NT* Ar_data = _lambdas.data();
                    NT* Av_data = _Av.data();
                    for(int i = 0; i < P.num_of_hyperplanes(); ++i) {
                        _distances_set.vec[i].first = ( *(b_data + i) - (*(Ar_data + i)) ) / (*(Av_data + i));
                    }
                    // rebuild the heap with the new values of (b - Ar) / Av
                    _distances_set.rebuild(_update_parameters.moved_dist);
                } else {
                    _p += (_lambda_prev * _v);
                }
                T -= _lambda_prev;
                
                T = T * coef;
                coef = P.compute_reflection(_v, _p, _AE, _AEA, vEv, _update_parameters);
                T = T / coef;

                it++;

                while (it < _rho)
                {
                    std::pair<NT, int> pbpair;
                    if constexpr (std::is_same<MT, Eigen::SparseMatrix<NT, Eigen::RowMajor>>::value) {
                        pbpair = P.line_positive_intersect(_p, _lambdas, _Av, _lambda_prev,
                                                           _distances_set, _AA, _update_parameters);
                    } else {
                        pbpair = P.line_positive_intersect(_p, _v, _lambdas, _Av, _lambda_prev,
                                                           _AA, _update_parameters);
                    }

                    if (T <= pbpair.first) {
                        _p += (T * _v);
                        _lambda_prev = T;
                        break;
                    }
                    _lambda_prev = dl * pbpair.first;
                    if constexpr (std::is_same<MT, Eigen::SparseMatrix<NT, Eigen::RowMajor>>::value) {
                        _update_parameters.moved_dist += _lambda_prev;
                    } else {
                        _p += (_lambda_prev * _v);
                    }
                    T -= _lambda_prev;
                    
                    T = T * coef;
                    coef = P.compute_reflection(_v, _p, _AE, _AEA, vEv, _update_parameters);
                    T = T / coef;

                    it++;
                }
                _p += _update_parameters.moved_dist * _v;
                _update_parameters.moved_dist = 0.0;
                if (it == _rho){
                    _p = p0;
                }
            }

            p = _p;
        }

    private :

        template <typename GenericPolytope>
        inline void initialize(GenericPolytope& P,
                               Point const& p,  // a point to start
                               RandomNumberGenerator &rng)
        {
            unsigned int n = P.dimension();
            const NT dl = 0.995;
            _lambdas.setZero(P.num_of_hyperplanes());
            _Av.setZero(P.num_of_hyperplanes());
            _p = p;
            {
                DenseMT temp_AE;
                temp_AE.noalias() = (DenseMT)(P.get_mat() * _E);
                _AEA = temp_AE.cwiseProduct((DenseMT)P.get_mat()).rowwise().sum();
                if constexpr (std::is_same_v<AE_type, DenseMT>)
                    _AE = temp_AE;
                else
                    _AE = temp_AE.sparseView();
            }
            _distances_set = BoundaryOracleHeap<NT>(P.num_of_hyperplanes());


            _v = GetDirection<Point>::apply(n, rng, false);
            _v = Point(_L_cov.template triangularView<Eigen::Lower>() * _v.getCoefficients());
            

            NT T = -std::log(rng.sample_urdist()) * _L;
            Point p0 = _p;
            int it = 0;
            NT coef = 1.0;
            NT vEv = (_v.getCoefficients().transpose() * _E.template selfadjointView<Eigen::Upper>()).dot(_v.getCoefficients());

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
            
            T = T * coef;
            coef = P.compute_reflection(_v, _p, _AE, _AEA, vEv, _update_parameters);
            T = T / coef;

            while (it <= _rho)
            {
                std::pair<NT, int> pbpair
                        = P.line_positive_intersect(_p, _v, _lambdas, _Av, _lambda_prev, _AA, _update_parameters);

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

                T = T * coef;
                coef = P.compute_reflection(_v, _p, _AE, _AEA, vEv, _update_parameters);
                T = T / coef;

                it++;
            }
        }

        NT _L;
        Point _p;
        Point _v;
        NT _lambda_prev;
        AA_type _AA;
        E_type _L_cov;   // LL' = inv(E)
        AE_type _AE;
        E_type _E;
        VT _AEA;
        unsigned int _rho;
        update_parameters _update_parameters;
        typename Point::Coeff _lambdas;
        typename Point::Coeff _Av;
        BoundaryOracleHeap<NT> _distances_set;
    };

};


#endif
