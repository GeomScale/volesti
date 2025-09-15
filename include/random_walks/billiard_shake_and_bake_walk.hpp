// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2025 Vissarion Fisikopoulos
// Copyright (c) 2018-2025 Apostolos Chalkis
// Copyright (c) 2025-2025 Iva Janković

// Contributed and/or modified by Iva Janković, as part of Google Summer of Code 2025 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef RANDOM_WALKS_BILLIARD_SHAKE_AND_BAKE_WALK_HPP
#define RANDOM_WALKS_BILLIARD_SHAKE_AND_BAKE_WALK_HPP

#include <Eigen/Eigen>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <iostream>

#include "sampling/sphere.hpp"
#include "convex_bodies/hpolytope.h"
#include "convex_bodies/convex_body.h"
#include <boost/random/exponential_distribution.hpp>

struct BilliardShakeAndBakeWalk
{

    struct update_parameters
    {
        update_parameters()
                :   facet_prev(-1), hit_ball(false), inner_vi_ak(0.0), ball_inner_norm(0.0), moved_dist(0.0) 
        {}
        int facet_prev;
        bool hit_ball;
        double inner_vi_ak;
        double ball_inner_norm;
        double moved_dist; 
    };

    template
    <
            typename Polytope,
            typename RandomNumberGenerator
    >
    struct Walk
    {
        using Point = typename Polytope::PointType;
        using VT = typename Polytope::VT;
        using NT = typename Point::FT;
        using MT = typename Polytope::MT;
        typedef typename Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> DenseMT;
        static constexpr bool SPARSE = std::is_same_v<MT, Eigen::SparseMatrix<NT, Eigen::RowMajor>>;
        using AA_type = std::conditional_t< SPARSE, typename Eigen::SparseMatrix<NT>, DenseMT >;

        static constexpr NT kDefaultEpsilon = NT(1e-10);

        template <typename GenericPolytope>
        Walk(GenericPolytope &P, 
            Point const& p, 
            RandomNumberGenerator &rng,
            int facet_idx, 
            int nr, //upper bound for reflections
            NT eps = kDefaultEpsilon): epsilon_{eps}
        {
            if(!P.is_normalized()) 
            {
                P.normalize();
            }
            _update_parameters = update_parameters();

            // if not given, square root of dimension
            _nr = (nr > 0)? nr : static_cast<unsigned int>(std::ceil(std::sqrt(P.dimension())));

            if constexpr (SPARSE) 
            {
                _AA = (P.get_mat() * P.get_mat().transpose());
            } 
            else 
            {
                _AA.noalias() = (DenseMT)(P.get_mat() * P.get_mat().transpose());
            }
            initialize(P, p, facet_idx, rng);
        }

        NT get_epsilon() const noexcept { return epsilon_; }

        void apply(Polytope& P, unsigned int walk_len, RandomNumberGenerator& rng)
        {
            typename Point::Coeff b;
            NT* b_data;
            if constexpr (SPARSE) 
            {
                b = P.get_vec();
                b_data = b.data();
            }

            for (unsigned int step = 0; step < walk_len; ++step)
            {
                _update_parameters.moved_dist = 0.0;
                // unsigned int r = (_nr == 1) ? 1 : 1 + static_cast<unsigned int>(rng.sample_urdist() * _nr); FOR UNIFORM
                double z = rng.sample_trunc_expdist();
                unsigned int r = static_cast<unsigned int>(std::floor((1.0 - z) * _nr)); // INVERSE EXPONENTIAL

                Point _v  = get_direction(P,rng);       
                auto pbair = P.line_first_positive_intersect(_p, _v,_Ar, _Av, _update_parameters);
                NT _lambda_prev = pbair.first;
                if (!std::isfinite(_lambda_prev) || _lambda_prev <= eps  || pbair.second < 0) 
                {
                    _lambda_prev = NT(0);
                    continue;
                }

                // from here same as accelerated billiard walk 
                if constexpr (SPARSE) 
                {
                    _update_parameters.moved_dist = _lambda_prev;
                    NT* Ar_data = _Ar.data();
                    NT* Av_data = _Av.data();

                    for(int i = 0; i < P.num_of_hyperplanes(); ++i) 
                    {
                        if (i == _update_parameters.facet_prev) 
                        {
                            continue; // Av_[i]=0
                        }

                         _distances_set.vec[i].first = ( *(b_data + i) - (*(Ar_data + i)) ) / (*(Av_data + i));
                    }

                    _distances_set.rebuild(_update_parameters.moved_dist);
                } 
                else 
                {
                    _p += (_lambda_prev * _v);
                }

                _A_row_k = P.get_row(pbair.second);
                _update_parameters.facet_prev = pbair.second;

                for (unsigned int k = 1; k < r; ++k) // from there we do reflections
                {
                    if constexpr (SPARSE)
                    {
                        P.compute_reflection_abw_sparse(_v, _p, _update_parameters);
                    }
                    else
                    {
                        P.compute_reflection(_v, _p, _update_parameters);
                    }


                    if constexpr (SPARSE)
                    {
                        pbair = P.line_positive_intersect(_p, _Ar, _Av, _lambda_prev,_distances_set, _AA,_update_parameters);
                    }
                    else
                    {
                        pbair = P.line_positive_intersect(_p, _v,_Ar, _Av, _lambda_prev,_AA,_update_parameters);

                    }

                    _lambda_prev = pbair.first;
                    if (!std::isfinite(_lambda_prev) || _lambda_prev <= eps  || pbair.second < 0) 
                    {
                        _lambda_prev = NT(0);
                        continue;
                    }
                
                    _update_parameters.moved_dist += _lambda_prev;

                    _p += _lambda_prev * _v;

                    _A_row_k = P.get_row(pbair.second);       
                    _update_parameters.facet_prev = pbair.second;
                }
            }
        }

        const Point& getCurrentPoint() const noexcept { return _p; }
        
    private :

        //From here same as Shake and Bake
        Point get_direction(Polytope& P, RandomNumberGenerator& rng)
        {
            int _dim = P.dimension();
            VT z = GetDirection<Point>::apply(_dim, rng).getCoefficients();
            MT I_cc = - _A_row_k * _A_row_k.transpose();
            I_cc.diagonal() += VT::Ones(_dim);
            NT U = rng.sample_urdist();               
            NT r = std::pow(U, NT(1)/NT(_dim-1)); 
            NT cz = _A_row_k.dot(z);
            VT z_tilde  = I_cc*z;
            z_tilde *= r;
            z_tilde /= std::sqrt(NT(1) - cz*cz);
            
            VT v = z_tilde - std::sqrt(NT(1) - r*r) * _A_row_k;
            return Point(v);
        }

        void initialize(Polytope& P,
                        const Point& boundary_pt,
                        int facet_idx,
                        RandomNumberGenerator& rng)
        {
            int _dim = P.dimension();
            int _m = P.num_of_hyperplanes();
            VT b=P.get_vec();
            int active_facet;

            NT kFacetEps = epsilon_;

            // Checking if boundary point belongs to facet_idx
            _p = boundary_pt;
            VT ai = P.get_row(facet_idx);
            NT dist = std::abs(ai.dot(_p.getCoefficients()) - b.coeff(facet_idx));
            if (dist > kFacetEps)
            {
                active_facet = -1;
                for (int i = 0; i < _m; ++i) 
                {
                    VT ai = P.get_row(i);
                    NT dist = std::abs(ai.dot(_p.getCoefficients()) - b.coeff(i));
                    if (dist < kFacetEps) 
                    {
                        active_facet = i;
                        break;
                    }
                }
                if (active_facet < 0)
                {
                    throw std::runtime_error("Boundary point not on any facet!");
                }
            }
            active_facet = facet_idx;

            //Normal of active facet
            _A_row_k = P.get_row(active_facet);

            //Initializing Ar and Av 
            _Ar.setZero(_m);
            _Av.setZero(_m);
            _lambda_hit = NT(0);
            
            //Calculating first Ar and initializing lambda 
            _Ar.noalias() = P.get_mat() * _p.getCoefficients();
            _lambda_hit = NT(0);

            _A_row_k = P.get_row(active_facet);

            _update_parameters.facet_prev = active_facet;

        }

        Point _p;
        NT _lambda_prev;
        AA_type _AA;
        update_parameters _update_parameters;
        typename Point::Coeff _Ar;
        typename Point::Coeff _Av;
        BoundaryOracleHeap<NT> _distances_set;
        NT epsilon_{kDefaultEpsilon};        
        NT _lambda_hit;
        VT _A_row_k;
        int _nr; 
    };

};

#endif