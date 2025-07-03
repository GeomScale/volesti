// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2025 Vissarion Fisikopoulos
// Copyright (c) 2018-2025 Apostolos Chalkis
// Copyright (c) 2025-2025 Iva Janković

// Contributed and/or modified by Iva Janković, as part of Google Summer of Code 2025 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef RANDOM_WALKS_SHAKE_AND_BAKE_WALK_HPP
#define RANDOM_WALKS_SHAKE_AND_BAKE_WALK_HPP

#include <Eigen/Eigen>
#include <Eigen/Dense>
#include <cmath>
#include <algorithm>
#include <stdexcept> 

#include "sampling/sphere.hpp"
#include "convex_bodies/hpolytope.h"
#include "convex_bodies/correlation_matrices/corre_matrix.hpp"

struct ShakeAndBakeWalk
{

    template <typename Polytope, typename RandomNumberGenerator>
    struct Walk
    {
        using Point = typename Polytope::PointType;
        using VT = typename Polytope::VT;
        using NT = typename Point::FT;
        using MT = typename Polytope::MT;

        struct update_parameters 
        {
            int   facet_prev   = -1;  
        };

        update_parameters params_;

        static constexpr NT kDefaultEpsilon = NT(1e-10);

        template <typename GenericPolytope>
        Walk(GenericPolytope&          P,
            const Point&             boundary_pt, 
            int                      facet_idx,    
            RandomNumberGenerator&   rng,
            NT                       eps = kDefaultEpsilon)
            : P_{P}, epsilon_{eps}
        {
            P_.normalize();
            initialize(boundary_pt, facet_idx, rng);
        }

        void set_epsilon(NT eps) noexcept { epsilon_ = eps; }
        NT   get_epsilon() const noexcept { return epsilon_; }

         void apply(unsigned int walk_len, RandomNumberGenerator& rng)
        {
            const NT eps = epsilon_; 

            for (unsigned step = 0; step < walk_len; ++step)
            {
                Point v = get_direction(rng);

                int facet_new;
                std::tie(lambda_hit_, facet_new) = P_.line_positive_intersect_skip(p_, v, Ar_, Av_, lambda_hit_, params_);

                if (!std::isfinite(lambda_hit_) || lambda_hit_ <= eps  || facet_new < 0) 
                {
                    lambda_hit_ = NT(0);
                    continue;
                }

                p_ += lambda_hit_ * v;
                facet_idx_ = facet_new;
                A_row_k_   = P_.get_facet_normal_vec(facet_idx_);
                params_.facet_prev  = facet_idx_;
            }
        }


        const Point& getCurrentPoint() const noexcept { return p_; }

    private:

        Point get_direction(RandomNumberGenerator& rng)
        {
            VT z = GetDirection<Point>::apply(dim_, rng).getCoefficients();
            MT I_cc = - A_row_k_ * A_row_k_.transpose();
            I_cc.diagonal() += VT::Ones(dim_);
            NT U = rng.sample_urdist();               
            NT r = std::pow(U, NT(1)/NT(dim_-1)); 
            NT cz = A_row_k_.dot(z);
            VT z_tilde  = I_cc*z;
            z_tilde *= r;
            z_tilde /= std::sqrt(NT(1) - cz*cz);
            
            VT v = z_tilde - std::sqrt(NT(1) - r*r) * A_row_k_;
            return Point(v);
        }

        void initialize(const Point& boundary_pt,
                        int   facet_idx,
                        RandomNumberGenerator& rng)
        {
            dim_ = P_.dimension();
            m_ = P_.num_of_hyperplanes();
            VT b=P_.get_vec();

            NT kFacetEps = epsilon_;

            // Checking if facet index belongs to the boundary point 
            p_ = boundary_pt;
            VT ai = P_.get_facet_normal_vec(facet_idx);
            NT dist = std::abs(ai.dot(p_.getCoefficients()) - b.coeff(facet_idx));
            if (dist > kFacetEps)
                facet_idx_ = -1;
                for (int i = 0; i < m_; ++i) {
                    VT ai = P_.get_facet_normal_vec(i);
                    NT dist = std::abs(ai.dot(p_.getCoefficients()) - b.coeff(i));
                    if (dist < kFacetEps) {
                        facet_idx_ = i;
                        break;
                    }
                }
                if (facet_idx_ < 0)
                    throw std::runtime_error("Boundary point not on any facet!");
                
            facet_idx_ = facet_idx;

            //Normal of active facet
            A_row_k_   = P_.get_facet_normal_vec(facet_idx_);

            //Calculating first Ar and initializing Av 
            Ar_.setZero(m_);
            Av_.setZero(m_);
            lambda_hit_ = NT(0);
            
            Ar_.noalias() = P_.get_mat() * p_.getCoefficients();
            lambda_hit_ = NT(0);

            A_row_k_   = P_.get_facet_normal_vec(facet_idx_);

            params_.facet_prev  = facet_idx_;

        }

        Polytope& P_;                

        NT   epsilon_{kDefaultEpsilon};

        int dim_{0};
        Point       p_;
        int         facet_idx_{-1};
        VT Ar_;            
        VT Av_;            
        NT lambda_hit_;
        int m_{0};
        VT A_row_k_;
    };
};

#endif // RANDOM_WALKS_SHAKE_AND_BAKE_WALK_HPP