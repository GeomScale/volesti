// VolEsti (volume computation and sampling library)

// Copyright (c) I am not sure what to put here :)

// Contributed and/or modified by Iva Janković, as part of Google Summer of Code 2025 program.

// Licensed under GNU LGPL-3.0; see LICENCE file

#ifndef RANDOM_WALKS_SHAKE_AND_BAKE_WALK_HPP
#define RANDOM_WALKS_SHAKE_AND_BAKE_WALK_HPP

#include <Eigen/Eigen>
#include <cmath>
#include <algorithm>
#include <stdexcept> 

#include "sampling/sphere.hpp"
#include "convex_bodies/hpolytope.h"
#include "convex_bodies/correlation_matrices/corre_matrix.hpp"

struct ShakeAndBakeWalk
{
    //enum Mode { Original, Limping, Running };

    template <typename Polytope, typename RandomNumberGenerator>
    struct Walk
    {
        using Point = typename Polytope::PointType;
        using VT = typename Polytope::VT;
        using NT = typename Point::FT;

        struct update_parameters {
            int   facet_prev   = -1;  
        };

        update_parameters params_;

        static constexpr NT kDefaultEpsilon = NT(1e-10);

        template <typename GenericPolytope>
        Walk(GenericPolytope&          P,
            const Point&             boundary_pt, 
            int                      facet_idx,    
            RandomNumberGenerator&   rng,
            //Mode                     m   = Mode::Original,
            NT                       eps = kDefaultEpsilon)
            : P_{P}, /*mode_{m}*/ epsilon_{eps}
        {
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

            VT ck = A_row_k_;
            ck.normalize();

            std::vector<NT> u(dim_);
            for (unsigned int i = 0; i < dim_; ++i) {
                u[i] = rng.sample_ndist();
            }

            NT dot = NT(0);
            for (unsigned int i = 0; i < dim_; ++i) dot += u[i] * ck[i];
            for (unsigned int i = 0; i < dim_; ++i) u[i] -= dot * ck[i];

            NT norm_u = NT(0);
            for (auto &x : u) norm_u += x * x;
            norm_u = std::sqrt(norm_u);
            for (auto &x : u) x /= norm_u;

            NT U = rng.sample_urdist();               
            NT r = std::pow(U, NT(1)/(dim_-1));   

            Point z(dim_);
            NT* zdata = z.pointerToData();
            for (unsigned i = 0; i < dim_; ++i) {
                zdata[i] = u[i] * r;
            }

            NT t = -std::sqrt(NT(1) - r*r);

            Point v(dim_);
            NT* vdata = v.pointerToData();
            for (unsigned i = 0; i < dim_; ++i) {
                vdata[i] = zdata[i] + t * ck[i];
            }

            NT check = NT(0);
            for (unsigned int i = 0; i < dim_; ++i) check += ck[i] * vdata[i];
            if (check > NT(0)) 
            {
                for (unsigned int i = 0; i < dim_; ++i) vdata[i] = -vdata[i];
            }

            return v;  
        }

        void initialize(const Point& boundary_pt,
                        int   facet_idx,
                        RandomNumberGenerator& rng)
        {
            dim_ = P_.dimension();
            m_ = P_.num_of_hyperplanes();
            VT b=P_.get_vec();

            NT kFacetEps = NT(1e-8);

            // Input values 
            p_         = boundary_pt;
            if (facet_idx < 0) //if the facet not given we calculate 
            {
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
            }
            else //we recheck if the facet and point compatible
            {
                VT ai = P_.get_facet_normal_vec(facet_idx);
                NT dist = std::abs(ai.dot(p_.getCoefficients()) - b.coeff(facet_idx));
                if (dist > kFacetEps)
                    throw std::runtime_error("That is not facet index of the boundary point!");
                facet_idx_ = facet_idx;
            }

            //Normal of active facet
            A_row_k_   = P_.get_facet_normal_vec(facet_idx_);

            //Calculating first Ar and initializing Av 
            Ar_.setZero(m_);
            Av_.setZero(m_);
            lambda_hit_ = NT(0);
            
            Ar_ = P_.get_mat() * p_.getCoefficients();
            lambda_hit_ = NT(0);

            A_row_k_   = P_.get_facet_normal_vec(facet_idx_);

            params_.facet_prev  = facet_idx_;

        }

        Polytope& P_;                

        //Mode mode_{Mode::Original};
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