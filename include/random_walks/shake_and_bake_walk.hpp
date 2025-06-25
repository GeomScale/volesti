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

struct ShakeAndBakeWalk
{
    enum Mode { Original, Limping, Running };

    template <typename Polytope, typename RandomNumberGenerator>
    struct Walk
    {
        using Point = typename Polytope::PointType;
        using VT = typename Polytope::VT;
        using NT = typename Point::FT;

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
            //const NT eps = epsilon_; not needed for Running

            for (unsigned step = 0; step < walk_len; ++step)
            {

                Point v = GetDirection<Point>::apply(dim_, rng);

                //Switching towards the inside of half-space
                NT dot_k = A_row_k_.dot(v.getCoefficients());
                if (dot_k > NT(0)) { v *= NT(-1); dot_k *= NT(-1); }

                int  facet_new;

                std::tie(lambda_hit_, facet_new) = P_.line_positive_intersect_skip(p_, v, Ar_, Av_, lambda_hit_, facet_idx_);           

                if (!std::isfinite(lambda_hit_) || lambda_hit_ <= NT(0) ||facet_new < 0)
                    continue;

                p_ +=(lambda_hit_ * v);
                
                facet_idx_ = facet_new;
                A_row_k_   = P_.get_facet_normal_vec(facet_idx_);


            }
        }


        const Point& getCurrentPoint() const noexcept { return p_; }

    private:

        void initialize(const Point& boundary_pt,
                        int   facet_idx,
                        RandomNumberGenerator& rng)
        {
            dim_ = P_.dimension();
            m_   = P_.num_of_hyperplanes();
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

            Point v = GetDirection<Point>::apply(dim_, rng);

            NT dot_k = A_row_k_.dot(v.getCoefficients());
            if (dot_k > NT(0)) { v *= NT(-1); dot_k *= NT(-1); }

            auto [lambda_hit_, facet_new] = P_.line_positive_intersect(p_, v, Ar_, Av_);
            if (!std::isfinite(lambda_hit_) || lambda_hit_ <= NT(0) || facet_new < 0)
                throw std::runtime_error("Shake-and-Bake init: неуспех првог пресека");

            p_ +=(lambda_hit_ * v);          // new boundary point
            facet_idx_ = facet_new;
            A_row_k_   = P_.get_facet_normal_vec(facet_idx_);


        }

        Polytope& P_;                

        //Mode mode_{Mode::Original};
        NT   epsilon_{kDefaultEpsilon};

        std::size_t dim_{0};
        Point       p_;
        int         facet_idx_{-1};
        VT Ar_;            
        VT Av_;            
        NT lambda_hit_;
        std::size_t m_{0};
        VT A_row_k_;
    };
};

#endif // RANDOM_WALKS_SHAKE_AND_BAKE_WALK_HPP