// VolEsti (volume computation and sampling library)

// Copyright (c) I am not sure what to put here :)

// Contributed and/or modified by Iva Janković, as part of Google Summer of Code 2025 program.

// Licensed under GNU LGPL-3.0; see LICENCE file

#ifndef RANDOM_WALKS_SHAKE_AND_BAKE_WALK_HPP
#define RANDOM_WALKS_SHAKE_AND_BAKE_WALK_HPP

#include <Eigen/Eigen>
#include <cmath>
#include <algorithm>

#include "sampling/sphere.hpp"
#include "preprocess/feasible_point.hpp"
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
        Walk(GenericPolytope&       P,
             RandomNumberGenerator& rng,
             Mode                   m   = Mode::Original,
             NT                     eps = kDefaultEpsilon)
            : P_{P}, mode_{m}, epsilon_{eps}
        {
            initialize(rng);
        }

        void set_epsilon(NT eps) noexcept { epsilon_ = eps; }
        NT   get_epsilon() const noexcept { return epsilon_; }

        void apply(unsigned int walk_len, RandomNumberGenerator& rng)
        {
            const NT eps = epsilon_;

            for (unsigned step = 0; step < walk_len; ++step)
            {

                Point v = GetDirection<Point>::apply(dim_, rng);

                //Switching towards the inside of half-space
                NT dot_k = A_row_k_.dot(v.getCoefficients());
                if (dot_k > NT(0)) { v *= NT(-1); dot_k *= NT(-1); }

                auto [lambda_hit, facet_new] = P_.line_positive_intersect(p_, v, Ar_, Av_);

                if (!std::isfinite(lambda_hit) || lambda_hit <= NT(0) || facet_new < 0)
                    continue;

                Point y = p_ + lambda_hit * v;
                if (!y.getCoefficients().allFinite())
                    continue;

                VT A_row_r = P_.get_facet_normal_vec(facet_new);
                //NT dot_r   = A_row_r.dot(v.getCoefficients());

                /* 3. Running: увек прихватамо */
                if (mode_ == Mode::Running) {
                    p_       = y;
                    facet_idx_ = facet_new;
                    A_row_k_   = A_row_r;
                    Ar_.noalias() -= lambda_hit * Av_;   
                    continue;
                }

                /* 4.  
                NT beta;
                if (mode_ == Mode::Original) {
                    NT den = dot_r - dot_k;
                    if (std::abs(den) < eps) continue;
                    beta = std::clamp(dot_r / den, NT(0), NT(1));
                } 
                else {
                    beta = -dot_k;
                }

                if (beta > NT(0) && beta <= NT(1) &&
                    rng.sample_urdist() < beta)
                {
                    p_         = y;
                    facet_idx_ = facet_new;
                    A_row_k_   = A_row_r;
                    Ar_.noalias() -= lambda_hit * Av_;  
                }*/
                /* ако тачка није прихваћена – Ar_ остаје исти */
            }
        }


        const Point& getCurrentPoint() const noexcept { return p_; }

    private:

        void initialize(RandomNumberGenerator& rng)
        {
            dim_ = P_.dimension();
            m_   = P_.num_of_hyperplanes();          

            // Boundary point vector, residual and facet index
            auto [x_vec, Ar_init, facet_idx] = compute_boundary_point<Point>(P_, rng, epsilon_);

            // Generating usable point
            p_ = Point(dim_);
            for (std::size_t i = 0; i < dim_; ++i)
                p_.set_coord(i, x_vec(i));

            //Caching the residual
            Ar_ = std::move(Ar_init);                
            //Allocating the size
            Av_.resize(m_);                          

            facet_idx_ = facet_idx;

            //Av for active facet
            A_row_k_   = P_.get_facet_normal_vec(facet_idx_);
        }


        Polytope& P_;                

        Mode mode_{Mode::Original};
        NT   epsilon_{kDefaultEpsilon};

        std::size_t dim_{0};
        Point       p_;
        int         facet_idx_{-1};
        VT Ar_;            
        VT Av_;            
        std::size_t m_{0};
        VT A_row_k_;
    };
};

#endif // RANDOM_WALKS_SHAKE_AND_BAKE_WALK_HPP