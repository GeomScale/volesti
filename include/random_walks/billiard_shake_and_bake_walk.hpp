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
#include "random_walks/accelerated_billiard_walk_utils.hpp"

struct BilliardShakeAndBakeWalk
{
    BilliardShakeAndBakeWalk(double L) : param(L, true) {}
    BilliardShakeAndBakeWalk()        : param(0, false) {}

    struct parameters {
        parameters(double L, bool set) : m_L(L), set_L(set) {}
        double m_L;  bool set_L;
    } param;

    template <typename Polytope, typename RandomNumberGenerator>
    struct Walk
    {
        using Point   = typename Polytope::PointType;
        using NT      = typename Point::FT;
        using VT      = typename Polytope::VT;
        using MT      = typename Polytope::MT;
        using DenseMT = Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic>;
        static constexpr bool SPARSE = std::is_same_v<MT, Eigen::SparseMatrix<NT, Eigen::RowMajor>>;

        struct update_parameters {
            int  facet_prev{-1};
            NT   inner_vi_ak{NT(0)};
            bool hit_ball{false};
            NT   ball_inner_norm{NT(1)};
        } params_;

        static constexpr NT kDefaultEps = NT(1e-10);

        template <typename GenericPolytope>
        Walk(GenericPolytope &P,
             const Point &boundary_pt,
             int          facet_idx,
             RandomNumberGenerator &rng,
             NT           eps = kDefaultEps)
          : P_(P), epsilon_(eps)
        {
            if(!P.is_normalized()) {
                P.normalize();
            }
            _Len = compute_diameter<GenericPolytope>::template compute<NT>(P);
            if constexpr (SPARSE)
                _AA = DenseMT(P.get_mat() * P.get_mat().transpose());
            else
                _AA.noalias() = P.get_mat() * P.get_mat().transpose();
            initialize(boundary_pt, facet_idx, rng);
        }

        template <typename GenericPolytope>
        Walk(GenericPolytope &P,
             const Point &boundary_pt,
             int          facet_idx,
             RandomNumberGenerator &rng,
             NT           eps,
             NT           manual_L)
          : P_(P), epsilon_(eps)
        {
            if (!P_.is_normalized()) P_.normalize();
            if constexpr (SPARSE)
                _AA = DenseMT(P.get_mat() * P.get_mat().transpose());
            else
                _AA.noalias() = P.get_mat() * P.get_mat().transpose();
            _Len = manual_L;
            initialize(boundary_pt, facet_idx, rng);
        }

        NT get_epsilon() const noexcept { return epsilon_; }

        void apply(unsigned walk_len, RandomNumberGenerator& rng)
        {
            for (unsigned step = 0; step < walk_len; ++step) {
                params_.facet_prev = -1;
                _lambda_prev      = NT(0);
                lambda_hit_       = NT(0);

                Point v = get_direction(rng);
                NT    T = -std::log(rng.sample_urdist()) * _Len;  

                Point p0 = p_;         
                int   reflections = 0;

                auto [first_lambda, first_facet] = P_.line_positive_intersect_skip(p_, v, Ar_, Av_, _lambda_prev,_AA, params_,params_.facet_prev);

                if (T <= first_lambda) {
                    p_ = p_ + T * v;
                    continue;
                }

                _lambda_prev = first_lambda;
                p_ += first_lambda * v;
                T  -= first_lambda;

                if constexpr (SPARSE) {
                    P_.compute_reflection_abw_sparse(v, p_, params_);
                } else {
                    P_.compute_reflection(v, p_, params_);
                }

                params_.facet_prev = first_facet;
                reflections++;

                while (T > epsilon_ && reflections < 50 * dim_) 
                {
                    auto [lam, facet] = P_.line_positive_intersect_skip(p_, v, Ar_, Av_, _lambda_prev,_AA, params_,params_.facet_prev);

                    if (!std::isfinite(lam) || lam <= epsilon_ || facet < 0)
                        break;

                    if (T <= lam) {
                        p_ += T * v;
                        _lambda_prev = T;
                        break;
                    }

                    _lambda_prev = lam;
                    p_ += lam * v;
                    T  -= lam;

                    if constexpr (SPARSE) {
                        P_.compute_reflection_abw_sparse(v, p_, params_);
                    } else {
                        P_.compute_reflection(v, p_, params_);
                    }
                    params_.facet_prev = facet;
                    reflections++;
                }

                if (reflections >= 50 * dim_) {
                    p_ = p0;
                }
            }
        }


        const Point& getCurrentPoint() const noexcept { return p_; }
        void update_delta(NT L) { _Len = L; }

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
                        int facet_idx,
                        RandomNumberGenerator& rng) {
            dim_ = P_.dimension();
            m_   = P_.num_of_hyperplanes();
            const VT &b = P_.get_vec();
            NT kFacetEps = epsilon_;

            // Checking if boundary point belongs to facet_idx
            p_ = boundary_pt;
            VT ai = P_.get_row(facet_idx);
            NT dist = std::abs(ai.dot(p_.getCoefficients()) - b.coeff(facet_idx));
            if (dist > kFacetEps)
            {
                facet_idx_ = -1;
                for (int i = 0; i < m_; ++i) {
                    VT ai = P_.get_row(i);
                    NT dist = std::abs(ai.dot(p_.getCoefficients()) - b.coeff(i));
                    if (dist < kFacetEps) {
                        facet_idx_ = i;
                        break;
                    }
                }
                if (facet_idx_ < 0)
                {
                    throw std::runtime_error("Boundary point not on any facet!");
                }
            }
            facet_idx_ = facet_idx;

            A_row_k_ = P_.get_row(facet_idx_);
            Ar_.setZero(m_);
            Av_.setZero(m_);
            Ar_.noalias() = P_.get_mat() * p_.getCoefficients();
            params_.facet_prev  = facet_idx_;
         
        }

        Polytope &P_;
        NT epsilon_{kDefaultEps};
        int dim_{0}, m_{0}, facet_idx_{-1};
        NT _lambda_prev{NT(0)}, lambda_hit_{NT(0)}, _Len{NT(1)};
        Point p_;
        VT Ar_, Av_, A_row_k_;
        DenseMT _AA;
    };
};

#endif 
