// VolEsti (volume computation and sampling library)

// Copyright (c) I don't know what to write here :)

// Licensed under GNU LGPL-3.0; see LICENCE file

#ifndef RANDOM_WALKS_SHAKE_AND_BAKE_WALK_HPP
#define RANDOM_WALKS_SHAKE_AND_BAKE_WALK_HPP

#include <Eigen/Eigen>
#include <limits>
#include <cmath>
#include <stdexcept>

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
        using NT = typename Point::FT;
        using VT = typename Polytope::VT;        

        Mode mode_{Mode::Original};

        static void set_epsilon(NT new_eps) noexcept { epsilon_ = new_eps; }
        static NT   get_epsilon()      noexcept      { return epsilon_;   }

        template <typename GenericPolytope>
        Walk(GenericPolytope&       P,
             RandomNumberGenerator& rng,
             Mode                   m = Mode::Original)
            : mode_{m}
        {
            initialize(P, rng);
        }

        template <typename GenericPolytope>
        inline void apply(GenericPolytope const& P,
                          Point&                out_p,
                          unsigned int          walk_len,
                          RandomNumberGenerator& rng)
        {
            const NT eps = epsilon_;

            for (unsigned t = 0; t < walk_len; ++t)
            {
                Point v = GetDirection<Point>::apply(P.dimension(), rng);

                NT dot_k = A_row_k_.dot(v.getCoefficients());
                if (dot_k > NT(0)) {
                    v     *= NT(-1);
                    dot_k = -dot_k;
                }

                const int m_fac = static_cast<int>(P.num_of_hyperplanes());
                VT Ar(m_fac), Av(m_fac);
                struct UP { NT inner_vi_ak; int facet_prev; } params;

                auto res        = P.line_first_positive_intersect(p_, v, Ar, Av, params);
                NT   lambda_hit = res.first;
                int  r          = res.second;

                if (!std::isfinite(lambda_hit) ||
                    lambda_hit <= NT(0) ||
                    r < 0)
                {
                    continue;
                }

                Point y = p_ + lambda_hit * v;

                for (std::size_t j = 0; j < dim_; ++j)
                    if (!std::isfinite(y[j])) goto next_iter;

                {
                    VT A_row_r = A_.row(r).transpose();
                    NT dot_r   = A_row_r.dot(v.getCoefficients());

                    if (mode_ == Mode::Running) {
                        p_       = y;
                        _k       = r;
                        A_row_k_ = A_row_r;
                        goto next_iter;
                    }

                    NT beta;
                    if (mode_ == Mode::Original) {
                        NT den = dot_r - dot_k;
                        if (std::abs(den) < eps) goto next_iter;
                        beta = dot_r / den;
                    } else {
                        beta = -dot_k;
                    }

                    if (beta > NT(0) && beta <= NT(1) &&
                        rng.sample_urdist() < beta)
                    {
                        p_       = y;
                        _k       = r;
                        A_row_k_ = A_row_r;
                    }
                }

            next_iter:
                continue;
            }

            out_p = p_;
        }

        const Point& getCurrentPoint() const noexcept { return p_; }

    private:
        template <typename GenericPolytope>
        void initialize(GenericPolytope const& P, RandomNumberGenerator& rng)
        {
            dim_        = P.dimension();
            num_facets_ = P.num_of_hyperplanes();
            A_          = P.get_mat();
            b_          = P.get_vec();

            VT x_vec = compute_boundary_point(A_, b_, rng);

            Point p0(dim_);
            for (std::size_t i = 0; i < dim_; ++i)
                p0.set_coord(i, x_vec(i));
            p_ = p0;

            _k = -1;
            for (std::size_t i = 0; i < num_facets_; ++i) {
                if (std::abs(A_.row(i).dot(x_vec) - b_(i)) < epsilon_) {
                    _k = int(i);
                    break;
                }
            }
            if (_k < 0)
                throw std::runtime_error("Boundary point is not on any facet");

            A_row_k_ = A_.row(_k).transpose();
        }

        /*--------------------------------------------------------+
        |  Members                                                |
        +--------------------------------------------------------*/
        std::size_t dim_{0}, num_facets_{0};
        Point       p_;
        int         _k{-1};

        VT A_row_k_;
        Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> A_;
        Eigen::Matrix<NT,Eigen::Dynamic,1>              b_;

        static inline NT epsilon_ = NT(1e-10);
    };
};

#endif // RANDOM_WALKS_SHAKE_AND_BAKE_WALK_HPP