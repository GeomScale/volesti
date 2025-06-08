#ifndef RANDOM_WALKS_SHAKE_AND_BAKE_WALK_HPP
#define RANDOM_WALKS_SHAKE_AND_BAKE_WALK_HPP

#include <Eigen/Eigen>
#include <limits>
#include "sampling/sphere.hpp"
#include "convex_bodies/hpolytope.h"

struct SBWalk
{
    enum Mode { Original, Limping, Running };

    template <
        typename Polytope,
        typename RandomNumberGenerator
    >
    struct Walk
    {
        using Point = typename Polytope::PointType;
        using NT    = typename Point::FT;

        Mode mode_{Mode::Original};

        template <typename GenericPolytope>
        Walk(GenericPolytope& P,
             Point const&      p0,
             RandomNumberGenerator& rng,
             Mode              m = Mode::Original)
            : mode_{m}
        {
            initialize(P, p0);
        }

        template <typename GenericPolytope>
        inline void apply(GenericPolytope const&  P,
                          Point&                  out_p,
                          unsigned int const      walk_len,
                          RandomNumberGenerator&  rng)
        {
            for (unsigned t = 0; t < walk_len; ++t)
            {
                Point v = GetDirection<Point>::apply(P.dimension(), rng);

                NT dot_k = A_row_k_.dot(v.getCoefficients());
                if (dot_k > NT(0)) {
                    for (size_t j = 0; j < dim_; ++j) {
                        v.set_coord(j, -v[j]);
                    }
                    dot_k = -dot_k;  
                }

                int m_fac = static_cast<int>(P.num_of_hyperplanes());
                Eigen::Matrix<NT, Eigen::Dynamic, 1> Ar(m_fac), Av(m_fac);
                struct UP { NT inner_vi_ak; int facet_prev; };
                UP params;

                auto result = P.line_first_positive_intersect(p_, v, Ar, Av, params);
                NT lambda_hit = result.first;
                int r = result.second;

                if (mode_ == Mode::Running) {
                    if (lambda_hit <= NT(0) || r < 0) {

                        continue;
                    }

                    for (size_t j = 0; j < dim_; ++j) {
                        p_.set_coord(j, p_[j] + lambda_hit * v[j]);
                    }

                    _k = r;
                    A_row_k_ = A_.row(_k);
                    continue;
                }



                if (lambda_hit <= NT(0) || r < 0) {

                    continue;
                }


                Point y = p_;
                for (size_t j = 0; j < dim_; ++j) {
                    y.set_coord(j, p_[j] + lambda_hit * v[j]);
                }


                Eigen::Matrix<NT,1,Eigen::Dynamic> A_row_r = A_.row(r);
                NT dot_r = A_row_r.dot(v.getCoefficients());

                NT beta;
                if (mode_ == Mode::Original) {
                    NT den = dot_r - dot_k;  
                    if (den <= NT(0)) {
                        continue;  
                    }
                    beta = dot_r / den;
                }
                else {
                    beta = -dot_k;
                }

                if (beta > NT(0) && beta <= NT(1)) {
                    if (rng.sample_urdist() < beta) {
                        p_ = y;
                        _k = r;
                        A_row_k_ = A_.row(_k);
                    }

                }

            }

            out_p = p_;
        }

        const Point& getCurrentPoint() const { return p_; }

    private:
        template <typename GenericPolytope>
        void initialize(GenericPolytope const& P, Point const& p0)
        {
            dim_        = P.dimension();
            num_facets_ = P.num_of_hyperplanes();

            A_ = P.get_mat();
            b_ = P.get_vec();

            NT eps = std::numeric_limits<NT>::epsilon() * NT(1e4);
            _k = -1;
            for (size_t i = 0; i < num_facets_; ++i) {
                if ( std::abs(A_.row(i).dot(p0.getCoefficients()) - b_(i)) < eps ) {
                    _k = int(i);
                    break;
                }
            }


            p_       = p0;
            A_row_k_ = A_.row(_k);
        }

        std::size_t dim_{0}, num_facets_{0};
        Point       p_;
        int         _k{-1};

        Eigen::Matrix<NT,1,Eigen::Dynamic>         A_row_k_;
        Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> A_;
        Eigen::Matrix<NT,Eigen::Dynamic,1>              b_;
    };
};

#endif  // RANDOM_WALKS_SHAKE_AND_BAKE_WALK_HPP
