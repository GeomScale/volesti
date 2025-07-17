// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2025 Vissarion Fisikopoulos
// Copyright (c) 2018-2025 Apostolos Chalkis
// Copyright (c) 2025 Vladimir Necula

// Contributed and/or modified by Vladimir Necula, as part of Google Summer of Code 2025 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef RANDOM_WALKS_SPARSE_BILLIARD_WALK_HPP
#define RANDOM_WALKS_SPARSE_BILLIARD_WALK_HPP

#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <optional>
#include "convex_bodies/hpolytope.h"
#include "sampling/sphere.hpp"
#include "generators/boost_random_number_generator.hpp"

struct SparseBilliardWalk {

    SparseBilliardWalk(double L)
            :   param(L, true)
    {}

    SparseBilliardWalk()
            :   param(0, false)
    {}

    struct parameters {
        parameters(double L = 0, bool set = false)
            : m_L(L), set_L(set) 
        {}
        double m_L;
        bool set_L;
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
    typedef typename Point::FT NT;
    typedef typename Point::Coeff VT;
    typedef Eigen::SparseMatrix<NT, Eigen::ColMajor> SparseMT;
    typedef Eigen::SparseMatrix<NT, Eigen::RowMajor> SparseRowMT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;

    SparseRowMT _A_original;
    
    template <typename GenericPolytope>
    Walk(GenericPolytope& P, const Point& p, RandomNumberGenerator& rng,
            parameters const& user_params,
            const SparseMT& Hessian)
    {
        _Len = user_params.set_L ?
                user_params.m_L :
                NT(6.0) * std::sqrt(static_cast<double>(P.dimension()));

        compute_cholesky_and_transformations(Hessian);
        
        _b = P.get_vec();
        _A_original = P.get_mat();
        _oracle_params.emplace(_L_inv, _A_original, _b);


        VT p_original = p.getCoefficients();
        VT p_rounded = _L_inv.transpose().template triangularView<Eigen::Lower>() * p_original; 
        Point p_rounded_point(p_rounded);

        initialize(P, p_rounded_point, rng);
    }

    template <typename GenericPolytope>
    void apply(GenericPolytope& P, 
            Point& p, 
            unsigned int const& walk_length,
            RandomNumberGenerator& rng)
    {
        unsigned int n = P.dimension();
        const NT dl = 0.995;
        
        for (auto j = 0u; j < walk_length; ++j)
        {
            NT T = rng.sample_urdist() * _Len;
            _v = GetDirection<Point>::apply(n, rng);

            Point p0 = _p;
            int it = 0;

            while (it < 50 * n)
            {
                std::pair<NT,int> pbpair;

                if (it == 0) {
                    pbpair = P.sparse_line_positive_intersect(_p, _v, _Ar, _Av, *_oracle_params);
                } else {
                    pbpair = P.sparse_line_positive_intersect(_p, _v, _Ar, _Av, _lambda_prev, *_oracle_params);
                }

                if (T <= pbpair.first) {
                    _p += T * _v;
                    _lambda_prev = T;
                    break;
                }

                _lambda_prev = dl * pbpair.first;
                _p += _lambda_prev * _v;
                T -= _lambda_prev;

                P.sparse_compute_reflection(_v, *_oracle_params);
                it++;
            }

            if (it == 50 * n)
                _p = p0;
        } 
        
        VT p_rounded = _p.getCoefficients();
        VT p_original = _L_inv.transpose().template triangularView<Eigen::Lower>().solve(p_rounded);
        p = Point(p_original);
    }

private:

    void compute_cholesky_and_transformations(const SparseMT &H)
    {
        Eigen::SimplicialLLT<SparseMT, Eigen::Lower> Chol(H);
        _L_inv = Chol.matrixL().transpose();
    }

    template <typename GenericPolytope>
    void initialize(GenericPolytope& P,
                    const Point& p_rounded,
                    RandomNumberGenerator& rng)
    {
        unsigned int n = P.dimension();
        const NT dl = 0.995;
        
        _p = p_rounded;
        _v = GetDirection<Point>::apply(n, rng);
                
        _Ar.setZero(P.num_of_hyperplanes());
        _Av.setZero(P.num_of_hyperplanes()); 
        _lambda_prev = 0;
        
        NT T = rng.sample_urdist() * _Len;
        
        auto pbpair = P.sparse_line_positive_intersect(_p, _v, _Ar, _Av, *_oracle_params);
        
        if (pbpair.second < 0) {
            _p += T * _v;
            _lambda_prev = T;
            return;
        }
        
        if (T <= pbpair.first) {
            _p += (T * _v);
            _lambda_prev = T;
            return;
        }
        
        _lambda_prev = dl * pbpair.first;
        _p += (_lambda_prev * _v);
        T -= _lambda_prev;
        
        P.sparse_compute_reflection(_v, *_oracle_params);
        
        int it = 0;
        while (it <= 50*n)
        {
            auto pbpair2 = P.sparse_line_positive_intersect(_p, _v, _Ar, _Av, _lambda_prev, *_oracle_params);
            
            if (T <= pbpair2.first) {
                _p += (T * _v);
                _lambda_prev = T;
                break;
            } else if (it == 50*n) {
                _lambda_prev = rng.sample_urdist() * pbpair2.first;
                _p += (_lambda_prev * _v);
                break;
            }
            
            _lambda_prev = dl * pbpair2.first;
            _p += (_lambda_prev * _v);
            T -= _lambda_prev;
            
            P.sparse_compute_reflection(_v, *_oracle_params);
            it++; 
        }
    }

    VT _b;
    SparseMT _L_inv;

    NT _Len;
    Point _p, _v;
    VT _Ar, _Av;
    NT _lambda_prev;

    struct OracleParams {
        const SparseMT& L_inv;
        const SparseRowMT& A_original;
        const VT& b_original;

        mutable std::vector<VT> A_rounded_rows;
        mutable std::vector<NT> row_norms;
        mutable std::vector<bool> computed;
        mutable std::vector<NT> b_rounded;

        NT  inner_vi_ak = NT(0);
        int facet_prev = -1;

        OracleParams(const SparseMT& L, const SparseRowMT& A, const VT& b)
            : L_inv(L), A_original(A), b_original(b),
            A_rounded_rows(A.rows()),
            row_norms (A.rows()),
            computed (A.rows(), false),
            b_rounded (A.rows()) 
        {}

        const VT& get_normalized_A_rounded_row(int facet) const {
            if (!computed[facet]) {
                VT A_row_dense = VT::Zero(A_original.cols());
                for (typename SparseRowMT::InnerIterator it(A_original, facet); it; ++it) 
                    A_row_dense[it.col()] = it.value();

                A_rounded_rows[facet] = L_inv.template triangularView<Eigen::Upper>().solve(A_row_dense);

                row_norms[facet] = A_rounded_rows[facet].norm();
                if (row_norms[facet] > NT(1e-12))
                    A_rounded_rows[facet] /= row_norms[facet];

                b_rounded[facet] =
                    (row_norms[facet] > NT(1e-12))
                    ? b_original(facet) / row_norms[facet]
                    : b_original(facet);

                computed[facet] = true;
            }
            return A_rounded_rows[facet];
        }

        NT get_b_rounded(int facet) const {
            if (!computed[facet])  get_normalized_A_rounded_row(facet);
            return b_rounded[facet];
        }
    };
    std::optional<OracleParams> _oracle_params; 

};
};

#endif // RANDOM_WALKS_SPARSE_BILLIARD_WALK_HPP