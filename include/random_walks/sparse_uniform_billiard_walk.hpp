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

    template <typename GenericPolytope>
    Walk(GenericPolytope& P, const Point& p, RandomNumberGenerator& rng,
            parameters const& user_params,
            const SparseMT& Hessian)
    {
        _Len = user_params.set_L ?
                user_params.m_L :
                NT(6.0) * std::sqrt(static_cast<double>(P.dimension()));

        auto A = P.get_mat();
        _b = P.get_vec();

        compute_cholesky_and_transformations(Hessian, A);
        _oracle_params.emplace(_L_inv, _A_rounded, _A_rounded_row_norms);

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

    void compute_cholesky_and_transformations(const SparseMT &H, const SparseMT &A)
    {
        Eigen::SimplicialLLT<SparseMT, Eigen::Lower> Chol(H);
        
        _L_inv = Chol.matrixL().transpose();
        
        MT A_dense = A.toDense();
        MT A_transposed = A_dense.transpose();
        MT temp = _L_inv.template triangularView<Eigen::Upper>().solve(A_transposed);
        _A_rounded = temp.transpose();
        
        _A_rounded_row_norms.setZero(_A_rounded.rows());
        NT* A_rounded_row_norms_data = _A_rounded_row_norms.data();
        for (int i = 0; i < _A_rounded.rows(); ++i) {
            NT row_norm = _A_rounded.row(i).norm();
            *A_rounded_row_norms_data = row_norm;
            _A_rounded.row(i) /= row_norm;
            A_rounded_row_norms_data++;
        }
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
                
        _Ar.setZero(_A_rounded.rows());
        _Av.setZero(_A_rounded.rows());
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
    MT _A_rounded;
    VT _A_rounded_row_norms;

    NT _Len;
    Point _p, _v;
    VT _Ar, _Av;
    NT _lambda_prev;

    struct OracleParams {
        const SparseMT& L_inv;
        const MT& A_rounded;
        const VT& row_norms;
        NT inner_vi_ak = NT(0);
        int facet_prev = -1;

        OracleParams(const SparseMT& L, const MT& A, const VT& r)
            : L_inv(L), A_rounded(A), row_norms(r) {}
    };
    std::optional<OracleParams> _oracle_params; 

};
};

#endif // RANDOM_WALKS_SPARSE_BILLIARD_WALK_HPP