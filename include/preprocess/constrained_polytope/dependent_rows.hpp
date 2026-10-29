// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2026 Vissarion Fisikopoulos
// Copyright (c) 2018-2026 Apostolos Chalkis
// Copyright (c) 2026      Dimitrios Pavlou

// Contributed and/or modified by Dimitrios Pavlou, as part of Google Summer of Code 2026 program

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef CONSTRAINED_POLYTOPE_DEPENDENT_ROWS_HPP
#define CONSTRAINED_POLYTOPE_DEPENDENT_ROWS_HPP

#include <algorithm>
#include <cmath>
#include <vector>
#include <stdexcept>
#include <cstdint>
#include <Eigen/Sparse>
#include "cholmod.h"
#include "SuiteSparseQR.hpp"
#include "convex_bodies/constrained_polytope.hpp"

// Configuration for removing linearly dependent equality rows.
struct DependentRowsRemovalConfig {
    // The rank tolerance of the QR factorization.
    double rank_tol = 1e-7;

    // The largest violation a dropped row may have.
    double residual_tol = 1e-7;
};

// Result of removing linearly dependent equality rows.
template <typename Point>
struct DependentRowsRemovalResult {
    // The simplified polytope.
    ConstrainedPolytope<Point> polytope;

    // The rank of A_eq, i.e. the number of rows kept.
    unsigned rank = 0;

    // The number of rows removed.
    unsigned removed_rows = 0;

    // The largest violation of a dropped row observed.
    double max_residual = 0.0;

    // False if a dropped row is not implied by the kept rows.
    bool valid = false;
};

namespace dependent_rows_removal_util {
    // Converts an Eigen sparse matrix into a CHOLMOD sparse matrix.
    // @tparam SparseMT the Eigen sparse matrix type
    // @param Ain the matrix to convert
    // @param cc the CHOLMOD workspace
    // @return a CHOLMOD sparse matrix
    template <typename SparseMT>
    cholmod_sparse* eigen_to_cholmod(SparseMT const& Ain, cholmod_common* cc) 
    {
        Eigen::SparseMatrix<double, Eigen::ColMajor> Mc = Ain.template cast<double>();
        Mc.makeCompressed();

        cholmod_sparse* Aout = cholmod_l_allocate_sparse(
            (size_t)Mc.rows(), (size_t)Mc.cols(), (size_t)Mc.nonZeros(),
            /*sorted*/ true, /*packed*/ true, /*stype*/ 0, CHOLMOD_REAL, cc);
        
        if (!Aout) throw std::runtime_error("eigen_to_cholmod: cholmod_l_allocate_sparse failed");

        int64_t* Ap = (int64_t*)Aout->p;
        int64_t* Ai = (int64_t*)Aout->i;
        double* Ax = (double*)Aout->x;

        for (Eigen::Index k = 0; k <= Mc.outerSize(); ++k) {
            Ap[k] = Mc.outerIndexPtr()[k];
        }
        for (Eigen::Index k = 0; k < Mc.nonZeros(); ++k) {
            Ai[k] = Mc.innerIndexPtr()[k];
            Ax[k] = Mc.valuePtr()[k];
        } 

        return Aout;
    }

    // Returns the indices of a maximal set of linearly independent rows of A, sorted.
    //
    // A rank revealing QR of A^T with column pivoting is computed. The first rank pivot
    // columns of A^T are independent columns of A^T.
    // @param A the matrix
    // @param tol the rank tolerance
    // @param cc the CHOLMOD workspace
    // @return the sorted indices of the independent rows
    inline std::vector<unsigned> independent_rows(Eigen::SparseMatrix<double> const& A,
                                                  double tol,
                                                  cholmod_common* cc)
    {
        cholmod_sparse* At = eigen_to_cholmod(A.transpose(), cc);

        cholmod_sparse* R = nullptr;
        int64_t* E = nullptr;
        int64_t const rank = SuiteSparseQR<double, int64_t>(SPQR_ORDERING_DEFAULT, tol, (int64_t)0,
                                                            At, &R, &E, cc);

        std::vector<unsigned> keep;
        if (rank >= 0) {
            keep.resize((std::size_t)rank);
            for (int64_t k = 0; k < rank; ++k)
                keep[k] = E ? (unsigned)E[k] : (unsigned)k;
            std::sort(keep.begin(), keep.end());
        }

        if (E) cholmod_l_free((size_t)A.rows(), sizeof(int64_t), E, cc);
        cholmod_l_free_sparse(&R, cc);
        cholmod_l_free_sparse(&At, cc);

        if (rank < 0) throw std::runtime_error("independent_rows: SPQR failed");
        return keep;
    }

    // Returns the minimum norm solution of A x = b.
    // @param A the matrix, full row rank
    // @param b the rhs
    // @param tol the rank tolerance
    // @param cc the CHOLMOD workspace
    // @return the solution x
    inline Eigen::VectorXd min_norm_solution(Eigen::SparseMatrix<double> const& A,
                                             Eigen::VectorXd const& b,
                                             double tol,
                                             cholmod_common* cc)
    {
        cholmod_sparse* S = eigen_to_cholmod(A, cc);
        cholmod_dense* B = cholmod_l_allocate_dense((size_t)b.size(), 1, (size_t)b.size(), CHOLMOD_REAL, cc);

        if (!B) {
            cholmod_l_free_sparse(&S, cc);
            throw std::runtime_error("min_norm_solution: cholmod_l_allocate_dense failed");
        }

        std::copy(b.data(), b.data()+b.size(), (double*)B->x);

        cholmod_dense* X = SuiteSparseQR_min2norm<double>(SPQR_ORDERING_DEFAULT, tol, S, B, cc);

        Eigen::VectorXd x;
        if (X) {
            x = Eigen::Map<Eigen::VectorXd>((double*)X->x, A.cols());
            cholmod_l_free_dense(&X, cc);
        }
        cholmod_l_free_dense(&B, cc);
        cholmod_l_free_sparse(&S, cc);

        if (x.size() == 0) throw std::runtime_error("min_norm_solution: SPQR failed");
        return x;
    }

    // Returns an orthonormal basis of the nullspace A.
    // @param A the matrix
    // @param tol the rank tolerance
    // @param cc the CHOLMOD workspace
    // @param rank receives the rank of A
    // @return the d x (d-rank) nullspace basis
    inline Eigen::MatrixXd nullspace_basis(Eigen::SparseMatrix<double> const& A,
                                        double tol,
                                        cholmod_common* cc,
                                        unsigned& rank)
    {
        int64_t const d = A.cols();
        cholmod_sparse* At = eigen_to_cholmod(A.transpose(), cc);

        SuiteSparseQR_factorization<double>* QR = SuiteSparseQR_factorize<double>(
            SPQR_ORDERING_DEFAULT, tol, At, cc);

        if (!QR) {
            cholmod_l_free_sparse(&At, cc);
            throw std::runtime_error("nullspace_basis: SPQR factorization failed");
        }

        int64_t const r = QR->rank;
        int64_t const k = d-r;
        rank = (unsigned)r;

        Eigen::MatrixXd N(d, k);

        if (k > 0) {
            cholmod_dense* X = cholmod_l_allocate_dense((size_t)d, (size_t)k, (size_t)d,
                                                        CHOLMOD_REAL, cc);

            if (!X) {
                SuiteSparseQR_free<double>(&QR, cc);
                cholmod_l_free_sparse(&At, cc);
                throw std::runtime_error("nullspace_basis: cholmod_l_allocate_dense failed");
            }

            double* Xx = (double*)X->x;
            std::fill(Xx, Xx+d*k, 0.0);
            for (int64_t c = 0; c < k; ++c)
                Xx[c*d+(r+c)] = 1.0;

            cholmod_dense* Y = SuiteSparseQR_qmult<double>(SPQR_QX, QR, X, cc);
            cholmod_l_free_dense(&X, cc);

            if (!Y) {
                SuiteSparseQR_free<double>(&QR, cc);
                cholmod_l_free_sparse(&At, cc);
                throw std::runtime_error("nullspace_basis: SPQR qmult failed");
            }

            std::copy((double*)Y->x, (double*)Y->x+d*k, N.data());
            cholmod_l_free_dense(&Y, cc);
        }

        SuiteSparseQR_free<double>(&QR, cc);
        cholmod_l_free_sparse(&At, cc);
        return N;
    }
};

// Removes the linearly dependent rows of A_eq x = b_eq.
//
// A maximal independent set of rows is kept, checking that the
// dropped rows are implied by the rest within residual_tol.
// @tparam Point the point type
// @param P the polytope
// @param config the configuration
// @return the result, valid is false if a dropped row is not implied
template <typename Point>
DependentRowsRemovalResult<Point>
remove_dependent_rows(ConstrainedPolytope<Point> const& P,
                      DependentRowsRemovalConfig const& config = DependentRowsRemovalConfig{})
{
    using VT = typename ConstrainedPolytope<Point>::VT;
    using MT = typename ConstrainedPolytope<Point>::MT;
    using Triplet = typename ConstrainedPolytope<Point>::Triplet;

    DependentRowsRemovalResult<Point> res;

    unsigned const d = P.getDimension();
    unsigned const m = P.getNumEqualities();
    MT const& A_eq = P.getEqualities();
    VT const& b_eq = P.getEqualityRHS();

    if (m == 0) {
        res.polytope = P;
        res.valid = true;
        return res;
    }

    Eigen::SparseMatrix<double> A = A_eq.template cast<double>();
    Eigen::VectorXd b = b_eq.template cast<double>();

    cholmod_common cc;
    cholmod_l_start(&cc);
    std::vector<unsigned> keep;
    MT A_keep;
    VT b_keep;
    Eigen::VectorXd x0;

    try {
        keep = dependent_rows_removal_util::independent_rows(A, config.rank_tol, &cc);

        std::vector<Triplet> triplets;
        b_keep.resize((Eigen::Index)keep.size());

        for (std::size_t r = 0; r < keep.size(); ++r) {
            for (typename MT::InnerIterator it(A_eq, keep[r]); it; ++it)
                triplets.emplace_back((int)r, (int)it.col(), it.value());
            b_keep((Eigen::Index)r) = b_eq(keep[r]);
        }

        A_keep.resize((Eigen::Index)keep.size(), d);
        A_keep.setFromTriplets(triplets.begin(), triplets.end());

        // Any x solves an empty system.
        if (keep.empty()) {
            x0 = Eigen::VectorXd::Zero(d);
        } else {
            x0 = dependent_rows_removal_util::min_norm_solution(
                A_keep.template cast<double>(), b_keep.template cast<double>(),
                config.rank_tol, &cc);
        }
        cholmod_l_finish(&cc);
    } catch (...) {
        cholmod_l_finish(&cc);
        throw;
    }

    std::vector<bool> kept(m, false);
    for (unsigned i : keep) kept[i] = true;

    Eigen::VectorXd r = A*x0-b;
    double worst = 0.0;
    for (unsigned i = 0; i < m; ++i)
        if (!kept[i]) worst = std::max(worst, std::abs(r(i)));

    res.rank = (unsigned)keep.size();
    res.removed_rows = m-res.rank;
    res.max_residual = worst;

    if (worst > config.residual_tol) return res;

    res.polytope = ConstrainedPolytope<Point>(d, A_keep, b_keep, P.getInequalities(),
                                              P.getInequalityRHS(), P.getLowerBounds(), P.getUpperBounds());

    res.valid = true;
    return res;
}

#endif