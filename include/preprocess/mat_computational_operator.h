// VolEsti (volume computation and sampling library)

// Copyright (c) 2024 Vissarion Fisikopoulos
// Copyright (c) 2024 Apostolos Chalkis
// Copyright (c) 2024 Elias Tsigaridas

// Licensed under GNU LGPL.3, see LICENCE file


#ifndef MAT_COMPUTATIONAL_OPERATOR_H
#define MAT_COMPUTATIONAL_OPERATOR_H

#include <memory>

#include "Spectra/include/Spectra/SymEigsSolver.h"
#include "Spectra/include/Spectra/MatOp/DenseSymMatProd.h"
#include "Spectra/include/Spectra/MatOp/SparseSymMatProd.h"


template<typename MT>
struct matrix_computational_operator {};


// Dense matrix operator
template<typename NT>
struct matrix_computational_operator<Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>>
{
    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;

    inline static std::unique_ptr<Eigen::LLT<MT>>
    initialize_chol(MT const&) 
    {
        return std::unique_ptr<Eigen::LLT<MT>>(new Eigen::LLT<MT>());
    }

    template <typename VT>
    inline static VT solve_vec(std::unique_ptr<Eigen::LLT<MT>> const& llt,
                               MT const& H, VT const& b)
    {
        llt->compute(H);
        return llt->solve(b);
    }

    inline static MT solve_mat(std::unique_ptr<Eigen::LLT<MT>> const& llt,
                               MT const& E, MT const& mat, NT &logdetE)
    {
        llt->compute(E);
        logdetE = llt->matrixL().toDenseMatrix().diagonal().array().log().sum();
        return llt->solve(mat);
    }

    template <typename diag_MT>
    inline static void update_Atrans_Diag_A(MT &H, MT const& A_trans,
                                            MT const& A, diag_MT const& D)
    {
        H.noalias() = A_trans * D * A;
    }

    template <typename diag_MT>
    inline static void update_Diag_A(MT &H, diag_MT const& D, MT const& A)
    {
        H.noalias() = D * A;
    }

    template <typename diag_MT>
    inline static void update_A_Diag(MT &H, MT const& A, diag_MT const& D)
    {
        H.noalias() = A * D;
    }

    inline static std::unique_ptr<Spectra::DenseSymMatProd<NT>> 
    get_mat_prod_op(MT const& E)
    {
        return std::make_unique<Spectra::DenseSymMatProd<NT>>(E);
    }
    
    inline static auto get_eigs_solver(std::unique_ptr<Spectra::DenseSymMatProd<NT>> const& op, int const n)
    {
        using SymDenseEigsSolver = Spectra::SymEigsSolver
          <
            NT, 
            Spectra::SELECT_EIGENVALUE::BOTH_ENDS, 
            Spectra::DenseSymMatProd<NT>
          >;
        // The value of ncv is chosen empirically
        return std::make_unique<SymDenseEigsSolver>(op.get(), 2, std::min(std::max(10, n/5), n));
    }

    inline static void init_Bmat(MT &B, int const n, MT const& , MT const& )
    {
        B.resize(n+1, n+1);
    }
    
    template <typename VT>
    inline static void update_Bmat(MT &B, VT const& AtDe, VT const& d,
                                   MT const& AtD, MT const& A)
    {
        const int n = A.cols();
        B.block(0, 0, n, n).noalias() = AtD * A;
        B.block(0, n, n, 1).noalias() = AtDe;
        B.block(n, 0, 1, n).noalias() = AtDe.transpose();
        B(n, n) = d.sum();
        B.noalias() += 1e-14 * MT::Identity(n + 1, n + 1);
    }
};


// Sparse matrix operator
template <typename NT>
struct matrix_computational_operator<Eigen::SparseMatrix<NT>>
{
    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;

    inline static std::unique_ptr<Eigen::SimplicialLLT<Eigen::SparseMatrix<NT>>>
    initialize_chol(Eigen::SparseMatrix<NT> const& mat) 
    {
        auto llt = std::make_unique<Eigen::SimplicialLLT<Eigen::SparseMatrix<NT>>>();
        llt->analyzePattern(mat);
        return llt;
    }

    template <typename VT>
    inline static VT solve_vec(std::unique_ptr<Eigen::SimplicialLLT<Eigen::SparseMatrix<NT>>> const& llt,
                               Eigen::SparseMatrix<NT> const& H, VT const& b)
    {
        llt->factorize(H);
        return llt->solve(b);
    }

    inline static MT solve_mat(std::unique_ptr<Eigen::SimplicialLLT<Eigen::SparseMatrix<NT>>> const& llt,
                               Eigen::SparseMatrix<NT> const& E, Eigen::SparseMatrix<NT> const& mat, NT &logdetE)
    {
        llt->factorize(E);
        logdetE = llt->matrixL().nestedExpression().diagonal().array().log().sum();
        return llt->solve(mat);
    }

    template <typename diag_MT>
    inline static void update_Atrans_Diag_A(Eigen::SparseMatrix<NT> &H,
                                            Eigen::SparseMatrix<NT> const& A_trans,
                                            Eigen::SparseMatrix<NT> const& A,
                                            diag_MT const& D)
    {
        H = A_trans * D * A;
    }

    template <typename diag_MT>
    inline static void update_Diag_A(Eigen::SparseMatrix<NT> &H,
                                     diag_MT const& D,
                                     Eigen::SparseMatrix<NT> const& A)
    {
        H = D * A;
    }

    template <typename diag_MT>
    inline static void update_A_Diag(Eigen::SparseMatrix<NT> &H,
                                     Eigen::SparseMatrix<NT> const& A,
                                     diag_MT const& D)
    {
        H = A * D;
    }

    inline static std::unique_ptr<Spectra::SparseSymMatProd<NT>>
    get_mat_prod_op(Eigen::SparseMatrix<NT> const& E)
    {
        return std::unique_ptr<Spectra::SparseSymMatProd<NT>>(new Spectra::SparseSymMatProd<NT>(E));
    }

    inline static auto get_eigs_solver(std::unique_ptr<Spectra::SparseSymMatProd<NT>> const& op, int const n)
    {
        using SymSparseEigsSolver = Spectra::SymEigsSolver
          <
            NT, 
            Spectra::SELECT_EIGENVALUE::BOTH_ENDS, 
            Spectra::SparseSymMatProd<NT>
          >;
        // The value of ncv is chosen empirically
        return std::make_unique<SymSparseEigsSolver>(op.get(), 2, std::min(std::max(10, n/5), n));
    }

    inline static void init_Bmat(Eigen::SparseMatrix<NT> &B, 
                                 int const n,
                                 Eigen::SparseMatrix<NT> const& A_trans,
                                 Eigen::SparseMatrix<NT> const& A)
    {
        // Initialize the structure of matrix B
        typedef Eigen::Triplet<NT> triplet;
        std::vector<triplet> trp;
        for (int i = 0; i < n; i++)
        {
            trp.push_back(triplet(i, i, NT(1)));
            trp.push_back(triplet(i, n, NT(1)));
            trp.push_back(triplet(n, i, NT(1)));
        }
        trp.push_back(triplet(n, n, NT(1)));
        
        Eigen::SparseMatrix<NT> ATA = A_trans * A;
        for (int k=0; k<ATA.outerSize(); ++k)
        {
            for (typename Eigen::SparseMatrix<NT>::InnerIterator it(ATA,k); it; ++it)
            {
                if (it.row() == it.col()) continue; // Diagonal element already allocated
                trp.push_back(triplet(it.row(), it.col(), NT(1)));
            }
        }
        B.resize(n+1, n+1);
        B.setFromTriplets(trp.begin(), trp.end());
    }

    template <typename VT>
    inline static void update_Bmat(Eigen::SparseMatrix<NT> &B,
                                   VT const& AtDe,
                                   VT const& d,
                                   Eigen::SparseMatrix<NT> const& AtD,
                                   Eigen::SparseMatrix<NT> const& A)
    {
        /* 
          B is (n+1)x(n+1) and AtD_A is nxn.
          We set B(1:n), 1:n) = AtD_A, B(n+1, :) = AtD^T, B(:, n+1) = AtD, B(n+1, n+1) = d.sum()
        */
        const int n = A.cols();
        Eigen::SparseMatrix<NT> AtD_A = AtD * A;
        int k = 0;
        while(k < B.outerSize())
        {
            typename Eigen::SparseMatrix<NT>::InnerIterator it2(AtD_A, k <= n-1 ? k : k-1);
            for (typename Eigen::SparseMatrix<NT>::InnerIterator it1(B, k); it1; ++it1)
            {                
                if (it1.row() <= n-1 && it1.col() <= n-1)
                {
                    it1.valueRef() = it2.value();
                }
                else if (it1.row() == n && it1.col() <= n-1)
                {
                    it1.valueRef() = AtDe.coeff(it1.col());
                }
                else if (it1.col() == n && it1.row() <= n-1)
                {
                    it1.valueRef() = AtDe.coeff(it1.row());
                }
                else // then, (it1.row() == n && it1.col() == n)
                {
                    it1.valueRef() = d.sum();
                }

                if (it1.row() == it1.col())
                {
                    it1.valueRef() += 1e-14;
                }
                if (it1.row()<n-1) ++it2;
            }
            k++;
        }
    }
};


#endif // MAT_COMPUTATIONAL_OPERATOR_H
