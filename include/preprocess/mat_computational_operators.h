// VolEsti (volume computation and sampling library)

// Copyright (c) 2024 Vissarion Fisikopoulos
// Copyright (c) 2024 Apostolos Chalkis
// Copyright (c) 2024 Elias Tsigaridas

// Licensed under GNU LGPL.3, see LICENCE file


#ifndef MAT_COMPUTATIONAL_OPERATORS_H
#define MAT_COMPUTATIONAL_OPERATORS_H

#include <memory>

#include "Spectra/include/Spectra/SymEigsSolver.h"
#include "Spectra/include/Spectra/MatOp/DenseSymMatProd.h"
#include "Spectra/include/Spectra/MatOp/SparseSymMatProd.h"


template <typename T>
struct AssertFalseType : std::false_type {};

template <typename NT, typename MT>
inline static auto
initialize_chol(MT const& mat) 
{
    using DenseMT = Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>;
    using SparseMT = Eigen::SparseMatrix<NT>;
    if constexpr (std::is_same<MT, DenseMT>::value)
    {
        return std::make_unique<Eigen::LLT<MT>>();
    } else if constexpr (std::is_same<MT, SparseMT>::value)  
    {
        auto llt = std::make_unique<Eigen::SimplicialLLT<MT>>();
        llt->analyzePattern(mat);
        return llt;
    } else 
    {
        static_assert(AssertFalseType<MT>::value,
            "Matrix type is not supported.");
    }
}

template <typename NT, typename MT>
inline static auto
initialize_chol(MT const& A_trans, MT const& A) 
{
    using DenseMT = Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>;
    using SparseMT = Eigen::SparseMatrix<NT>;
    if constexpr (std::is_same<MT, DenseMT>::value)
    {
        return std::make_unique<Eigen::LLT<MT>>();
    } else if constexpr (std::is_same<MT, SparseMT>::value)  
    {
        MT mat = A_trans * A;
        return initialize_chol<NT>(mat);
    } else 
    {
        static_assert(AssertFalseType<MT>::value,
            "Matrix type is not supported.");
    }
}

template <typename NT, typename MT, typename Eigen_lltMT,  typename VT>
inline static VT solve_vec(std::unique_ptr<Eigen_lltMT> const& llt,
                           MT const& H, VT const& b)
{
    using DenseMT = Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>;
    using SparseMT = Eigen::SparseMatrix<NT>;
    if constexpr (std::is_same<MT, DenseMT>::value)
    {
        llt->compute(H);
        return llt->solve(b);
    } else if constexpr (std::is_same<MT, SparseMT>::value)  
    {
        llt->factorize(H);
        return llt->solve(b);
    } else
    {
        static_assert(AssertFalseType<MT>::value,
            "Matrix type is not supported.");
    }
}

template <typename Eigen_lltMT, typename MT, typename NT>
inline static Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>
solve_mat(std::unique_ptr<Eigen_lltMT> const& llt,
          MT const& H, MT const& mat, NT &logdetE)
{
    using DenseMT = Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>;
    using SparseMT = Eigen::SparseMatrix<NT>;
    if constexpr (std::is_same<MT, DenseMT>::value)
    {
        llt->compute(H);
        logdetE = llt->matrixL().toDenseMatrix().diagonal().array().log().sum();
        return llt->solve(mat);
    } else if constexpr (std::is_same<MT, SparseMT>::value)  
    {
        llt->factorize(H);
        logdetE = llt->matrixL().nestedExpression().diagonal().array().log().sum();
        return llt->solve(mat);
    } else 
    {
        static_assert(AssertFalseType<MT>::value,
            "Matrix type is not supported.");
    }
}

template <typename NT, typename MT, typename diag_MT>
inline static void update_Atrans_Diag_A(MT &H, MT const& A_trans,
                                        MT const& A, diag_MT const& D)
{
    using DenseMT = Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>;
    using SparseMT = Eigen::SparseMatrix<NT>;
    if constexpr (std::is_same<MT, DenseMT>::value)
    {
        H.noalias() = A_trans * D * A;
    } else if constexpr (std::is_same<MT, SparseMT>::value)  
    {
        H = A_trans * D * A;
    } else 
    {
        static_assert(AssertFalseType<MT>::value,
            "Matrix type is not supported.");
    }
}

template <typename NT, typename MT, typename diag_MT>
inline static void update_Diag_A(MT &H, diag_MT const& D, MT const& A)
{
    using DenseMT = Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>;
    using SparseMT = Eigen::SparseMatrix<NT>;
    if constexpr (std::is_same<MT, DenseMT>::value)
    {
        H.noalias() = D * A;
    } else if constexpr (std::is_same<MT, SparseMT>::value)  
    {
        H = D * A;
    } else 
    {
        static_assert(AssertFalseType<MT>::value,
            "Matrix type is not supported.");
    }
}

template <typename NT, typename MT, typename diag_MT>
inline static void update_A_Diag(MT &H, MT const& A, diag_MT const& D)
{
    using DenseMT = Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>;
    using SparseMT = Eigen::SparseMatrix<NT>;
    if constexpr (std::is_same<MT, DenseMT>::value)
    {
        H.noalias() = A * D;
    } else if constexpr (std::is_same<MT, SparseMT>::value)  
    {
        H = A * D;
    } else 
    {
        static_assert(AssertFalseType<MT>::value,
            "Matrix type is not supported.");
    }
}

template <typename NT, typename MT>
inline static auto
get_mat_prod_op(MT const& E)
{
    using DenseMT = Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>;
    using SparseMT = Eigen::SparseMatrix<NT>;
    if constexpr (std::is_same<MT, DenseMT>::value)
    {
        return std::make_unique<Spectra::DenseSymMatProd<NT>>(E);
    } else if constexpr (std::is_same<MT, SparseMT>::value)  
    {
        return std::make_unique<Spectra::SparseSymMatProd<NT>>(E);
    } else 
    {
        static_assert(AssertFalseType<MT>::value,
            "Matrix type is not supported.");
    }
}

template<typename NT, typename SpectraMatProdNT>
inline static auto get_eigs_solver(std::unique_ptr<SpectraMatProdNT> const& op, int const n)
{
    using DenseMatProd = Spectra::DenseSymMatProd<NT>;
    using SparseMatProd = Spectra::SparseSymMatProd<NT>;
    if constexpr (std::is_same<SpectraMatProdNT, DenseMatProd>::value)
    {
        using SymDenseEigsSolver = Spectra::SymEigsSolver
          <
            NT, 
            Spectra::SELECT_EIGENVALUE::BOTH_ENDS, 
            DenseMatProd
          >;
        // The value of ncv is chosen empirically
        return std::make_unique<SymDenseEigsSolver>(op.get(), 2, std::min(std::max(10, n/5), n));
    } else if constexpr (std::is_same<SpectraMatProdNT, SparseMatProd>::value)  
    {
        using SymSparseEigsSolver = Spectra::SymEigsSolver
          <
            NT, 
            Spectra::SELECT_EIGENVALUE::BOTH_ENDS, 
            SparseMatProd
          >;
        // The value of ncv is chosen empirically
        return std::make_unique<SymSparseEigsSolver>(op.get(), 2, std::min(std::max(10, n/5), n));
    } else 
    {
        static_assert(AssertFalseType<SpectraMatProdNT>::value,
            "Matrix-vector multiplication multiplication is not supported.");
    }
}

template <typename NT, typename MT>
inline static void
init_Bmat(MT &B, int const n, MT const& A_trans, MT const& A) 
{
    using DenseMT = Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>;
    using SparseMT = Eigen::SparseMatrix<NT>;
    if constexpr (std::is_same<MT, DenseMT>::value)
    {
        B.resize(n+1, n+1);
    } else if constexpr (std::is_same<MT, SparseMT>::value)  
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
        
        MT ATA = A_trans * A;
        for (int k=0; k<ATA.outerSize(); ++k)
        {
            for (typename MT::InnerIterator it(ATA,k); it; ++it)
            {
                if (it.row() == it.col()) continue; // Diagonal elements are already allocated
                trp.push_back(triplet(it.row(), it.col(), NT(1)));
            }
        }
        B.resize(n+1, n+1);
        B.setFromTriplets(trp.begin(), trp.end());
    } else 
    {
        static_assert(AssertFalseType<MT>::value,
            "Matrix type is not supported.");
    }
}

template <typename NT, typename MT, typename VT>
inline static void
update_Bmat(MT &B, VT const& AtDe, VT const& d,
            MT const& AtD, MT const& A) 
{
    using DenseMT = Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>;
    using SparseMT = Eigen::SparseMatrix<NT>;
    const int n = A.cols();
    /* 
        B is (n+1)x(n+1) and AtD_A is nxn.
        We set B(1:n), 1:n) = AtD_A, B(n+1, :) = AtD^T, B(:, n+1) = AtD, B(n+1, n+1) = d.sum()
    */
    if constexpr (std::is_same<MT, DenseMT>::value)
    {
        B.block(0, 0, n, n).noalias() = AtD * A;
        B.block(0, n, n, 1).noalias() = AtDe;
        B.block(n, 0, 1, n).noalias() = AtDe.transpose();
        B(n, n) = d.sum();
        B.noalias() += 1e-14 * MT::Identity(n + 1, n + 1);
    } else if constexpr (std::is_same<MT, SparseMT>::value)  
    {
        MT AtD_A = AtD * A;
        int k = 0;
        while(k < B.outerSize())
        {
            typename MT::InnerIterator it2(AtD_A, k <= n-1 ? k : k-1);
            for (typename MT::InnerIterator it1(B, k); it1; ++it1)
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
    } else 
    {
        static_assert(AssertFalseType<MT>::value,
            "Matrix type is not supported.");
    }
}


#endif // MAT_COMPUTATIONAL_OPERATORS_H
