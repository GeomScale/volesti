// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2024 Vissarion Fisikopoulos
// Copyright (c) 2018-2024 Apostolos Chalkis

//Contributed and/or modified by Alexandros Manochis, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file


#ifndef ROUNDING_UTIL_FUNCTIONS_HPP
#define ROUNDING_UTIL_FUNCTIONS_HPP

#include <memory>

#include "Spectra/include/Spectra/SymEigsSolver.h"
#include "Spectra/include/Spectra/MatOp/DenseSymMatProd.h"
#include "Spectra/include/Spectra/MatOp/SparseSymMatProd.h"


enum EllipsoidType
{
  MAX_ELLIPSOID = 1,
  LOG_BARRIER = 2,
  VOLUMETRIC_BARRIER = 3,
  VAIDYA_BARRIER = 4
};

template <int T>
struct AssertBarrierFalseType : std::false_type {};

template <typename T>
struct AssertFalseType : std::false_type {};

template <typename VT>
auto get_max_step(VT const& Ad, VT const& b_Ax)
{
    using NT = typename VT::Scalar;
    const int m = Ad.size();
    NT max_element = std::numeric_limits<NT>::lowest(), max_element_temp;
    for (int i = 0; i < m; i++) {
        max_element_temp = Ad.coeff(i) / b_Ax.coeff(i);
        if (max_element_temp > max_element) {
            max_element = max_element_temp;
        }
    }

    return NT(1) / max_element;
}

template <typename NT, typename MT>
inline static auto
initialize_chol(MT const& mat) 
{
    using DenseMT = Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>;
    using SparseMT = Eigen::SparseMatrix<NT>;
    if constexpr (std::is_same<MT, DenseMT>::value)
    {
        return std::make_unique<Eigen::LLT<MT>>();
    } else if constexpr (std::is_base_of<Eigen::SparseMatrixBase<MT>, MT >::value)  
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
    } else if constexpr (std::is_base_of<Eigen::SparseMatrixBase<MT>, MT >::value)  
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
    } else if constexpr (std::is_base_of<Eigen::SparseMatrixBase<MT>, MT >::value)  
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
    } else if constexpr (std::is_base_of<Eigen::SparseMatrixBase<MT>, MT >::value)  
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
    } else if constexpr (std::is_base_of<Eigen::SparseMatrixBase<MT>, MT >::value)  
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
    } else if constexpr (std::is_base_of<Eigen::SparseMatrixBase<MT>, MT >::value)  
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
    } else if constexpr (std::is_base_of<Eigen::SparseMatrixBase<MT>, MT >::value)  
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
    } else if constexpr (std::is_base_of<Eigen::SparseMatrixBase<MT>, MT >::value)  
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
    } else if constexpr (std::is_base_of<Eigen::SparseMatrixBase<MT>, MT >::value)  
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
    } else if constexpr (std::is_base_of<Eigen::SparseMatrixBase<MT>, MT >::value)  
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

template <int BarrierType, typename NT>
std::tuple<NT, NT> init_step()
{
  if constexpr (BarrierType == EllipsoidType::LOG_BARRIER)
  {
    return {NT(1), NT(0.99)};
  } else if constexpr (BarrierType == EllipsoidType::VOLUMETRIC_BARRIER ||
                       BarrierType == EllipsoidType::VAIDYA_BARRIER)
  {
    return {NT(0.5), NT(0.4)};
  } else {
    static_assert(AssertBarrierFalseType<BarrierType>::value,
            "Barrier type is not supported.");
  }
}

template <typename MT_dense, int BarrierType, typename MT, typename VT, typename llt_type, typename NT>
void get_barrier_hessian_grad(MT const& A, MT const& A_trans, VT const& b, 
                              VT const& x, VT const& Ax, llt_type const& llt,
                              MT &H, VT &grad, VT &b_Ax, NT &obj_val)
{
  b_Ax.noalias() = b - Ax;
  VT s = b_Ax.cwiseInverse();
  VT s_sq = s.cwiseProduct(s);
  VT sigma;
  // Hessian of the log-barrier function
  update_Atrans_Diag_A<NT>(H, A_trans, A, s_sq.asDiagonal());

  if constexpr (BarrierType == EllipsoidType::VOLUMETRIC_BARRIER ||
                BarrierType == EllipsoidType::VAIDYA_BARRIER)
  {
    // Computing sigma(x)_i = (a_i^T H^{-1} a_i) / (b_i - a_i^Tx)^2
    MT_dense HA = solve_mat(llt, H, A_trans, obj_val);
    MT_dense aiHai = HA.transpose().cwiseProduct(A);
    sigma = (aiHai.rowwise().sum()).cwiseProduct(s_sq);
  }

  if constexpr (BarrierType == EllipsoidType::LOG_BARRIER)
  {
    grad.noalias() = A_trans * s;
  } else if constexpr (BarrierType == EllipsoidType::VOLUMETRIC_BARRIER)
  {
    // Gradient of the volumetric barrier function
    grad.noalias() = A_trans * (s.cwiseProduct(sigma));
    // Hessian of the volumetric barrier function
    update_Atrans_Diag_A<NT>(H, A_trans, A, s_sq.cwiseProduct(sigma).asDiagonal());
  } else if constexpr (BarrierType == EllipsoidType::VAIDYA_BARRIER)
  {
    const int m = b.size(), d = x.size();
    NT const d_m = NT(d) / NT(m);
    // Weighted gradient of the log barrier function
    grad.noalias() = A_trans * s;
    grad *= d_m;
    // Add the gradient of the volumetric function
    grad.noalias() += A_trans * (s.cwiseProduct(sigma));
    // Weighted Hessian of the log barrier function
    H *= d_m;
    // Add the Hessian of the volumetric function
    MT Hvol(d, d);
    update_Atrans_Diag_A<NT>(Hvol, A_trans, A, s_sq.cwiseProduct(sigma).asDiagonal());
    H += Hvol;
    obj_val -= s.array().log().sum();
  } else {
    static_assert(AssertBarrierFalseType<BarrierType>::value,
            "Barrier type is not supported.");
  }
}

template <int BarrierType, typename NT>
void get_step_next_iteration(NT const obj_val_prev, NT const obj_val,
                             NT const tol_obj, NT &step_iter)
{
  if constexpr (BarrierType == EllipsoidType::LOG_BARRIER)
  {
    step_iter = NT(1);
  } else if constexpr (BarrierType == EllipsoidType::VOLUMETRIC_BARRIER)
  {
    step_iter *= (obj_val_prev <= obj_val - tol_obj) ? NT(0.9) : NT(0.999);
  } else if constexpr (BarrierType == EllipsoidType::VAIDYA_BARRIER)
  {
    step_iter *= NT(0.999);
  } else {
    static_assert(AssertBarrierFalseType<BarrierType>::value,
            "Barrier type is not supported.");
  }
}

#endif // ROUNDING_UTIL_FUNCTIONS_HPP