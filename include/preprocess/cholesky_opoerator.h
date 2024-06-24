// VolEsti (volume computation and sampling library)

// Copyright (c) 2024 Vissarion Fisikopoulos
// Copyright (c) 2024 Apostolos Chalkis
// Copyright (c) 2024 Elias Tsigaridas

// Licensed under GNU LGPL.3, see LICENCE file


#ifndef CHOLESKY_OPERATOR_H
#define CHOLESKY_OPERATOR_H

#include <memory>

template<typename MT>
struct cholesky_operator
{
    inline static std::unique_ptr<Eigen::LLT<MT>>
    initialize(MT const&) 
    {
        return std::unique_ptr<Eigen::LLT<MT>>(new Eigen::LLT<MT>());
    }

    template <typename VT>
    inline static VT solve(std::unique_ptr<Eigen::LLT<MT>> const& llt, MT const& H, VT const& b)
    {
        llt->compute(H);
        return llt->solve(b);
    }

    template <typename diag_MT>
    inline static void update_hessian_Atrans_D_A(MT &H, MT const& A_trans, MT const& A, diag_MT const& D)
    {
        H.noalias() = A_trans * D * A;
    }

    inline static void init_Bmat(MT &B, int const n, MT const& , MT const& )
    {
        B.resize(n+1, n+1);
    }
    
    template <typename VT>
    inline static void update_Bmat(MT &B, VT const& AtDe, VT const& d, MT const& AtD, MT const& A)
    {
        const int n = A.cols();
        B.block(0, 0, n, n).noalias() = AtD * A;
        B.block(0, n, n, 1).noalias() = AtDe;
        B.block(n, 0, 1, n).noalias() = AtDe.transpose();
        B(n, n) = d.sum();
        B.noalias() += 1e-14 * MT::Identity(n + 1, n + 1);
    }
};

template <typename NT>
struct cholesky_operator<Eigen::SparseMatrix<NT>>
{
    inline static std::unique_ptr<Eigen::SimplicialLLT<Eigen::SparseMatrix<NT>>>
    initialize(Eigen::SparseMatrix<NT> const& mat) 
    {
        std::unique_ptr<Eigen::SimplicialLLT<Eigen::SparseMatrix<NT>>> llt = std::unique_ptr<Eigen::SimplicialLLT<Eigen::SparseMatrix<NT>>>(new Eigen::SimplicialLLT<Eigen::SparseMatrix<NT>>());
        llt->analyzePattern(mat);
        return llt;
    }

    template <typename VT>
    inline static VT solve(std::unique_ptr<Eigen::SimplicialLLT<Eigen::SparseMatrix<NT>>> const& llt, Eigen::SparseMatrix<NT> const& H, VT const& b)
    {
        llt->factorize(H);
        return llt->solve(b);
    }

    template <typename diag_MT>
    inline static void update_hessian_Atrans_D_A(Eigen::SparseMatrix<NT> &H, Eigen::SparseMatrix<NT> const& A_trans, Eigen::SparseMatrix<NT> const& A, diag_MT const& D)
    {
        H = A_trans * D * A;
    }

    inline static void init_Bmat(Eigen::SparseMatrix<NT> &B, int const n, Eigen::SparseMatrix<NT> const& A_trans, Eigen::SparseMatrix<NT> const& A)
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
                if (it.row() == it.col()) continue;
                trp.push_back(triplet(it.row(), it.col(), NT(1)));
            }
        }
        B.resize(n+1, n+1);
        B.setFromTriplets(trp.begin(), trp.end());
    }

    template <typename VT>
    inline static void update_Bmat(Eigen::SparseMatrix<NT> &B, VT const& AtDe, VT const& d, Eigen::SparseMatrix<NT> const& AtD, Eigen::SparseMatrix<NT> const& A)
    {
        const int n = A.cols();
        Eigen::SparseMatrix<NT> AtD_A = AtD * A;
        //std::cout<<"B indeces: "<<std::endl;
        int k = 0;
        //for (int k=0; k<AtD_A.outerSize(); ++k)
        while(k < B.outerSize())
        {
            //std::cout<<"k: "<<k<<std::endl;
            //typename Eigen::SparseMatrix<NT>::InnerIterator it1(B,k);
            typename Eigen::SparseMatrix<NT>::InnerIterator it2(AtD_A,k <= n-1 ? k : k-1);
            for (typename Eigen::SparseMatrix<NT>::InnerIterator it1(B,k); it1; ++it1)
            {
                //std::cout<<"row1: "<<it1.row()<<", col1: "<<it1.col()<<std::endl;
                //std::cout<<"row2: "<<it2.row()<<", col2: "<<it2.col()<<std::endl;
                
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
                else if (it1.row() == n && it1.col() == n)
                {
                    it1.valueRef() = d.sum();
                } 
                else
                {
                    std::cout<<"error in row,col"<<std::endl;
                    exit(-1);
                }

                if (it1.row() == it1.col())
                {
                    it1.valueRef() += 1e-14;
                    //it1.valueRef() += 1.0;
                }
                if (it1.row()<n-1) ++it2;
            }
            k++;
        }
        //std::cout<<"B:\n"<<Eigen::MatrixXd(B)<<std::endl;
        //std::cout<<"AtD_A:\n"<<Eigen::MatrixXd(AtD_A)<<std::endl;
        //std::cout<<"AtDe: "<<AtDe.transpose()<<std::endl;
        //std::cout<<"d.sum(): "<<d.sum()<<std::endl;
    }
};



#endif
