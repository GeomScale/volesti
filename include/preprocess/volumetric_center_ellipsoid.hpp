// VolEsti (volume computation and sampling library)

// Copyright (c) 2024 Vissarion Fisikopoulos
// Copyright (c) 2024 Apostolos Chalkis
// Copyright (c) 2024 Elias Tsigaridas

// Licensed under GNU LGPL.3, see LICENCE file


#ifndef VOLUMETRIC_CENTER_ELLIPSOID_HPP
#define VOLUMETRIC_CENTER_ELLIPSOID_HPP

#include <tuple>

#include "preprocess/max_inscribed_ball.hpp"
#include "preprocess/feasible_point.hpp"
#include "preprocess/mat_computational_operators.hpp"

template <typename MT_dense, typename MT, typename VT, typename llt_type, typename NT>
void get_hessian_grad_volumetric_barrier(MT const& A, MT const& A_trans, VT const& b, 
                                         VT const& x, VT const& Ax, llt_type const& llt,
                                         MT &H, VT &grad, VT &b_Ax, NT &obj_val)
{
    b_Ax.noalias() = b - Ax;
    VT s = b_Ax.cwiseInverse();
    VT s_sq = s.cwiseProduct(s);
    // Hessian of the log-barrier function
    update_Atrans_Diag_A<NT>(H, A_trans, A, s_sq.asDiagonal());
    // Computing sigma(x)_i = (a_i^T H^{-1} a_i) / (b_i - a_i^Tx)^2
    MT_dense HA = solve_mat(llt, H, A_trans, obj_val);
    MT_dense aiHai = HA.transpose().cwiseProduct(A);
    VT sigma = (aiHai.rowwise().sum()).cwiseProduct(s_sq);
    // Gradient of the volumetric barrier function
    grad.noalias() = A_trans * (s.cwiseProduct(sigma));
    // Hessian of the volumetric barrier function
    update_Atrans_Diag_A<NT>(H, A_trans, A, s_sq.cwiseProduct(sigma).asDiagonal());
    //std::cout<<"H:\n"<<H<<std::endl;
}

/*
    This implementation computes the volumetric center of a polytope given 
    as a set of linear inequalities P = {x | Ax <= b}. The volumetric center
    is the minimizer of the volumetric barrier function i.e., the optimal solution
    of the following optimization problem,

    \min -, 
    
    The function solves the problem by using the Newton method.

    Input: (i)   Matrix A, vector b such that the polytope P = {x | Ax<=b}
           (ii)  The number of maximum iterations, max_iters
           (iii) Tolerance parameter grad_err_tol to bound the L2-norm of the gradient
           (iv)  Tolerance parameter rel_pos_err_tol to check the relative progress in each iteration

    Output: (i)   The Hessian computed on the volumetric center
            (ii)  The volumetric center of the polytope
            (iii) A boolean variable that declares convergence
    
    Note: Using MT as to deal with both dense and sparse matrices, MT_dense will be the type of result matrix
*/
template <typename MT_dense, typename MT, typename VT, typename NT>
std::tuple<MT_dense, VT, bool>  volumetric_center_ellipsoid_linear_ineq(MT const& A, VT const& b, VT const& x0,
                                                                        unsigned int const max_iters = 500,
                                                                        NT const grad_err_tol = 1e-08,
                                                                        NT const rel_pos_err_tol = 1e-12)
{
    // Initialization
    VT x = x0;
    VT Ax = A * x;
    const int n = A.cols(), m = A.rows();
    MT H(n, n), A_trans = A.transpose();
    VT grad(n), d(n), Ad(m), b_Ax(m), step_d(n), x_prev;
    NT grad_err, rel_pos_err, rel_pos_err_temp, step, obj_val, obj_val_prev;
    unsigned int iter = 0;
    bool converged = false;
    const NT tol_bnd = NT(0.01), tol_obj = NT(1e-06);
    NT step_iter = NT(0.5);

    auto llt = initialize_chol<NT>(A_trans, A);
    get_hessian_grad_volumetric_barrier<MT_dense>(A, A_trans, b, x, Ax, llt,
                                                  H, grad, b_Ax, obj_val);
    do {
        iter++;
        // Compute the direction
        d.noalias() = - solve_vec<NT>(llt, H, grad);
        Ad.noalias() = A * d;
        // Compute the step length
        step = std::min(NT(0.4) * get_max_step<VT, NT>(Ad, b_Ax), step_iter);
        step_d.noalias() = step*d;
        x_prev = x;
        x += step_d;
        Ax.noalias() += step*Ad;

        // Compute the max_i\{ |step*d_i| ./ |x_i| \} 
        rel_pos_err = std::numeric_limits<NT>::lowest();
        for (int i = 0; i < n; i++)
        {
            rel_pos_err_temp = std::abs(step_d.coeff(i) / x_prev.coeff(i));
            if (rel_pos_err_temp > rel_pos_err)
            {
                rel_pos_err = rel_pos_err_temp;
            }
        }
        
        obj_val_prev = obj_val;
        get_hessian_grad_volumetric_barrier<MT_dense>(A, A_trans, b, x, Ax, llt,
                                                      H, grad, b_Ax, obj_val);
        grad_err = grad.norm();
        std::cout<<"iter: "<<iter<<", grad_err: "<<grad_err<<", obj_val: "<<obj_val<<"\n------------\n"<<std::endl;

        if (iter >= max_iters || grad_err <= grad_err_tol || rel_pos_err <= rel_pos_err_tol)
        {
            converged = true;
            break;
        }
        step_iter *= NT(0.999);
    } while (true);
    
    return std::make_tuple(MT_dense(H), x, converged);
}

template <typename MT_dense, typename MT, typename VT, typename NT>
std::tuple<MT_dense, VT, bool>  volumetric_center_ellipsoid_linear_ineq(MT const& A, VT const& b,
                                                                        unsigned int const max_iters = 500,
                                                                        NT const grad_err_tol = 1e-08,
                                                                        NT const rel_pos_err_tol = 1e-12) 
{
    VT x0 = compute_feasible_point(A, b);
    std::cout<<"x0: "<<x0.transpose()<<std::endl;
    return volumetric_center_ellipsoid_linear_ineq<MT_dense, MT, VT, NT>(A, b, x0, max_iters, grad_err_tol, rel_pos_err_tol);
}

#endif // VOLUMETRIC_CENTER_ELLIPSOID_HPP
