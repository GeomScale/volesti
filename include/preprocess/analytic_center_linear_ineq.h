// VolEsti (volume computation and sampling library)

// Copyright (c) 2024 Vissarion Fisikopoulos
// Copyright (c) 2024 Apostolos Chalkis
// Copyright (c) 2024 Elias Tsigaridas

// Licensed under GNU LGPL.3, see LICENCE file


#ifndef ANALYTIC_CENTER_H
#define ANALYTIC_CENTER_H

#include <tuple>

#include "preprocess/max_inscribed_ball.hpp"
#include "preprocess/feasible_point.hpp"
#include "preprocess/mat_computational_operators.h"

template <typename VT, typename NT>
NT get_max_step(VT const& Ad, VT const& b_Ax)
{
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

template <typename MT, typename VT, typename NT>
void get_hessian_grad_logbarrier(MT const& A, MT const& A_trans, VT const& b, 
                                 VT const& x, VT const& Ax, MT &H, VT &grad, VT &b_Ax)
{
    b_Ax.noalias() = b - Ax;
    VT s = b_Ax.cwiseInverse();
    VT s_sq = s.cwiseProduct(s);
    // Gradient of the log-barrier function
    grad.noalias() = A_trans * s;
    // Hessian of the log-barrier function
    update_Atrans_Diag_A<NT>(H, A_trans, A, s_sq.asDiagonal());
}

/*
    This implementation computes the analytic center of a polytope given 
    as a set of linear inequalities P = {x | Ax <= b}. The analytic center
    is the minimizer of the log barrier function i.e., the optimal solution
    of the following optimization problem (Convex Optimization, Boyd and Vandenberghe, Section 8.5.3),

    \min -\sum \log(b_i - a_i^Tx), where a_i is the i-th row of A.
    
    The function solves the problem by using the Newton method.

    Input: (i)   Matrix A, vector b such that the polytope P = {x | Ax<=b}
           (ii)  The number of maximum iterations, max_iters
           (iii) Tolerance parameter grad_err_tol to bound the L2-norm of the gradient
           (iv)  Tolerance parameter rel_pos_err_tol to check the relative progress in each iteration

    Output: (i)   The Hessian computed on the analytic center
            (ii)  The analytic center of the polytope
            (iii) A boolean variable that declares convergence
    
    Note: Using MT as to deal with both dense and sparse matrices, MT_dense will be the type of result matrix
*/
template <typename MT_dense, typename MT, typename VT, typename NT>
std::tuple<MT_dense, VT, bool>  analytic_center_linear_ineq(MT const& A, VT const& b, VT const& x0,
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
    NT grad_err, rel_pos_err, rel_pos_err_temp, step;
    unsigned int iter = 0;
    bool converged = false;
    const NT tol_bnd = NT(0.01);

    auto llt = initialize_chol<NT>(A_trans, A);
    get_hessian_grad_logbarrier<MT, VT, NT>(A, A_trans, b, x, Ax, H, grad, b_Ax);
    
    do {
        iter++;
        // Compute the direction
        d.noalias() = - solve_vec<NT>(llt, H, grad);
        Ad.noalias() = A * d;
        // Compute the step length
        step = std::min((NT(1) - tol_bnd) * get_max_step<VT, NT>(Ad, b_Ax), NT(1));
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
        
        get_hessian_grad_logbarrier<MT, VT, NT>(A, A_trans, b, x, Ax, H, grad, b_Ax);
        grad_err = grad.norm();

        if (iter >= max_iters || grad_err <= grad_err_tol || rel_pos_err <= rel_pos_err_tol)
        {
            converged = true;
            break;
        }
    } while (true);
    
    return std::make_tuple(MT_dense(H), x, converged);
}

template <typename MT_dense, typename MT, typename VT, typename NT>
std::tuple<MT_dense, VT, bool>  analytic_center_linear_ineq(MT const& A, VT const& b,
                                                            unsigned int const max_iters = 500,
                                                            NT const grad_err_tol = 1e-08,
                                                            NT const rel_pos_err_tol = 1e-12) 
{
    VT x0 = compute_feasible_point(A, b);
    return analytic_center_linear_ineq<MT_dense, MT, VT, NT>(A, b, x0, max_iters, grad_err_tol, rel_pos_err_tol);
}

#endif
