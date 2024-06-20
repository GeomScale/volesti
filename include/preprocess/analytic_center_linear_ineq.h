// VolEsti (volume computation and sampling library)

// Copyright (c) 2024 Vissarion Fisikopoulos
// Copyright (c) 2024 Apostolos Chalkis
// Copyright (c) 2024 Elias Tsigaridas

// Licensed under GNU LGPL.3, see LICENCE file


#ifndef ANALYTIC_CENTER_H
#define ANALYTIC_CENTER_H

#include <tuple>

#include "max_inscribed_ball.hpp"

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
    const int m = A.rows();
    VT s(m);

    b_Ax.noalias() = b - Ax;
    NT *s_data = s.data();
    for (int i = 0; i < m; i++) {
        *s_data = NT(1) / b_Ax.coeff(i);
        s_data++;
    }

    VT s_sq = s.cwiseProduct(s);
    // Gradient of the log-barrier function
    grad.noalias() = A_trans * s;
    // Hessian of the log-barrier function
    H.noalias() = A_trans * s_sq.asDiagonal() * A;    
}

/*
    This implementation computes the analytic center of a polytope given 
    as a set of linear inequalities P = {x | Ax <= b}. The analytic center
    is the optimal solution of the following optimization problem,
    \min - \sum \log(b_i - a_i^Tx), where a_i is the i-th row of A.
    The function implements the Newton method.

    Input: (i)   Matrix A, vector b such that the polytope P = {x | Ax<=b}
           (ii)  The number of maximum iterations, max_iters
           (iii) Tolerance parameter grad_err_tol to bound the L2-norm of the gradient
           (iv)  Tolerance parameter rel_pos_err_tol to check the relative progress in each iteration

    Output: (i)  The analytic center of the polytope
            (ii) A boolean variable that declares convergence
*/
template <typename MT, typename VT, typename NT>
std::tuple<VT, bool>  analytic_center_linear_ineq(MT const& A, VT const& b, 
                                                  unsigned int const max_iters = 500,
                                                  NT const grad_err_tol = 1e-08,
                                                  NT const rel_pos_err_tol = 1e-12) 
{
    VT x;
    bool feasibility_only = true, converged;
    // Compute a feasible point
    std::tie(x, std::ignore, converged) =  max_inscribed_ball(A, b, max_iters, 1e-08, feasibility_only);
    VT Ax = A * x;
    if(!converged || (Ax.array() > b.array()).any()) {
        std::runtime_error("The computation of the analytic center failed.");
    }
    // Initialization
    const int n = A.cols(), m = A.rows();
    MT H(n, n), A_trans = A.transpose();
    VT grad(n), d(n), Ad(m), b_Ax(m), step_d(n), x_prev;
    NT grad_err, rel_pos_err, rel_pos_err_temp, step;
    unsigned int iter = 0;
    converged = false;

    get_hessian_grad_logbarrier<MT, VT, NT>(A, A_trans, b, x, Ax, H, grad, b_Ax);
    
    do {
        iter++;
        // Compute the direction
        d.noalias() = - H.lu().solve(grad);
        Ad.noalias() = A * d;
        // Compute the step length
        step = std::min(NT(0.99) * get_max_step<VT, NT>(Ad, b_Ax), NT(1));
        step_d.noalias() = step*d;
        x_prev = x;
        x += step_d;
        Ax.noalias() += step*Ad;

        // Compute the max_i\{ |step*d_i| ./ |x_i| \} 
        rel_pos_err = std::numeric_limits<NT>::lowest();
        for (int i = 0; i < n; i++) {
            rel_pos_err_temp = std::abs(step_d.coeff(i) / x_prev.coeff(i));
            if (rel_pos_err_temp > rel_pos_err) {
                rel_pos_err = rel_pos_err_temp;
            } 
        }
        
        get_hessian_grad_logbarrier<MT, VT, NT>(A, A_trans, b, x, Ax, H, grad, b_Ax);
        grad_err = grad.norm();

        if(iter >= max_iters || grad_err <= grad_err_tol || rel_pos_err <= rel_pos_err_tol) {
            converged = true;
            break;
        }
    } while(true);
    
    return std::make_tuple(x, converged);
}

#endif
