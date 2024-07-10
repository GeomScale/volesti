// VolEsti (volume computation and sampling library)

// Copyright (c) 2024 Vissarion Fisikopoulos
// Copyright (c) 2024 Apostolos Chalkis
// Copyright (c) 2024 Elias Tsigaridas

// Licensed under GNU LGPL.3, see LICENCE file


#ifndef BARRIER_CENTER_ELLIPSOID_HPP
#define BARRIER_CENTER_ELLIPSOID_HPP

#include <tuple>

#include "preprocess/max_inscribed_ball.hpp"
#include "preprocess/feasible_point.hpp"
#include "preprocess/rounding_util_functions.hpp"

/*
    This implementation computes the analytic or the volumetric center of a polytope given 
    as a set of linear inequalities P = {x | Ax <= b}. The analytic center is the tminimizer
    of the log barrier function, i.e., the optimal solution
    of the following optimization problem (Convex Optimization, Boyd and Vandenberghe, Section 8.5.3),

    \min -\sum \log(b_i - a_i^Tx), where a_i is the i-th row of A.
    
    The volumetric center is the minimizer of the volumetric barrier function, i.e., the optimal
    solution of the following optimization problem,

    \min logdet \nabla^2 f(x), where f(x) the log barrier function

    The Vaidya center is the minimizer of the Vaidya barrier function, i.e., the optimal
    solution of the following optimization problem,

    \min logdet \nabla^2 f(x) + d/m f(x), where f(x) the log barrier function.
    
    The function solves the problems by using the Newton method.

    Input: (i)   Matrix A, vector b such that the polytope P = {x | Ax<=b}
           (ii)  The number of maximum iterations, max_iters
           (iii) Tolerance parameter grad_err_tol to bound the L2-norm of the gradient
           (iv)  Tolerance parameter rel_pos_err_tol to check the relative progress in each iteration

    Output: (i)   The Hessian of the barrier function
            (ii)  The analytic/volumetric center of the polytope
            (iii) A boolean variable that declares convergence
    
    Note: Using MT as to deal with both dense and sparse matrices, MT_dense will be the type of result matrix
*/
template <typename MT_dense, int BarrierType, typename NT, typename MT, typename VT>
std::tuple<MT_dense, VT, bool>  barrier_center_ellipsoid_linear_ineq(MT const& A, VT const& b, VT const& x0,
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

    auto [step_iter, max_step_multiplier] = init_step<BarrierType, NT>();
    auto llt = initialize_chol<NT>(A_trans, A);
    get_barrier_hessian_grad<MT_dense, BarrierType>(A, A_trans, b, x, Ax, llt,
                                                    H, grad, b_Ax, obj_val);
    do {
        iter++;
        // Compute the direction
        d.noalias() = - solve_vec<NT>(llt, H, grad);
        Ad.noalias() = A * d;
        // Compute the step length
        step = std::min(max_step_multiplier * get_max_step(Ad, b_Ax), step_iter);
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
        get_barrier_hessian_grad<MT_dense, BarrierType>(A, A_trans, b, x, Ax, llt,
                                                        H, grad, b_Ax, obj_val);
        grad_err = grad.norm();

        if (iter >= max_iters || grad_err <= grad_err_tol || rel_pos_err <= rel_pos_err_tol)
        {
            converged = true;
            break;
        }
        get_step_next_iteration<BarrierType>(obj_val_prev, obj_val, tol_obj, step_iter);
    } while (true);
    
    return std::make_tuple(MT_dense(H), x, converged);
}

template <typename MT_dense, int BarrierType, typename NT, typename MT, typename VT>
std::tuple<MT_dense, VT, bool>  barrier_center_ellipsoid_linear_ineq(MT const& A, VT const& b,
                                                                     unsigned int const max_iters = 500,
                                                                     NT const grad_err_tol = 1e-08,
                                                                     NT const rel_pos_err_tol = 1e-12) 
{
    VT x0 = compute_feasible_point(A, b);
    return barrier_center_ellipsoid_linear_ineq<MT_dense, BarrierType>(A, b, x0, max_iters, grad_err_tol, rel_pos_err_tol);
}

#endif // BARRIER_CENTER_ELLIPSOID_HPP
