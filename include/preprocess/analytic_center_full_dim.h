// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2024 Vissarion Fisikopoulos
// Copyright (c) 2018-2024 Apostolos Chalkis
// Copyright (c) 2020-2024 Elias Tsigaridas

//Contributed by Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file


#ifndef ANALYTIC_CENTER_H
#define ANALYTIC_CENTER_H

#include <tuple>

#include "max_inscribed_ball.hpp"

template <typename VT, typename NT>
NT get_max_step(VT const& Ad, VT const& Ax, VT const& b)
{
    const int m = Ad.size();
    NT max_element = std::numeric_limits<NT>::lowest(), max_element_temp;
    for (int i = 0; i < m; i++)
    {
        max_element_temp = Ad.coeff(i) / (b.coeff(i) - Ax.coeff(i));
        if (max_element_temp > max_element) {
            max_element = max_element_temp;
        }
    }

    return NT(1) / max_element;
}

template <typename MT, typename VT, typename NT>
void get_hessian_grad_logbarrier(MT const& A, MT const& A_trans, VT const& b, 
                                 VT const& x, VT const& Ax, MT &H, VT &grad)
{
    const int m = A.rows();
    VT s(m);

    VT b_Ax = b - Ax;
    NT *s_data = s.data();
    for (int i = 0; i < m; i++)
    {
        *s_data = NT(1) / b_Ax.coeff(i);
        s_data++;
    }

    VT s_sq = s.cwiseProduct(s);
    grad.noalias() = A_trans * s;
    H.noalias() = A_trans * s_sq.asDiagonal() * A;    
}

template <typename MT, typename VT, typename NT>
std::tuple<VT, bool>  analytic_center_linear_ineq(MT const& A, VT const& b, 
                                                  unsigned int const maxiter = 500,
                                                  NT const grad_err_tol = 1e-08,
                                                  NT const rel_pos_err_tol = 1e-12) 
{
    VT x;
    bool feasibility_only = true, converged;
    std::tie(x, std::ignore, converged) =  max_inscribed_ball(A, b, maxiter, 1e-08, feasibility_only);
    VT Ax = A * x;
    if(!converged || (Ax.array() > b.array()).any()){
        std::runtime_error("The computation of the analytic center failed.");
    }
    const int n = A.cols(), m = A.rows();
    MT H(n, n);
    MT A_trans = A.transpose();
    VT grad(n);
    VT d(n);
    NT grad_err, rel_pos_err, rel_pos_err_temp, step;
    int iter = 0;
    VT Ad(m);
    VT step_d;
    VT x_prev;
    get_hessian_grad_logbarrier<MT, VT, NT>(A, A_trans, b, x, Ax, H, grad);
    converged = false;
    
    do {
        iter++;
        d.noalias() = - H.lu().solve(grad);
        Ad.noalias() = A * d;
        step = std::min(NT(0.99) * get_max_step<VT, NT>(Ad, Ax, b), NT(1));
        step_d.noalias() = step*d;
        x_prev = x;
        x.noalias() += step_d;
        Ax.noalias() += step*Ad;

        rel_pos_err = std::numeric_limits<NT>::lowest();
        for (int i = 0; i < n; i++) {
            rel_pos_err_temp = std::abs(step_d.coeff(i) / x_prev.coeff(i));
            if (rel_pos_err_temp > rel_pos_err) {
                rel_pos_err = rel_pos_err_temp;
            } 
        }
        
        get_hessian_grad_logbarrier<MT, VT, NT>(A, A_trans, b, x, Ax, H, grad);
        grad_err = grad.norm();

        if(iter > maxiter || grad_err <= grad_err_tol || rel_pos_err <= rel_pos_err_tol){
            converged = true;
            break;
        }
    } while(true);
    
    return std::make_tuple(x, converged);
}

#endif
