// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

//Contributed and/or modified by Alexandros Manochis, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file


#ifndef MAX_INNER_BALL
#define MAX_INNER_BALL

/*
    This implmentation computes the largest inscribed ball in a given convex polytope P.
    The polytope has to be given in H-representation P = {x | Ax <= b} and the rows of A
    has to be normalized. It solves the Linear program: max t, s.t. Ax + t*e <= b, where 
    e is the vector of ones.

    The implementation is based on Yin Zhang's Matlab implementation in https://github.com/Bounciness/Volume-and-Sampling/blob/1c7adfb46c2c01037e625db76ff00e73616441d4/external/mve11/mve_cobra/mve_presolve_cobra.m

    Input: matrix A, vector b such that the polytope P = {x | Ax<=b}
           tolerance parameter tol

    Output: center of the ball x
            radius r
*/

template <typename MT, typename VT, typename NT>
void calcstep(MT const& A, MT const& A_trans, MT const& R, VT &s, 
              VT &y, VT &r1, VT const& r2, NT const& r3, VT &r4,
              VT &dx, VT &ds, NT &dt, VT &dy, VT &tmp, VT &rhs)
{
    int m = A.rows(), n = A.cols();
    NT *vec_iter1 = tmp.data(), *vec_iter2 = y.data(), *vec_iter3 = s.data(),
       *vec_iter4 = r1.data(), *vec_iter5 = r4.data();
    for (int i = 0; i < m; ++i) {
        *vec_iter1 = ((*vec_iter4) / (*vec_iter2) - (*vec_iter5)) / (*vec_iter3);
        vec_iter1++; vec_iter2++; vec_iter3++; vec_iter4++; vec_iter5++;
    }

    rhs.block(0,0,n,1).noalias() = r2 + A_trans * tmp;
    rhs(n) = r3 + tmp.sum();
    VT dxdt = R.colPivHouseholderQr().solve(R.transpose().colPivHouseholderQr().solve(rhs));

    dx = dxdt.block(0,0,n,1);
    dt = dxdt(n);
    ds.noalias() = r1 - A*dx - VT::Ones(m) * dt;
    vec_iter1 = dy.data(); vec_iter2 = r4.data(); vec_iter3 = y.data();
    vec_iter4 = ds.data(); vec_iter5 = s.data();

    for (int i = 0; i < m; ++i) {
        *vec_iter1 = ((*vec_iter2) - (*vec_iter3) * (*vec_iter4)) / (*vec_iter5);
        vec_iter1++; vec_iter2++; vec_iter3++; vec_iter4++; vec_iter5++;
    }
}


template <typename MT, typename VT, typename NT>
std::tuple<VT, NT, bool>  max_inscribed_ball(MT const& A, VT const& b, unsigned int maxiter, NT tol) 
{
    int m = A.rows(), n = A.cols();
    bool converge;

    NT bnrm = b.norm();
    VT o_m = VT::Zero(m), o_n = VT::Zero(n), e_m = VT::Ones(m);

    VT x = o_n, y = e_m / m;
    NT t = b.minCoeff() - 1.0;
    VT s = b - e_m * t;

    VT dx = o_n;
    VT dxc = dx, ds = o_m;
    VT dsc = ds, dy = o_m, mu_ds_dy(m), tmp(m), rhs(n + 1);
    VT dyc = dy, r1(m), r2(n), r4(m), r23(n + 1), AtDe(n), d(m);

    NT dt = NT(0), dtc = NT(0), tau, sigma0 = 0.2, r3, gap, 
       prif, drif, rgap, total_err, alphap, alphad,
       ratio, sigma, mu, t_prev = 1000.0 * t + 100.0;
    
    NT const tau0 = 0.995, power_num = 5.0 * std::pow(10.0, 15.0);
    NT *vec_iter1, *vec_iter2, *vec_iter3, *vec_iter4;

    MT B(n + 1, n + 1), AtD(n, m), R(n + 1, n + 1), 
       eEye = std::pow(10.0, -14.0) * MT::Identity(n + 1, n + 1),
       A_trans = A.transpose();

    for (unsigned int i = 0; i < maxiter; ++i) {

        // KKT residuals
        r1.noalias() = b - (A * x + s + t * e_m);
        r2.noalias() = -A_trans * y;
        r3 = 1.0 - y.sum();
        r4 = -s.cwiseProduct(y);

        r23.block(0, 0, n, 1) = r2;
        r23(n) = r3;
        gap = -r4.sum();

        // relative residual norms and gap
        prif = r1.norm() / (1.0 + bnrm);
        drif = r23.norm() / 10.0;
        rgap = std::abs(b.dot(y) - t) / (1.0 + std::abs(t));
        total_err = std::max(prif, drif);
        total_err = std::max(total_err, rgap);

        // progress output & check stopping
        if (total_err < tol || ( t > 0 && ( (std::abs(t - t_prev) <= tol * t && std::abs(t - t_prev) <= tol * t_prev) 
                    || (t_prev >= (1.0 - tol) * t && i > 0) 
                    || (t <= (1.0 - tol) * t_prev && i > 0) ) ) ) 
        {
            //converged
            converge = true;
            break;
        }

        if (dt > 1000.0 * bnrm || t > 1000000.0 * bnrm) 
        {
            //unbounded
            converge = false;
            break;
        }

        // Shur complement matrix
        vec_iter1 = d.data();
        vec_iter3 = s.data();
        vec_iter2 = y.data();
        for (int j = 0; j < m; ++j) {
            *vec_iter1 = std::min(power_num, (*vec_iter2) / (*vec_iter3));
            AtD.col(j).noalias() = A_trans.col(j) * (*vec_iter1);
            vec_iter1++;
            vec_iter3++;
            vec_iter2++;
        }

        AtDe.noalias() = AtD * e_m;
        B.block(0, 0, n, n).noalias() = AtD * A;
        B.block(0, n, n, 1).noalias() = AtDe;
        B.block(n, 0, 1, n).noalias() = AtDe.transpose();
        B(n, n) = d.sum();
        B.noalias() += eEye;

        // Cholesky decomposition
        Eigen::LLT <MT> lltOfB(B);
        R = lltOfB.matrixL().transpose();

        // predictor step & length
        calcstep(A, A_trans, R, s, y, r1, r2, r3, r4, dx, ds, dt, dy, tmp, rhs);

        alphap = -1.0;
        alphad = -1.0;
        vec_iter1 = ds.data();
        vec_iter2 = s.data();
        vec_iter3 = dy.data();
        vec_iter4 = y.data();

        for (int j = 0; j < m; ++j) {
            alphap = std::min(alphap, (*vec_iter1) / (*vec_iter2));
            alphad = std::min(alphad, (*vec_iter3) / (*vec_iter4));
            vec_iter1++;
            vec_iter2++;
            vec_iter3++;
            vec_iter4++;
        }
        alphap = -1.0 / alphap;
        alphad = -1.0 / alphad;

        // determine mu
        ratio = (s + alphap * ds).dot((y + alphad * dy)) / gap;
        sigma = std::min(sigma0, ratio * ratio);
        mu = (sigma * gap) / NT(m);

        // corrector and combined step & length
        mu_ds_dy.noalias() = e_m * mu - ds.cwiseProduct(dy);
        calcstep(A, A_trans, R, s, y, o_m, o_n, 0.0, mu_ds_dy, dxc, dsc, dtc, dyc, tmp, rhs);

        dx += dxc;
        ds += dsc;
        dt += dtc;
        dy += dyc;

        alphap = -0.5;
        alphad = -0.5;
        vec_iter1 = ds.data();
        vec_iter2 = s.data();
        vec_iter3 = dy.data();
        vec_iter4 = y.data();

        for (int j = 0; j < m; ++j) {
            alphap = std::min(alphap, (*vec_iter1) / (*vec_iter2));
            alphad = std::min(alphad, (*vec_iter3) / (*vec_iter4));
            vec_iter1++;
            vec_iter2++;
            vec_iter3++;
            vec_iter4++;
        }
        alphap = -1.0 / alphap;
        alphad = -1.0 / alphad;

        // update iterates
        tau = std::max(tau0, 1.0 - gap / NT(m));
        alphap = std::min(1.0, tau * alphap);
        alphad = std::min(1.0, tau * alphad);

        x += alphap * dx;
        s += alphap * ds;
        t_prev = t;
        t += alphap * dt;
        y += alphad * dy;
    }

    std::tuple<VT, NT, bool> result = std::make_tuple(x, t, converge);
    return result;
}

#endif
