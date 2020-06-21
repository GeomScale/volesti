//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.
//Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.

// VolEsti is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or (at
// your option) any later version.
//
// VolEsti is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
//
// See the file COPYING.LESSER for the text of the GNU Lesser General
// Public License.  If you did not receive this file along with HeaDDaCHe,
// see <http://www.gnu.org/licenses/>.


#ifndef MVE_COMPUTATION_HPP
#define MVE_COMPUTATION_HPP


template <typename MT, typename VT, typename NT>
std::pair<std::pair<MT, VT>, bool> mve_computation(MT A, VT b, VT const& x0,
             unsigned int const& maxiter, NT const& tol, NT const& reg)
{
    int m = A.rows(), n = A.cols();

    bool converged = false;

    NT bnrm = b.norm(), minmu = std::pow(10.0, -8.0), tau0 = 0.75, tau, last_r1 = std::numeric_limits<NT>::lowest(),
          last_r2 = std::numeric_limits<NT>::lowest(), prev_obj = std::numeric_limits<NT>::lowest(), *vec_iter1,
          *vec_iter2, *vec_iter3, gap, rmu, res, objval, r1, r2 ,r3, rel, Rel, astep, ax, ay, az,
          reg_lim = std::pow(10.0, -10.0);

    VT bmAx0 = b - A * x0, x = VT::Zero(n), y = VT::Ones(m), bmAx = VT::Ones(m), h(m), z(m), yz(m), yh(m), R1(n), R2(m),
       R3(m), y2h(m), y2h_z(m), h_z(m), R3Dy(m), R23(m), dx(n), Adx(m), dyDy(m), dy(m), dz(m),
       ones_m = VT::Ones(m);
    MT Q(m, m), Y(m, m), E2(n, n), YQ(m,m), YA(m, n), G(m,m), T(m,n), ATP(n,m), ATP_A(n,n);

    A = (ones_m.cwiseProduct(bmAx0.cwiseInverse())).asDiagonal() * A, b = ones_m;

    MT A_trans = A.transpose();

    int i = 1;
    while (i <= maxiter) {

        Y = y.asDiagonal();
        E2.noalias() = (A_trans * Y * A).inverse();

        Q.noalias() = A * E2 * A_trans;
        h = Q.diagonal();
        h = h.cwiseSqrt();

        if (i == 1) {
            NT t = bmAx.cwiseProduct(h.cwiseInverse()).minCoeff();
            y *= (1.0 / (t * t));
            h *= t;
            vec_iter1 = bmAx.data();
            vec_iter2 = h.data();
            vec_iter3 = z.data();
            for (int j = 0; j < m; ++j) {
                *vec_iter3 = std::max(0.1, (*vec_iter1 - (*vec_iter2)));
                vec_iter1++;
                vec_iter2++;
                vec_iter3++;
            }
            Q *= (t * t);
            Y *= (1.0 / (t * t));
        }

        yz = y.cwiseProduct(z);
        yh = y.cwiseProduct(h);

        gap = yz.sum() / NT(m);
        rmu = std::min(0.5, gap) * gap;
        rmu = std::max(rmu, minmu);

        R1.noalias() = - A_trans * yh;
        R2 = bmAx - h - z;
        R3.noalias() = rmu * ones_m - yz;

        r1 = R1.template lpNorm<Eigen::Infinity>();
        r2 = R2.template lpNorm<Eigen::Infinity>();
        r3 = R3.template lpNorm<Eigen::Infinity>();

        res = std::max(r1, r2);
        res = std::max(res, r3);
        objval = std::log(E2.determinant()) / 2.0;

        Eigen::SelfAdjointEigenSolver <MT> eigensolver(E2);
        rel = eigensolver.eigenvalues().minCoeff();
        Rel = eigensolver.eigenvalues().maxCoeff();

        if (i % 10 == 0) {

            if (std::abs((last_r1 - r1) / std::min(NT(std::abs(last_r1)), NT(std::abs(r1)))) < 0.01 &&
                std::abs((last_r2 - r2) / std::min(NT(abs(last_r2)), NT(std::abs(r2)))) < 0.01 &&
                Rel / rel > 100.0 &&
                reg > reg_lim) {
                converged = false;
                //std::cout << "Stopped making progress, stopping and restarting." << std::endl;
                break;
            }
            last_r2 = r2;
            last_r1 = r1;
        }

        if ((res < tol * (1.0 + bnrm) && rmu <= minmu) || (i > 100 && prev_obj != std::numeric_limits<NT>::lowest() &&
                                                           (prev_obj >= (1.0 - tol) * objval ||
                                                            objval <= (1.0 - tol) * prev_obj))) {
            //std::cout<<"iteration = "<<i<<std::endl;
            //std::cout << "converged!" << std::endl;
            x += x0;
            converged = true;
            break;
        }

        prev_obj = objval;
        YQ.noalias() = Y * Q;
        G = YQ.cwiseProduct(YQ.transpose());
        y2h = 2.0 * yh;
        YA.noalias() = Y * A;

        vec_iter1 = y2h.data();
        vec_iter2 = z.data();
        vec_iter3 = y2h_z.data();
        for (int j = 0; j <m; ++j) {
            *vec_iter3 = std::max(reg, (*vec_iter1) * (*vec_iter2));
            vec_iter1++;
            vec_iter2++;
            vec_iter3++;
        }

        G.diagonal() += y2h_z;
        h_z = h + z;

        for (int j = 0; j < n; ++j) {
            T.col(j) = G.colPivHouseholderQr().solve(YA.col(j).cwiseProduct(h_z));
        }
        ATP.noalias() = (y2h.asDiagonal()*T - YA).transpose();

        vec_iter1 = R3.data();
        vec_iter2 = y.data();
        vec_iter3 = R3Dy.data();
        for (int j = 0; j <m; ++j) {
            *vec_iter3 = (*vec_iter1) / (*vec_iter2);
            vec_iter1++;
            vec_iter2++;
            vec_iter3++;
        }

        R23 = R2 - R3Dy;
        ATP_A.noalias() = ATP * A;
        ATP_A.diagonal() += ones_m * reg;
        dx = ATP_A.colPivHouseholderQr().solve(R1 + ATP * R23);

        Adx.noalias() = A * dx;
        dyDy = G.colPivHouseholderQr().solve(y2h.cwiseProduct(Adx-R23));

        dy = y.cwiseProduct(dyDy);
        dz = R3Dy - z.cwiseProduct(dyDy);

        vec_iter1 = Adx.data();
        vec_iter2 = bmAx.data();
        ax = -0.5;
        for (int j = 0; j <m; ++j) {
            ax = std::min(ax, - (*vec_iter1) / (*vec_iter2));
            vec_iter1++;
            vec_iter2++;
        }

        ax = -1.0 / ax;
        ay = -1.0 / std::min(dyDy.minCoeff(), -0.5);

        vec_iter1 = dz.data();
        vec_iter2 = z.data();
        az = -0.5;
        for (int j = 0; j <m; ++j) {
            az = std::min(az, (*vec_iter1) / (*vec_iter2));
            vec_iter1++;
            vec_iter2++;
        }

        az = -1.0 / az;
        tau = std::max(tau0, 1.0 - res);
        astep = tau*std::min(std::min(1.0, ax), std::min(ay, az));

        x += astep*dx;
        y += astep*dy;
        z += astep*dz;

        bmAx -= astep * Adx;

        if (reg > std::pow(10,-6.0) && i >= 10){
            break;
        }

        i++;
    }

    return std::pair<std::pair<MT, VT>, bool>(std::pair<MT, VT>(E2, x), converged);

}


#endif
