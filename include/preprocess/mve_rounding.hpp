//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

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


#ifndef MVE_ROUNDING_HPP
#define MVE_ROUNDING_HPP

template <
        typename MT,
        typename VT,
        typename Polytope,
        typename Point,
        typename NT,
        typename RandomNumberGenerator
>
std::pair< std::pair<MT, VT>, NT > mve_rounding(Polytope &P, std::pair<Point, NT> InnerBall,
                                                RandomNumberGenerator &rng)
{
    std::pair<std::pair<MT, VT>, bool> iter_res;
    iter_res.second = false;

    VT x0 = InnerBall.first.getCoefficients();//, b = P.get_vec();
    MT E, L;
    unsigned int maxiter = 150, iter = 1;

    NT R = 100.0, r = 1.0, tol = std::pow(10, -6.0), reg = std::pow(10, -4.0), round_val = 1.0;

    MT T = MT::Identity(d,d);
    VT shift = VT::Zero(d);

    while ((R > 6.0 * r && iter < 20) || !iter_res.second )
    {
        iter_res = mve_computation(P.get_mat(), P.get_vec(), x0, maxiter, tol, reg);
        E = iter_res.first.first;

        Eigen::LLT<MT> lltOfA(E); // compute the Cholesky decomposition of E
        L = lltOfA.matrixL();

        P.shift(iter_res.first.second);
        MT L_1 = L.inverse();
        shift += T * iter_res.first.second;
        T = T * L_1.transpose();
        round_val *= L_1.determinant();

        reg = std::max(reg / 10.0, std::pow(10, -10.0));
        P.normalize();

        Eigen::EigenSolver<MT> eigensolver(E);
        r = std::real(eigensolver.eigenvalues()[0]);
        R = std::real(eigensolver.eigenvalues()[0]);
        for(i=1; i<d; i++) {
            if (std::real(eigensolver.eigenvalues()[i]) < r) r = std::real(eigensolver.eigenvalues()[i]);
            if (std::real(eigensolver.eigenvalues()[i]) > R) R = std::real(eigensolver.eigenvalues()[i]);
        }
        x0 = VT::Zero(P.dimension());
    }

    return std::pair< std::pair<MT, VT>, NT > (std::pair<MT, VT>(T, shift), round_val);
}

#endif
