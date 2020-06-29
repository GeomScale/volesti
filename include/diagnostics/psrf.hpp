//Contributed and/or modified by Alexandros Manochis, as part of Google Summer of Code 2020 program.

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


#ifndef PSRF_HPP
#define PSRF_HPP

template <typename NT, typename VT, typename MT>
NT perform_psrf(MT const& points)
{
    //typedef double NT;

    unsigned int N = points.cols(), d = points.rows();
    unsigned int N1 = N/2;
    unsigned int N2 = N - N1;

    MT chain1(d, N1), chain2(d, N2), S1 = MT::Zero(d,d), S2 = MT::Zero(d,d), S(d,d), B = MT::Zero(d,d);
    VT mean1 = VT::Zero(d), mean2 = VT::Zero(d);

    for (int i = 0; i < N1; ++i) {
        mean1 += points.col(i);
        chain1.col(i) = points.col(i);
    }

    for (int i = 0; i < N2; ++i) {
        mean2 += points.col(i + N1);
        chain2.col(i) = points.col(i + N1);
    }

    mean1 = mean1 / NT(N1); mean2 = mean2 / NT(N2);

    for (int i = 0; i < N1; ++i) {
        S1 = S1 + (chain1.col(i) - mean1) * (chain1.col(i) - mean1).transpose();
    }

    for (int i = 0; i < N2; ++i) {
        S2 = S2 + ((chain2.col(i) - mean2) * (chain2.col(i) - mean2).transpose());
    }
    S1 = S1 / (NT(N1) - 1.0); S2 = S2 / (NT(N2) - 1.0);

    S = (S1 + S2) / 2.0;

    VT mean00 = (mean1 + mean2) / 2.0;

    B = B + (mean1 - mean00) * (mean1 - mean00).transpose();
    B = B + (mean2 - mean00) * (mean2 - mean00).transpose();
    B = B * NT(N/2.0);

    MT SB = S.inverse() * B;
    Eigen::SelfAdjointEigenSolver <MT> eigensolver(SB);
    //rel = eigensolver.eigenvalues().minCoeff();
    NT l_max = eigensolver.eigenvalues().maxCoeff();


    NT R = std::sqrt((NT(N1) - 1.0)/NT(N1) + 1.5*(l_max)/NT(N1));

    return R;

}


#endif
