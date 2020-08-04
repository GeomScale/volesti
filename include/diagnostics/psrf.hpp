// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

//Contributed and/or modified by Alexandros Manochis, as part of Google Summer of Code 2020 program.

//Licensed under GNU LGPL.3, see LICENCE file

/*
    This function implements a multivariate version of the Rubin & Gelman diagnostic.
    It is based on "Inference from iterative simulation using multiple sequences, 1992" by D. B. Rubin and A. Gelman
    and "General Methods for Monitoring Convergence of Iterative Simulations, 2012" by S. Brooks and A. Gelman

    The sample is splitted into two parts. Then a multivariate psrf is computed as proposed by S. Brooks and A. Gelman
*/

#ifndef PSRF_HPP
#define PSRF_HPP

template <typename NT, typename VT, typename MT>
NT perform_psrf(MT const& points)
{
    unsigned int N = points.cols(), d = points.rows();
    unsigned int N1 = N / 2;
    unsigned int N2 = N - N1;

    MT chain1(d, N1), chain2(d, N2), S1 = MT::Zero(d, d), S2 = MT::Zero(d, d), S(d, d), B = MT::Zero(d, d);
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

    MT SB = S.inverse() * B;
    Eigen::SelfAdjointEigenSolver <MT> eigensolver(SB);
    NT l_max = eigensolver.eigenvalues().maxCoeff();

    NT R = (NT(N1) - NT(1))/NT(N1) + 1.5 * l_max;
    return R;
}


#endif
