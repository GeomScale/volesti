// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

//Contributed and/or modified by Alexandros Manochis, as part of Google Summer of Code 2020 program.

//Licensed under GNU LGPL.3, see LICENCE file

/*
    This function implements a multivariate version of the Rubin & Gelman diagnostic.
    It is based on "Inference from iterative simulation using multiple sequences, 1992" by D. B. Rubin and A. Gelman
    and "General Methods for Monitoring Convergence of Iterative Simulations, 1998" by S. Brooks and A. Gelman

    The sample is splitted into two parts. Then a multivariate psrf is computed as proposed by S. Brooks and A. Gelman

    Inputs: samples, a matrix that contains sample points column-wise

    Output: The value of multivariate PSRF by S. Brooks and A. Gelman
*/

#ifndef DIAGNOSTICS_PSRF_HPP
#define DIAGNOSTICS_PSRF_HPP

template <typename NT, typename VT, typename MT>
NT multivariate_psrf(MT const& samples)
{
    unsigned int N = samples.cols(), d = samples.rows();
    unsigned int N1 = N / 2;
    unsigned int N2 = N - N1;

    VT mean1 = samples.block(0, 0, d, N1).rowwise().mean();
    VT mean2 = samples.block(0, N1, d, N - N1).rowwise().mean();

    MT norm_chain1 = samples.block(0, 0, d, N1).colwise() - mean1;
    MT norm_chain2 = samples.block(0, N1, d, N - N1).colwise() - mean2;

    MT W = ((norm_chain1 * norm_chain1.transpose()) / (NT(N1) - 1.0) +
            (norm_chain2 * norm_chain2.transpose()) / (NT(N2) - 1.0)) / NT(2);

    VT mean00 = (mean1 + mean2) / 2.0;

    MT B = (mean1 - mean00) * (mean1 - mean00).transpose() +
           (mean2 - mean00) * (mean2 - mean00).transpose();

    MT WB = W.inverse() * B;
    Eigen::SelfAdjointEigenSolver <MT> eigensolver(WB);
    NT l_max = eigensolver.eigenvalues().maxCoeff();

    NT R = (NT(N1) - NT(1))/NT(N1) + 1.5 * l_max;
    return R;
}


#endif
