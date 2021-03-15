// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

//Contributed and/or modified by Alexandros Manochis, as part of Google Summer of Code 2020 program.

//Licensed under GNU LGPL.3, see LICENCE file

/*
    This function implements the Rubin & Gelman diagnostic.
    It is based on "Inference from iterative simulation using multiple sequences, 1992" by D. B. Rubin and A. Gelman

    For each coordinate the sample is splitted into two parts.
    Then the psrf of D.B. Rubin and A. Gelman is computed for each coordinate

    Inputs: samples, a matrix that contains sample points column-wise

    Output: The value of PSRF of D.B. Rubin and A. Gelman for each coordinate
*/

#ifndef DIAGNOSTICS_MARGINAL_PSRF_HPP
#define DIAGNOSTICS_MARGINAL_PSRF_HPP

template <typename NT, typename VT, typename MT>
VT univariate_psrf(MT const& samples)
{
    MT runs = samples.transpose();
    unsigned int N = samples.cols(), d = samples.rows();
    unsigned int N1 = N / 2;
    unsigned int N2 = N - N1;
    VT coord_samples(N), results(d);
    NT mean1, mean2, mean00, sum, R, W, B, sigma;

    for (int i = 0; i < d; i++)
    {
        coord_samples = runs.col(i);
        mean1 = coord_samples.block(0,0,N1,1).mean();
        mean2 = coord_samples.block(N1,0,N2,1).mean();

        sum = NT(0);
        for (int j = 0; j < N1; j++)
        {
            sum += (coord_samples(j) - mean1) * (coord_samples(j) - mean1);
        }
        W = sum / (NT(N1) - NT(1));

        sum = NT(0);
        for (int j = N1; j < N; j++)
        {
            sum += (coord_samples(j) - mean2) * (coord_samples(j) - mean2);
        }
        W += (sum / (NT(N2) - NT(1)));
        W = W / NT(2);

        mean00 = coord_samples.mean();

        B = (mean1 - mean00) * (mean1 - mean00) + (mean2 - mean00) * (mean2 - mean00);
        sigma = ((NT(N1) - NT(1)) / NT(N1)) * W + B;
        R = std::sqrt(sigma / W);

        results(i) = R;
    }
    return results;
}

#endif
