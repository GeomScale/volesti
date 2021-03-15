// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

//Contributed and/or modified by Alexandros Manochis, as part of Google Summer of Code 2020 program.

//Licensed under GNU LGPL.3, see LICENCE file

/*
    This function implements a multivariate version of the Geweke diagnostic.
    It is reduced to Hotelling's Two Sample test, which is a multivariate
    extension of the common two sample Student's t-test. The null hypothesis
    is that there is no difference between sample means.

    It is based on "Evaluating the accuracy of sampling-based approaches
                    to the calculation of posterior moments, 1992" by J. Geweke

    Inputs: samples, a matrix that contains sample points column-wise
            frac_first, the portion of the first in order points in matrix samples
            frac_last, the portion of the last in order points in matrix samples
            alpha, the confidence level for the statistical test

    Output: A boolean to denote the result of Geweke diagnostic:
            (i)  false if the null hypothesis is rejected
            (ii) true if the null hypothesis is not rejected
*/


#ifndef DIAGNOSTICS_GEWEKE_HPP
#define DIAGNOSTICS_GEWEKE_HPP

#include <boost/math/distributions/fisher_f.hpp>

template <typename VT, typename MT, typename NT>
bool perform_geweke(MT const& samples,
                    NT frac_first = 0.1,
                    NT frac_last = 0.5,
                    NT alpha = 0.05)
{
    unsigned int d = samples.rows(), N = samples.cols();
    unsigned int N1 = N * frac_first;
    unsigned int N2 = N * frac_last;

    // Compute sample means and covariances
    VT mean1 = samples.block(0, 0, d, N1).rowwise().mean();
    VT mean2 = samples.block(0, N - N2, d, N2).rowwise().mean();

    MT norm_chain1 = samples.block(0, 0, d, N1).colwise() - mean1;
    MT norm_chain2 = samples.block(0, N - N2, d, N2).colwise() - mean2;

    MT sigma1 = (norm_chain1 * norm_chain1.transpose()) / (NT(N1) - 1.0);
    MT sigma2 = (norm_chain2 * norm_chain2.transpose()) / (NT(N2) - 1.0);

    // Compute the pooled covariance matrix
    MT S_pl = ((NT(N1) - NT(1)) * sigma1 + (NT(N2) - 1.0) * sigma2) / (NT(N1) + NT(N2) - NT(2));

    // T2 follows Hotelling's T-squared distribution under the assumption of
    // equal covariances and when the null hypothesis is true
    NT T2 = (mean1 - mean2).transpose() * S_pl.inverse() * (mean1 - mean2);
    T2 = ((NT(N1) * NT(N2)) / (NT(N1) + NT(N2))) * T2;

    // U follows Fischer distribution
    // We use this transformation to check the null hypothesis more easily
    NT U = ((NT(N1) + NT(N2) - NT(d) - 1.0) / ((NT(N1) + NT(N2) - 2.0) * NT(d))) * T2;

    boost::math::fisher_f dist(d, int(N1) + int(N2) - d - 1);

    NT F1 = boost::math::quantile(dist, alpha / 2.0);
    NT F2 = boost::math::quantile(boost::math::complement(dist, alpha / 2.0));

    if (U <= F1 || U > F2) { // reject null hypothesis
        return false;
    }
    return true; // do not reject null hypothesis

}

#endif
