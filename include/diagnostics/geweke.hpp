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


#ifndef GEWEKE_HPP
#define GEWEKE_HPP

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

    MT sigma1 = MT::Zero(d, d), sigma2 = MT::Zero(d, d);
    VT mean1 = VT::Zero(d), mean2 = VT::Zero(d);

    // Compute sample means and covariances
    for (int i = 0; i < N1; ++i) {
        mean1 += samples.col(i);
    }
    mean1 = mean1 / NT(N1);

    for (int i = 0; i < N1; ++i) {
        sigma1 = sigma1 + (samples.col(i) - mean1) * (samples.col(i) - mean1).transpose();
    }
    sigma1 = sigma1 / (NT(N1) - 1.0);


    for (int i = N-N2; i < N; ++i) {
        mean2 += samples.col(i);
    }
    mean2 = mean2 / NT(N2);

    for (int i = N-N2; i < N; ++i) {
        sigma2 = sigma2 + (samples.col(i) - mean2) * (samples.col(i) - mean2).transpose();
    }
    sigma2 = sigma2 / (NT(N2) - 1.0);

    // Compute the pooled covariance matrix
    MT S_pl = ((NT(N1) - 1.0)*sigma1 + (NT(N2) - 1.0)*sigma2) / (NT(N1) + NT(N2) - 2.0);

    // T2 follows Hotelling's T-squared distribution under the assumption of 
    // equal covariances and when the null hypothesis is true
    NT T2 = (mean1 - mean2).transpose() * S_pl.inverse() * (mean1 - mean2);
    T2 = ((NT(N1) * NT(N2))/(NT(N1) + NT(N2))) * T2;

    // U follows Fischer distribution
    // We use this transformation to check the null hypothesis more easily
    NT U = ((NT(N1) + NT(N2) - NT(d) - 1.0) / ((NT(N1) + NT(N2) - 2.0) * NT(d))) * T2;

    boost::math::fisher_f dist(d, int(N1) + int(N2) - d - 1);
    
    NT F1 = boost::math::quantile(dist, alpha / 2.0);
    NT F2 = boost::math::quantile(boost::math::complement(dist, alpha/2.0));

    if (U <= F1 || U > F2) { // reject null hypothesis
        return false;
    }
    return true; // do not reject null hypothesis

}

#endif
