// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

//Contributed and/or modified by Alexandros Manochis, as part of Google Summer of Code 2020 program.

//Licensed under GNU LGPL.3, see LICENCE file


#ifndef GEWEKE_HPP
#define GEWEKE_HPP

#include <boost/math/distributions/fisher_f.hpp>

template <typename VT, typename MT, typename NT>
bool perform_geweke(MT const& samples, NT const& frac1, NT const& frac2)
{
    unsigned int d = samples.rows(), N = samples.cols();
    unsigned int N1 = N * frac1;
    unsigned int N2 = N * frac2;

    MT sigma1 = MT::Zero(d, d), sigma2 = MT::Zero(d, d);
    VT mean1 = VT::Zero(d), mean2 = VT::Zero(d);
    NT alpha = 0.05;

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

    MT S_pl = ((NT(N1) - 1.0)*sigma1 + (NT(N2) - 1.0)*sigma2) / (NT(N1) + NT(N2) - 2.0);

    NT T2 = (mean1 - mean2).transpose() * S_pl.inverse() * (mean1 - mean2);
    T2 = ((NT(N1) * NT(N2))/(NT(N1) + NT(N2))) * T2;

    NT U = ((NT(N1) + NT(N2) - NT(d) - 1.0) / ((NT(N1) + NT(N2) - 2.0) * NT(d))) * T2;

    boost::math::fisher_f dist(d, int(N1) + int(N2) - d - 1);
    
    NT F1 = boost::math::quantile(dist, alpha / 2.0);
    NT F2 = boost::math::quantile(boost::math::complement(dist, alpha/2.0));

    if (U <= F1 || U > F2) {
        return false;
    }
    return true;

}

#endif
