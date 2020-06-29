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


#ifndef GEWEKE_HPP
#define GEWEKE_HPP

#include <boost/math/distributions/fisher_f.hpp>

template <typename VT, typename MT, typename NT>
bool perform_geweke(MT const& samples, NT const& frac1, NT const& frac2)
{
    unsigned int d = samples.rows(), N = samples.cols();
    unsigned int N1 = N*frac1;
    unsigned int N2 = N*frac2;

    MT sigma1 = MT::Zero(d,d), sigma2 = MT::Zero(d,d);
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
    T2 = ((NT(N1)*NT(N2))/(NT(N1) + NT(N2))) * T2;

    NT U = ((NT(N1) + NT(N2) - NT(d) - 1.0) / ((NT(N1) + NT(N2) - 2.0)*NT(d))) * T2;

    boost::math::fisher_f dist(d, N1 + N2 - d - 1);

    if (U > (boost::math::quantile(boost::math::complement(dist, alpha)))) {
        return false;
    }
    return true;

}

#endif
