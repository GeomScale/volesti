// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

//Contributed and/or modified by Alexandros Manochis, as part of Google Summer of Code 2020 program.

//Licensed under GNU LGPL.3, see LICENCE file

/*
    This function implements a multivariate version of the raftery & Lewis diagnostic.
    It is based on Matlab version of coda package in http://www.spatial-econometrics.com/gibbs/
    and "How many iterations in the Gibbs sampler?, 1992" by A. Raftery and S. Lewis

    Inputs: samples, a matrix that contains sample points column-wise
            q, the quantile of the quantity of interest. The default value is 0.025.
            r, the level of precision desired. The default value is 0.01.
            s, the probability associated with r. The default value is 0.95.

    Outputs: (i)   The number of draws required for burn-in
             (ii)  The skip parameter for 1st-order Markov chain
             (iii) The skip parameter sufficient to get independence chain
             (iv)  The number of draws required to achieve r precision
             (v)   The number of draws if the chain is white noise
             (vi)  The I-statistic from Raftery and Lewis (1992)
*/

#ifndef DIAGNOSTICS_RAFTERY_HPP
#define DIAGNOSTICS_RAFTERY_HPP

template <typename NT>
NT round_to_zero(NT x)
{
    return (x > 0.0) ? std::floor(x) : std::ceil(x);
}

#include "raftery_subroutines/empquant.hpp"
#include "raftery_subroutines/indtest.hpp"
#include "raftery_subroutines/mctest.hpp"
#include "raftery_subroutines/mcest.hpp"
#include "raftery_subroutines/thin.hpp"
#include "raftery_subroutines/ppnd.hpp"


template <typename VT, typename MT, typename NT>
MT perform_raftery(MT const& samples, NT const& q, NT const& r, NT const& s)
{
    MT runs = samples.transpose();

    typedef Eigen::Matrix<int,Eigen::Dynamic,Eigen::Dynamic> MTint;
    typedef Eigen::Matrix<int,Eigen::Dynamic,1> VTint;

    unsigned int n = runs.rows(), d = runs.cols(), kthin, kmind;
    MT results(d, 6);
    MTint work = MTint::Zero(n, d);
    VTint tmp = VTint::Zero(n);
    std::pair<int, VTint> xy;
    std::pair<NT, NT> g2bic;

    NT cutpt, alpha, beta, g2, bic, epss;
    int tcnt;

    MT sorted_samples(n, d);
    VT a(n);
    std::vector<NT> temp_col(n);

    for (int i = 0; i < d; i++)
    {
        a = runs.col(i);
        temp_col = std::vector<NT>(&a[0], a.data() + a.cols() * a.rows());
        std::sort(temp_col.begin(), temp_col.end());
        sorted_samples.col(i) = Eigen::Map<VT>(&temp_col[0], temp_col.size());
    }

    for (int i = 0; i < d; i++)
    {
        cutpt = empquant<VT>(sorted_samples.col(i), q);
        for (int j = 0; j < n; j++)
        {
            if (runs(j, i) <= cutpt) work(j, i) = 1;
        }
        kthin = 1; bic = 1.0; epss = 0.001;

        while(bic > 0.0)
        {
            xy = thin<VTint>(work.col(i), n, kthin);
            tcnt = xy.first;
            tmp = xy.second;
            g2bic = mctest<MTint, NT>(tmp, tcnt);
            g2 = g2bic.first;
            bic = g2bic.second;
            kthin++;
            if (kthin > n / 2) {
                break;
            }
        }

        kthin--;
        g2bic = mcest<MTint, NT>(tmp, tcnt);
        alpha = g2bic.first;
        beta = g2bic.second;
        kmind = kthin;
        g2bic = indtest<MTint, NT>(tmp, tcnt);
        g2 = g2bic.first;
        bic = g2bic.second;

        while (bic > 0.0)
        {
            xy = thin<VTint>(work.col(i), n, kmind);
            tcnt = xy.first;
            tmp = xy.second;
            g2bic = indtest<MTint, NT>(tmp, tcnt);
            g2 = g2bic.first;
            bic = g2bic.second;
            kmind++;
            if (kmind > n) {
                break;
            }
        }

        NT psum = alpha + beta;
        NT tmp1  = std::log(psum * epss / std::max(alpha, beta)) / std::log(std::abs(1.0 - psum));
        NT nburn = round_to_zero((tmp1 + 1.0) * NT(kthin));
        NT phi   = ppnd((s + 1.0) / 2.0);
        NT tmp2  = (2.0 - psum) * alpha * beta * (phi * phi) / (psum * psum * psum * (r * r));
        NT nprec = round_to_zero(tmp2 + 1.0) * kthin;
        NT nmin  = round_to_zero(((1.0 - q) * q * (phi * phi) / (r * r)) + 1.0);
        NT irl   = (nburn + nprec) / nmin;
        NT kind  = std::max(round_to_zero(irl + 1.0), NT(kmind));

        results(i, 0) = NT(kthin);
        results(i, 1) = NT(nburn);
        results(i, 2) = kind;
        results(i, 3) = NT(nburn) + nprec;
        results(i, 4) = nmin;
        results(i, 5) = irl;
    }
    return results;
}


#endif
