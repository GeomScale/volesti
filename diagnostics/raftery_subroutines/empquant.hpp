// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

//Contributed and/or modified by Alexandros Manochis, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file

//Based on Matlab version of coda package in http://www.spatial-econometrics.com/gibbs/

#ifndef EMPQUANT_HPP
#define EMPQUANT_HPP


template <typename VT, typename NT>
NT empquant(VT const& sorted_samples, NT const& q)
{
    unsigned int n = sorted_samples.rows();

    NT order = (n - 1) * q + 1.0;
    NT fract = order - NT(int(order));
    int low = std::max(round_to_zero(order), 1.0);
    int high = std::min(low + 1.0, NT(n));

    NT y = (1.0 - fract) * sorted_samples(low - 1) + fract * sorted_samples(high - 1);

    return y;
}


#endif
