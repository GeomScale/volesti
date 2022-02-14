// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis
// Copyright (c) 2021 Vaibhav Thakkar

// Contributed and/or modified by Vaibhav Thakkar, as part of Google Summer of Code 2021 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef MATH_HELPERS_HPP
#define MATH_HELPERS_HPP

#include <vector>
#include <boost/math/special_functions/gamma.hpp>


//An implementation of Welford's algorithm for mean and variance.
template <typename NT>
std::pair<NT, NT> get_mean_variance(std::vector<NT>& vec)
{
    NT mean = 0;
    NT M2 = 0;
    NT variance = 0;
    NT delta;

    unsigned int i=0;
    for (auto vecit = vec.begin(); vecit!=vec.end(); vecit++, i++)
    {
        delta = *vecit - mean;
        mean += delta / (i + 1);
        M2 += delta * (*vecit - mean);
        variance = M2 / (i + 1);
    }
    return std::pair<NT, NT> (mean, variance);
}


template <typename NT>
static NT log_gamma_function(NT x)
{
    if (x <= NT(100)) return std::log(tgamma(x));
    return (std::log(x - NT(1)) + log_gamma_function(x - NT(1)));
}

#endif  // MATH_HELPERS_HPP