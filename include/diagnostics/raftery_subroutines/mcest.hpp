// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

//Contributed and/or modified by Alexandros Manochis, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file

//Based on Matlab version of coda package in http://www.spatial-econometrics.com/gibbs/

#ifndef MCEST_HPP
#define MCEST_HPP

template<typename MT, typename NT, typename VT>
std::pair<NT,NT> mcest(VT const& d, int const& n)
{
    MT t = MT::Zero(2, 2);
    for (int i1 = 1; i1 < n; i1++){
        t(d(i1 - 1), d(i1)) = t(d(i1 - 1), d(i1)) + 1;
    }
    NT alpha = NT(t(0, 1)) / NT((t(0, 0) + t(0, 1))), beta = NT(t(1, 0)) / NT((t(1, 0)+t(1, 1)));

    return std::pair<NT, NT>(alpha, beta);
}


#endif
