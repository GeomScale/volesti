// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

//Contributed and/or modified by Alexandros Manochis, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file

//Based on Matlab version of coda package in http://www.spatial-econometrics.com/gibbs/

#ifndef THIN_HPP
#define THIN_HPP

template <typename VT>
std::pair<int,VT>  thin(VT const& work, unsigned int const& n, unsigned int const& kthin)
{
    VT y((n-1) / kthin + 1);

    int i = 0, j = 0;
    while (i < n)
    {
        y(j) = work(i);
        j++;
        i += kthin;
    }

    return std::pair<int, VT>((n-1) / kthin + 1, y);
}


#endif
