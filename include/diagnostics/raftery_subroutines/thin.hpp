// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

//Contributed and/or modified by Alexandros Manochis, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef THIN_HPP
#define THIN_HPP

template <typename VT, typename MT>
std::pair<int,VT>  thin(MT work, unsigned int n, unsigned int kthin)
{
    VT y((n-1) / kthin + 1);

    int i = 0, j = 0;
    while (i < n)
    {
        y(j) = work(i,0);
        j++;
        i += kthin;
    }

    return std::pair<int,VT>((n-1) / kthin + 1, y);
}



#endif
