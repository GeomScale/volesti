// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

//Contributed and/or modified by Alexandros Manochis, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef MCEST_HPP
#define MCEST_HPP

template<typename MT, typename NT, typename VT>
std::pair<NT,NT> mcest(VT d, int n)
{
    MT t = MT::Zero(2,2);
    for (int i1 = 1; i1<n; i1++){// 2:n;
        t(d(i1-1)+1,d(i1)+1) = t(d(i1-1)+1,d(i1)+1)+1;
    }//end;
    NT alpha = NT(t(1,2)) / NT((t(1,1)+t(1,2))), beta = NT(t(2,1)) / NT((t(2,1)+t(2,2)));

    return std::pair<NT,NT>(alpha, beta);
}


#endif
