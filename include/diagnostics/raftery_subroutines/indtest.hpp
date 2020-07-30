// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

//Contributed and/or modified by Alexandros Manochis, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file

//Based on Matlab version of coda package in http://www.spatial-econometrics.com/gibbs/

#ifndef INDTEST_HPP
#define INDTEST_HPP


template<typename MT, typename NT, typename VT>
std::pair<NT,NT> indtest(VT const& d, int const& n)
{
    MT t = MT::Zero(2, 2);
    int t1, t2;
    NT fitted, focus;
    for (int i1 = 1; i1 < n; i1++){
        t(d(i1 - 1), d(i1))=t(d(i1 - 1), d(i1)) + 1;
    }
    NT dcm1 = NT(n) - 1.0,  g2 = 0.0;
    for (int i1 = 0; i1 < 2; i1++){
        for (int i2 = 0; i2 < 2; i2++){
            if (t(i1, i2) != 0){
                t1 = t(i1, 0) + t(i1, 1); 
                t2 = t(0, i2) + t(1, i2);
                fitted = (NT(t1) * NT(t2)) / dcm1; 
                focus = NT(t(i1, i2));
                g2 = g2 + std::log(focus / fitted) * focus;
            }
        } 
    } 
    g2 = g2 * 2.0; 
    NT bic = g2 - std::log(dcm1);
    return std::pair<NT,NT>(g2, bic);
}


#endif
