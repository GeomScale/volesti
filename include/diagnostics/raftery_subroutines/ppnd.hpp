// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

//Contributed and/or modified by Alexandros Manochis, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file

//Based on Matlab version of coda package in http://www.spatial-econometrics.com/gibbs/

#ifndef PPND_HPP
#define PPND_HPP

template <typename NT>
NT ppnd(NT const& p)
{
    NT split1 = 0.425,  split2 = 5.0, const1 = 0.180625, const2 = 1.6, a0=3.3871327179e+00,
    a1 = 5.0434271938e+01, a2 = 1.5929113202e+02, a3 = 5.9109374720e+01, b1 = 1.7895169469e+01,
    b2 = 7.8757757664e+01, b3 = 6.7187563600e+01, c0 = 1.4234372777e+00, c1 = 2.7568153900e+00,
    c2 = 1.3067284816e+00, c3 = 1.7023821103e-01, d1 = 7.3700164250e-01, d2 = 1.2021132975e-01,
    e0 = 6.6579051150e+00, e1 = 3.0812263860e+00, e2 = 4.2868294337e-01, e3 = 1.7337203997e-02,
    f1 = 2.4197894225e-01, f2 = 1.2258202635e-02; 
    
    NT q = p - 0.5, r, y;
  
    if (std::abs(q) <= split1){
        r = const1 - q * q;
        y = q * (((a3 * r + a2) * r + a1) * r + a0) / (((b3 * r + b2) * r + b1) * r + 1.0);
        return y;
    } else if (q < 0.0){
        r = p; 
    } else{
        r = 1 - p;
    }
    if (r <= 0.0){
        return 0.0;
    }
 
    r = std::sqrt(-1.0 * std::log(r));
 
    if (r <= split2){
        r = r - const2;
        y = (((c3 * r + c2) * r + c1) * r + c0) / ((d2 * r + d1) * r + 1.0);
    } else{
        r = r - split2;
        y = (((e3 * r + e2) * r + e1) * r + e0)/((f2 * r + f1) * r + 1.0);
    }
 
    if (q < 0.0) return -y;

    return y;
}

#endif

