// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

//Contributed and/or modified by Alexandros Manochis, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef EMPQUANT_HPP
#define EMPQUANT_HPP


template <typename VT, typename MT, typename NT>
NT empquant(MT samples, NT q)
{

    //[n junk] = size(runs);  
    
    unsigned int  n = samples.rows(), d = samples.cols();
    VT a(n);
    std::vector<NT> temp_col(n);
    MT work(n, d);

    for (int i = 0; i<d; i++)
    {
        a = samples.col(i);
        temp_col = std::vector<NT>(&a[0], a.data()+a.cols()*a.rows());
        std::sort(temp_col.begin(), temp_col.end());
        work.col(i) = Eigen::Map<VT>(&temp_col[0], temp_col.size());
    }

   

    NT order = (n - 1)*q + 1.0;
    NT fract = order - NT(int(order));// % 1.0;

    int low = std::max(std::floor(order), 1.0);
    int high = std::min(low + 1.0, NT(n));

    NT y = (1.0 - fract) * work(low, 0) + fract*work(high, 0);

    return y;
}


#endif
