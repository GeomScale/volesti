// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef ZONOTOPE_EXACT_VOL_H
#define ZONOTOPE_EXACT_VOL_H

#include <algorithm>
#include <iostream>
#include <string>

// From rosetta code at http://rosettacode.org/wiki/Combinations#C.2B.2B
// We made some adjustments to vectorize the output
// Compute all the N combinations from N elements
inline std::vector< std::vector<int> > comb(int N, int K)
{
    std::string bitmask(K, 1); // K leading 1's
    bitmask.resize(N, 0); // N-K trailing 0's
    std::vector<int> iter_comb(K,0);
    std::vector<std::vector<int> > combs;
    int count;

    // print integers and permute bitmask
    do {
        count = 0;
        for (int i = 0; i < N; ++i) // [0..N-1] integers
        {
            if (bitmask[i]){
                iter_comb[count] = i;
                count++;
            }
        }
        combs.push_back(iter_comb);
    } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
    return combs;
}


template <typename NT, typename Polytope>
NT exact_zonotope_vol(const Polytope &ZP){

    typedef typename Polytope::MT 	MT;
    typedef std::vector< std::vector<int> >::iterator  IntMatIt;
    typedef std::vector<int>::iterator  IntIt;

    int n = ZP.dimension(), col, k = ZP.num_of_generators();
    MT V1 = ZP.get_mat().transpose(), SubV(n,n), V(n, 2*k);
    V << V1, -V1;
    NT vol = 0.0;

    std::vector< std::vector<int> > combs = comb(2*k, n);
    IntMatIt iter_combs;
    IntIt it;

    iter_combs = combs.begin();
    for ( ;  iter_combs!=combs.end(); ++iter_combs) {
        it = (*iter_combs).begin();
        col = 0;
        // construct all the nxn submatrices
        for ( ; it!=(*iter_combs).end(); ++it, ++col) {
            SubV.col(col) = V.col(*it);
        }
        vol += std::abs(SubV.determinant());
    }
    return vol;
}

template <typename NT>
NT vol_Ali(std::vector<NT> &plane, const NT &zit, const unsigned int dim) {

    unsigned int i, J = 0, counter = 0, K = 0, k;
    std::vector <NT> Y(dim + 2, 0.0), X(dim + 2, 0.0), a(dim + 2, 0.0);

    if (zit < 0) {
        X[0] = zit;
        J++;
    } else {
        Y[0] = zit;
        counter++;
    }

    for (i = 0; i < dim; i++) {

        a[i] = 0.0;

        if (plane[i] + zit < 0) {
            X[J] = plane[i] + zit;
            J++;
        } else {
            Y[counter] = plane[i] + zit;
            counter++;
        }
    }
    K = dim + 1 - J;
    a[0] = 1.0;
    a[dim] = 0.0;
    a[dim + 1] = 0.0;

    for (i = 0; i < J; i++) {
        for (k = 1; k < K + 1; k++) {
            a[k] = (Y[k - 1] * a[k] - X[i] * a[k - 1]) / (Y[k - 1] - X[i]);
        }
    }
    return a[K];
}

#endif
