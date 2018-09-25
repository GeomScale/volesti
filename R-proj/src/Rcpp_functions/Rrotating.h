// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

#ifndef RROTATING_H
#define RROTATING_H

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
Rcpp::NumericMatrix rotating (Rcpp::NumericMatrix A) {

    unsigned int m = A.nrow() - 1;
    unsigned int n = A.ncol() - 1;
    std::vector <std::vector<NT>> Pin(m + 1, std::vector<NT>(n + 1));

    for (unsigned int i = 0; i < m + 1; i++) {
        for (unsigned int j = 0; j < n + 1; j++) {
            Pin[i][j] = A(i, j);
        }
    }
    // construct polytope
    if (Zono) {
        ZP.init(Pin);
    } else if (!Vpoly) {
        HP.init(Pin);
    } else {
        VP.init(Pin);
    }

    Rcpp::NumericMatrix Mat;
    if (Zono) {
        rotating<NT>(ZP);
        Mat = extractMatPoly(ZP);
    }else if (!Vpoly) {
        rotating<NT>(HP);
        Mat = extractMatPoly(HP);
    } else {
        rotating<NT>(VP);
        Mat = extractMatPoly(VP);
    }
    return Mat;

}

#endif
