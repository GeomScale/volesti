// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

//#include <Rcpp.h>

template <class VT, class MT>
MT sampleTr(VT l, VT u, MT sig, int N, Rcpp::Function rv, MT G){

    //arma::mat X = Rcpp::as<arma::mat>(rv(l, u, sig, N));
    MT X2 = G * Rcpp::as<MT>(rv(l, u, sig, N));
    return X2;

}