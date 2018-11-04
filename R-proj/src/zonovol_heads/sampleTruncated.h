// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

//#include <Rcpp.h>

template <class VT, class MT>
MT sampleTr(VT l, VT u, MT sig, int N, Rcpp::Function rv, MT G){

    MT X2 = G * Rcpp::as<MT>(rv(Rcpp::wrap(l), Rcpp::wrap(u), Rcpp::wrap(sig), N));
    return X2;

}

template <typename NT, class VT, class MT>
NT test_botev(VT l, VT u, MT sig, int N, Rcpp::Function mvNcdf){

    Rcpp::List prob = mvNcdf(Rcpp::wrap(l), Rcpp::wrap(u), Rcpp::wrap(sig), N);

    return prob[0];

}


template <class VT, class MT>
MT sampleTr_gibbs(VT l, VT u, MT sig, int N, int Wstep, Rcpp::Function rtmv, MT G){

    //Rcpp::String algo("gibbs");
    int k = G.cols();
    VT mu = VT::Zero(G.cols());
    int burn = 10*k;
    //std::cout<<u<<"\n"<<l<<std::endl;
    MT X2 = G * Rcpp::as<MT>(rtmv(N, Rcpp::wrap(mu), Rcpp::wrap(sig), Rcpp::wrap(u), Rcpp::wrap(l), burn, Wstep));
    return X2;

}


template <class VT, class MT>
std::pair<MT,MT> sample_cube(VT l, VT u, MT sig, int N, Rcpp::Function rv, MT G){

    MT X2 = Rcpp::as<MT>(rv(Rcpp::wrap(l), Rcpp::wrap(u), Rcpp::wrap(sig), N));
    MT X1 = G * X2;
    return std::pair<MT,MT> (X1,X2);

}


template <class VT, class MT>
MT sampleTr(VT l, VT u, MT sig, int N, Rcpp::Function rv, MT G, int &count){

    int k = sig.cols();
    MT X1 = Rcpp::as<MT>(rv(Rcpp::wrap(l), Rcpp::wrap(u), Rcpp::wrap(sig), N));
    MT X2(k,N);
    VT c(k);
    bool append = false;
    count = 0;

    for (int i = 0; i < N; ++i) {
        c = X1.col(i);
        for (int j = 0; j < k; ++j) {
            if(c(j)>1.0 || c(j)< -1.0) {
                append = true;
                break;
            }
        }
        if (append) {
            X2.col(count) = c;
            count++;
            append = false;
        }
    }

    return G * X2;

}
