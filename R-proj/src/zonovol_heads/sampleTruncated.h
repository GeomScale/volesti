// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

#include <Rcpp.h>

template <class PointList>
PointList sampleTr(Rcpp::NumericVector l, Rcpp::NumericVector u,
                   Rcpp::NumericMatrix sig, Rcpp::Function rv, PointList &randPoints){

    Rcpp::NumericMatrix X = rv(l, u, sig, 100);
    //std::vector<std::vector<double> > Pin(std::vector<std::vector<double> >::iterator(X.begin()), std::vector<std::vector<double> >::iterator(X.begin()) );
    //std::cout<<Pin.size()<<std::endl;
    randPoints.clear();
    return randPoints;

}