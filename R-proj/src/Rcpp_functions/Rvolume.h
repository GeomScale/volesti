// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

#ifndef RVOLUME_H
#define RVOLUME_H

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
double Rvolume (Rcpp::NumericMatrix A, unsigned int walk_len, double e,
                Rcpp::NumericVector InnerVec, bool CG, unsigned int win_len,
                unsigned int N, double C, double ratio, double frac,
                bool ball_walk, double delta, bool Vpoly, bool Zono,
                bool coord, bool rounding) {

    bool rand_only=false,
            NN=false,
            birk=false,
            verbose =false,
            coordinate=coord;
    Hpolytope HP;
    Vpolytope VP;
    Zonotope ZP;

    unsigned int m=A.nrow()-1;
    unsigned int n=A.ncol()-1;
    unsigned int rnum = std::pow(e,-2) * 400 * n * std::log(n);
    Rcpp::NumericMatrix vol_res(1,1);

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    // the random engine with this seed
    RNGType rng(seed);
    boost::random::uniform_real_distribution<>(urdist);
    boost::random::uniform_real_distribution<> urdist1(-1,1);

    std::vector<std::vector<NT> > Pin(m+1, std::vector<NT>(n+1));

    for (i=0; i<m+1; i++){
        for(j=0; j<n+1; j++){
            Pin[i][j]=A(i,j);
        }
    }
    // construct polytope
    if (Zono) {
        ZP.init(Pin);
        if (exact_zono) {
            vol_res(0,0) = exact_zonotope_vol<NT>(ZP);
            return vol_res;
        }
    } else if (!Vpoly) {
        HP.init(Pin);
    } else {
        VP.init(Pin);
    }

    std::pair<Point,NT> InnerBall;

    if(InnerVec.size()==n+1) { //if it is given as an input
        // store internal point hat is given as input
        std::vector<NT> temp_p;
        for (unsigned int j=0; j<n; j++){
            temp_p.push_back(InnerVec[j]);
        }
        InnerBall.first = Point( n , temp_p.begin() , temp_p.end() );
        // store the radius of the internal ball that is given as input
        InnerBall.second = InnerVec[n];
    } else {
        // no internal ball or point is given as input
        if (Zono) {
            InnerBall = ZP.ComputeInnerBall();
        } else if (!Vpoly) {
            InnerBall = HP.ComputeInnerBall();
        } else {
            InnerBall = VP.ComputeInnerBall();
        }
    }


    // initialization
    vars<NT, RNGType> var(rnum,n,walk_len,n_threads,0.0,0.0,0,0.0,0, InnerBall.second,rng,urdist,urdist1,
                          delta,verbose,rand_only,rounding,NN,birk,ball_walk,coordinate);
    NT vol;
    if (CG) {
        vars<NT, RNGType> var2(rnum, n, 10 + n / 10, n_threads, 0.0, e, 0, 0.0, 0, InnerBall.second, rng,
                               urdist, urdist1, delta, verbose, rand_only, rounding, NN, birk, ball_walk, coordinate);
        vars_g<NT, RNGType> var1(n, walk_len, N, win_len, 1, e, InnerBall.second, rng, C, frac, ratio, delta, false, verbose,
                                 rand_only, rounding, NN, birk, ball_walk, coordinate);
        if (Zono) {
            vol = volume_gaussian_annealing(ZP, var1, var2, InnerBall);
        } else if (!Vpoly) { // if the input is a H-polytope
            vol = volume_gaussian_annealing(HP, var1, var2, InnerBall);
        } else {  // if the input is a V-polytope
            vol = volume_gaussian_annealing(VP, var1, var2, InnerBall);
        }
    } else {
        if (Zono) {
            vol = volume(ZP, var, var, InnerBall);
        } else if (!Vpoly) { // if the input is a H-polytope
            vol = volume(HP, var, var, InnerBall);
        } else { // if the input is a V-polytope
            vol = volume(VP, var, var, InnerBall);
        }
    }

    return vol;

}

#endif
