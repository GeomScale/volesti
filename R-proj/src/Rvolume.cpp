// [[Rcpp::depends(BH)]]

// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

#include <Rcpp.h>
#include <RcppEigen.h>
#include "../../include/volume/volume.h"
//#include "vpolyintersectvpoly.h"

//#include <Rcpp.h>



template <class Point, class NT, class Polytope>
double generic_volume(Polytope& P, unsigned int walk_len, double e,
                      Rcpp::Nullable<Rcpp::NumericVector> InnerBall, bool CG, unsigned int win_len,
                      unsigned int N, double C, double ratio, double frac,
                      bool ball_walk, double delta, bool coord, bool rounding)
{
    bool rand_only=false,
         NN=false,
         birk=false,
         verbose =true,
         coordinate=coord;
    unsigned int n_threads=1;

    //unsigned int m;//=A.nrow()-1;
    unsigned int n = P.dimension();//=A.ncol()-1;
    unsigned int rnum = std::pow(e,-2) * 400 * n * std::log(n);

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    // the random engine with this seed
    typedef boost::mt19937    RNGType;
    RNGType rng(seed);
    boost::random::uniform_real_distribution<>(urdist);
    boost::random::uniform_real_distribution<> urdist1(-1,1);

    std::pair<Point,NT> InnerB;

    if(InnerBall.isNotNull()) { //if it is given as an input
        // store internal point hat is given as input
        Rcpp::NumericVector InnerVec = Rcpp::as<Rcpp::NumericVector>(InnerBall);
        std::vector<NT> temp_p;
        for (unsigned int j=0; j<n; j++){
            temp_p.push_back(InnerVec[j]);
        }
        InnerB.first = Point( n , temp_p.begin() , temp_p.end() );
        // store the radius of the internal ball that is given as input
        InnerB.second = InnerVec[n];
    } else {
        // no internal ball or point is given as input
        InnerB = P.ComputeInnerBall();
    }

    // initialization
    vars<NT, RNGType> var(rnum,n,walk_len,n_threads,0.0,e,0,0.0,0, InnerB.second,rng,urdist,urdist1,
                          delta,verbose,rand_only,rounding,NN,birk,ball_walk,coordinate);
    NT vol;
    if (CG) {
        vars<NT, RNGType> var2(rnum, n, 10 + n / 10, n_threads, 0.0, e, 0, 0.0, 0, InnerB.second, rng,
                               urdist, urdist1, delta, verbose, rand_only, rounding, NN, birk, ball_walk, coordinate);
        vars_g<NT, RNGType> var1(n, walk_len, N, win_len, 1, e, InnerB.second, rng, C, frac, ratio, delta, false, verbose,
                                 rand_only, rounding, NN, birk, ball_walk, coordinate);
        vol = volume_gaussian_annealing(P, var1, var2, InnerB);
    } else {
        vol = volume(P, var, var, InnerB);
    }

     return vol;
}

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
double Rvolume (Rcpp::Reference P,  Rcpp::Nullable<unsigned int> walk_len = R_NilValue,
                Rcpp::Nullable<double> error = R_NilValue,
                Rcpp::Nullable<Rcpp::NumericVector> InnerBall = R_NilValue,
                Rcpp::Nullable<std::string> Algo = R_NilValue,
                Rcpp::Nullable<std::string> WalkType = R_NilValue, Rcpp::Nullable<bool> rounding = R_NilValue,
                Rcpp::Nullable<Rcpp::List> Parameters = R_NilValue) {

    typedef double NT;
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef boost::mt19937 RNGType;
    typedef HPolytope <Point> Hpolytope;
    typedef VPolytope <Point, RNGType> Vpolytope;
    typedef Zonotope <Point> zonotope;
    typedef IntersectionOfVpoly<Vpolytope> InterVP;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    int n = P.field("dimension"), walkL;

    bool CG, coordinate, ball_walk, round;
    unsigned int win_len = 4*n*n+500, N = 500 * 2 +  n * n / 2;
    double C = 2.0, ratio = 1.0-1.0/(NT(n)), frac = 0.1, e, delta = -1.0;

    if(!rounding.isNotNull()){
        round = false;
    } else {
        round = Rcpp::as<bool>(rounding);
    }

    if(!WalkType.isNotNull() || Rcpp::as<std::string>(WalkType).compare(std::string("CDHR"))==0){
        coordinate = true;
        ball_walk = false;
    } else if (Rcpp::as<std::string>(WalkType).compare(std::string("RDHR"))==0) {
        coordinate = false;
        ball_walk = false;
    } else if (Rcpp::as<std::string>(WalkType).compare(std::string("BW"))==0) {
        coordinate = false;
        ball_walk = true;
    } else {
        throw std::range_error("Unknown walk type!");
    }

    if(!Algo.isNotNull() || Rcpp::as<std::string>(Algo).compare(std::string("SOB"))==0){

        CG = false;

        if(!walk_len.isNotNull()){
            walkL= 10+n/10;
        } else {
            walkL = Rcpp::as<unsigned int>(walk_len);
        }

        if(!error.isNotNull()){
            e = 1.0;
        } else {
            e = Rcpp::as<NT>(error);
        }

    } else if (Rcpp::as<std::string>(Algo).compare(std::string("CG"))==0) {

        CG = true;

        if (!error.isNotNull()) {
            e = 0.1;
        } else {
            e = Rcpp::as<NT>(error);
        }

        if (!walk_len.isNotNull()) {
            walkL = 1;
        } else {
            walkL = Rcpp::as<int>(walk_len);
        }

    } else {
        throw std::range_error("Unknown method!");
    }

    if(Parameters.isNotNull()) {

        if (Rcpp::as<Rcpp::List>(Parameters).containsElementNamed("BW_rad")) {
            delta = Rcpp::as<NT>(Rcpp::as<Rcpp::List>(Parameters)["BW_rad"]);
        }
        if (Rcpp::as<Rcpp::List>(Parameters).containsElementNamed("C")) {
            C = Rcpp::as<NT>(Rcpp::as<Rcpp::List>(Parameters)["C"]);
            N = 500 * ((int) C) + n * n / 2;
        }
        if (Rcpp::as<Rcpp::List>(Parameters).containsElementNamed("N")) {
            N = Rcpp::as<int>(Rcpp::as<Rcpp::List>(Parameters)["N"]);
        }
        if (Rcpp::as<Rcpp::List>(Parameters).containsElementNamed("Window")) {
            win_len = Rcpp::as<int>(Rcpp::as<Rcpp::List>(Parameters)["Window"]);
        }
        if (Rcpp::as<Rcpp::List>(Parameters).containsElementNamed("frac")) {
            frac = Rcpp::as<NT>(Rcpp::as<Rcpp::List>(Parameters)["frac"]);
        }
        if (Rcpp::as<Rcpp::List>(Parameters).containsElementNamed("ratio")) {
            ratio = Rcpp::as<NT>(Rcpp::as<Rcpp::List>(Parameters)["ratio"]);
        }
    }


    //std::cout<<"n = "<<n<<" Algo null = "<<!Algo.isNotNull()<<" CG = "<<CG<<" N = "<<N<<" C = "<<C<<" win_len = "<<win_len<<" frac = "<<frac<<" ratio = "<<ratio<<" walk_len = "<<walkL<<" error = "<<e<<" coordinate = "<<coordinate<< " ball_walk = "<<ball_walk<<std::endl;

    int type = P.field("type");
    if (type==1) {
        // Hpolytope
        Hpolytope HP;
        HP.init(n, Rcpp::as<MT>(P.field("A")), Rcpp::as<VT>(P.field("b")));
        return generic_volume<Point,NT>(HP, walkL, e, InnerBall, CG, win_len, N, C, ratio, frac, ball_walk, delta,
                                        coordinate, round);
    } else if(type==2) {
        // Vpolytope
        Vpolytope VP;
        VP.init(n, Rcpp::as<MT>(P.field("V")), VT::Ones(Rcpp::as<MT>(P.field("V")).rows()));
        return generic_volume<Point,NT>(VP, walkL, e, InnerBall, CG, win_len, N, C, ratio, frac, ball_walk, delta,
                                        coordinate, round);

    } else if(type==3){
        // Zonotope
        zonotope ZP;
        ZP.init(n, Rcpp::as<MT>(P.field("G")), VT::Ones(Rcpp::as<MT>(P.field("G")).rows()));
        return generic_volume<Point,NT>(ZP, walkL, e, InnerBall, CG, win_len, N, C, ratio, frac, ball_walk, delta,
                                        coordinate, round);
    } else {
        // Intersection of two V-polytopes
        Vpolytope VP1;
        Vpolytope VP2;
        InterVP VPcVP;
        VP1.init(n, Rcpp::as<MT>(P.field("V1")), VT::Ones(Rcpp::as<MT>(P.field("V1")).rows()));
        VP2.init(n, Rcpp::as<MT>(P.field("V2")), VT::Ones(Rcpp::as<MT>(P.field("V2")).rows()));
        VPcVP.init(VP1, VP2);
        return generic_volume<Point,NT>(VPcVP, walkL, e, InnerBall, CG, win_len, N, C, ratio, frac, ball_walk, delta,
                                        coordinate, round);
    }


    return 0;
}
