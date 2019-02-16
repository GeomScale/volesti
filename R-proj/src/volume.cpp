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
double generic_volume(Polytope& P, unsigned int walk_length, double e,
                      Rcpp::Nullable<Rcpp::NumericVector> InnerBall, bool CG, unsigned int win_len,
                      unsigned int N, double C, double ratio, double frac,
                      bool ball_walk, double delta, bool coord, bool rounding)
{
    bool rand_only=false,
         NN=false,
         birk=false,
         verbose =false,
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
    vars<NT, RNGType> var(rnum,n,walk_length,n_threads,0.0,e,0,0.0,0, InnerB.second,rng,urdist,urdist1,
                          delta,verbose,rand_only,rounding,NN,birk,ball_walk,coordinate);
    NT vol;
    if (CG) {
        vars<NT, RNGType> var2(rnum, n, 10 + n / 10, n_threads, 0.0, e, 0, 0.0, 0, InnerB.second, rng,
                               urdist, urdist1, delta, verbose, rand_only, rounding, NN, birk, ball_walk, coordinate);
        vars_g<NT, RNGType> var1(n, walk_length, N, win_len, 1, e, InnerB.second, rng, C, frac, ratio, delta, false, verbose,
                                 rand_only, rounding, NN, birk, ball_walk, coordinate);
        vol = volume_gaussian_annealing(P, var1, var2, InnerB);
    } else {
        vol = volume(P, var, var, InnerB);
    }

     return vol;
}

//' The main function for volume approximation of a convex Polytope (H-polytope, V-polytope or a zonotope)
//'
//' For the volume approximation can be used two algorithms. Either SequenceOfBalls or CoolingGaussian. A H-polytope with \eqn{m} facets is described by a \eqn{m\times d} matrix \eqn{A} and a \eqn{m}-dimensional vector \eqn{b}, s.t.: \eqn{Ax\leq b}. A V-polytope is described as a set of \eqn{d}-dimensional points. A zonotope is desrcibed by the Minkowski sum of \eqn{d}-dimensional segments.
//'
//' @param P A convex polytope. It is an object from class (a) Hpolytope or (b) Vpolytope or (c) Zonotope.
//' @param walk_length Optional. The number of the steps for the random walk. Default value is \eqn{\lfloor 10 + d/10\rfloor} for SequenceOfBalls and \eqn{1} for CoolingGaussian.
//' @param error Optional. Declare the upper bound for the approximation error. Default value is \eqn{1} for SequenceOfBalls and \eqn{0.1} for CoolingGaussian.
//' @param InnerBall Optional. A \eqn{d+1} vector that containes an inner ball. The first \eqn{d} coordinates corresponds to the center and the last one to the radius of the ball. If it is not given then for H-polytopes the Chebychev ball is computed, for V-polytopes \eqn{d+1} vertices are picked randomly and the Chebychev ball of the defined simplex is computed. For a zonotope that is defined by the Minkowski sum of \eqn{m} segments we compute the maximal \eqn{r} s.t.: \eqn{re_i\in Z} for all \eqn{i=1,\dots ,d}, then the ball centered at the origin with radius \eqn{r/\sqrt{d}} is an inscribed ball.
//' @param Algo Optional. A string that declares which algorithm to use: a) 'SoB' for SequenceOfBalls or b) 'CG' for CoolingGaussian.
//' @param WalkType Optional. A string that declares the random walk method: a) 'CDHR' for Coordinate Directions Hit-and-Run, b) 'RDHR' for Random Directions Hit-and-Run or c) 'BW' for Ball Walk. The default walk is  'CDHR'.
//' @param rounding Optional. A boolean parameter for rounding. Default value is FALSE.
//' @param Parameters Optional. A list for the parameters of the algorithms:
//' \itemize{
//' \item{Window }{The length of the sliding window for CG algorithm. The default value is \eqn{500+4dimension^2}.}
//'  \item{C }{A constant for the lower bound of \eqn{variance/mean^2} in schedule annealing of CG algorithm. The default value is \eqn{2}.}
//'  \item{N }{The number of points we sample in each step of schedule annealing in CG algorithm. The default value is \eqn{500C + dimension^2 / 2}.}
//'  \item{ratio }{Parameter of schedule annealing of CG algorithm, larger ratio means larger steps in schedule annealing. The default value is \eqn{1 - 1/dimension}.}
//'  \item{frac }{The fraction of the total error to spend in the first gaussian in CG algorithm. The default value is \eqn{0.1}.}
//'  \item{BW_rad}{The radius for the ball walk. The default value is \eqn{4r/dimension}, where \eqn{r} is the radius of the inscribed ball of the polytope.}
//' }
//'
//' @references \cite{I.Z.Emiris and V. Fisikopoulos,
//' \dQuote{Practical polytope volume approximation,} \emph{ACM Trans. Math. Soft.,} 2014.},
//' @references \cite{B. Cousins and S. Vempala, \dQuote{A practical volume algorithm,} \emph{Springer-Verlag Berlin Heidelberg and The Mathematical Programming Society,} 2015.}
//'
//'
//' @return The approximation of the volume of a convex polytope.
//' @examples
//' # calling SOB algorithm for a H-polytope (2d unit simplex)
//' P = GenSimplex(2,'H')
//' vol = volume(P)
//'
//' # calling CG algorithm for a V-polytope (3d cube)
//' P = GenSimplex(2,'V')
//' vol = volume(P, Algo = list("CG"=TRUE))
//'
//' # calling CG algorithm for a 5-dimensional zonotope defined as the Minkowski sum of 10 segments
//' Z = GenZonotope(2, 4)
//' vol = volume(Z, WalkType = list("method"="hnr", "coordinate"=FALSE, "W"=5), rounding=TRUE)
//' @export
// [[Rcpp::export]]
double volume (Rcpp::Reference P,  Rcpp::Nullable<unsigned int> walk_length = R_NilValue,
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
    unsigned int n = P.field("dimension"), walkL;

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
        throw Rcpp::exception("Unknown walk type!");
    }

    if(!Algo.isNotNull() || Rcpp::as<std::string>(Algo).compare(std::string("SOB"))==0){

        CG = false;

        if(!walk_length.isNotNull()){
            walkL= 10+n/10;
        } else {
            walkL = Rcpp::as<unsigned int>(walk_length);
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

        if (!walk_length.isNotNull()) {
            walkL = 1;
        } else {
            walkL = Rcpp::as<int>(walk_length);
        }

    } else {
        throw Rcpp::exception("Unknown method!");
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

    int type = P.field("type");
    switch(type) {
        case 1: {
            // Hpolytope
            Hpolytope HP;
            HP.init(n, Rcpp::as<MT>(P.field("A")), Rcpp::as<VT>(P.field("b")));
            return generic_volume<Point, NT>(HP, walkL, e, InnerBall, CG, win_len, N, C, ratio, frac, ball_walk, delta,
                                             coordinate, round);
        }
        case 2: {
            // Vpolytope
            Vpolytope VP;
            VP.init(n, Rcpp::as<MT>(P.field("V")), VT::Ones(Rcpp::as<MT>(P.field("V")).rows()));
            return generic_volume<Point, NT>(VP, walkL, e, InnerBall, CG, win_len, N, C, ratio, frac, ball_walk, delta,
                                             coordinate, round);
        }
        case 3: {
            // Zonotope
            zonotope ZP;
            ZP.init(n, Rcpp::as<MT>(P.field("G")), VT::Ones(Rcpp::as<MT>(P.field("G")).rows()));
            return generic_volume<Point, NT>(ZP, walkL, e, InnerBall, CG, win_len, N, C, ratio, frac, ball_walk, delta,
                                             coordinate, round);
        }
        case 4: {
            // Intersection of two V-polytopes
            Vpolytope VP1;
            Vpolytope VP2;
            InterVP VPcVP;
            VP1.init(n, Rcpp::as<MT>(P.field("V1")), VT::Ones(Rcpp::as<MT>(P.field("V1")).rows()));
            VP2.init(n, Rcpp::as<MT>(P.field("V2")), VT::Ones(Rcpp::as<MT>(P.field("V2")).rows()));
            VPcVP.init(VP1, VP2);
            return generic_volume<Point, NT>(VPcVP, walkL, e, InnerBall, CG, win_len, N, C, ratio, frac, ball_walk,
                                             delta,
                                             coordinate, round);
        }
    }

    return 0;
}
