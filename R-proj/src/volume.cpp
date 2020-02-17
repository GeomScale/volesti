// [[Rcpp::depends(BH)]]

// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

#include <Rcpp.h>
#include <RcppEigen.h>
#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include "volume.h"
#include "cooling_balls.h"
#include "cooling_hpoly.h"


template <class Point, class NT, class Polytope>
double generic_volume(Polytope& P, unsigned int walk_step, double e,
                      Rcpp::Nullable<Rcpp::NumericVector> InnerBall, bool CG, bool CB, bool hpoly, unsigned int win_len,
                      unsigned int N, double C, double ratio, double frac,  NT lb, NT ub, NT p, NT alpha,
                      unsigned int NN, unsigned int nu, bool win2, bool ball_walk, double delta, bool cdhr,
                      bool rdhr, bool billiard, double diam, bool rounding, int type)
{
    bool rand_only=false,
         NNN=false,
         birk=false,
         verbose = false;
    unsigned int n_threads=1;
    NT round_val = 1.0, rmax = 0.0;

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
    } else if(type == 2 && CB) {
        if (rounding) {

            //std::cout<<"hello2"<<std::endl;
            InnerB.first = P.get_mean_of_vertices();
            InnerB.second = 0.0;
            vars <NT, RNGType> var2(1, n, 1, n_threads, 0.0, e, 0, 0.0, 0, InnerB.second, 2 * P.get_max_vert_norm(),
                                    rng, urdist, urdist1, -1, verbose, rand_only, rounding, NNN, birk, ball_walk,
                                    cdhr, rdhr, billiard);
            std::pair <NT, NT> res_round = rounding_min_ellipsoid(P, InnerB, var2);
            round_val = res_round.first;

            rounding = false;
            InnerB.second = 0.0;
            InnerB.first = Point(n);
            get_vpoly_center(P);
            rmax = P.get_max_vert_norm();
        } else {
            InnerB.second = 0.0;
            InnerB.first = Point(n);
            get_vpoly_center(P);
            rmax = P.get_max_vert_norm();
        }
    } else {
        // no internal ball or point is given as input
        InnerB = P.ComputeInnerBall();
    }

    //set parameters for billiard and ball walk
    if (billiard && diam < 0.0) {
        P.comp_diam(diam, InnerB.second);
    } else if (ball_walk && delta < 0.0) {
        delta = 4.0 * InnerB.second / NT(n);
    }

    // initialization
    vars<NT, RNGType> var(rnum,n,walk_step,n_threads,0.0,e,0,0.0,0, InnerB.second, diam, rng,urdist,urdist1,
                          delta,verbose,rand_only,rounding,NNN,birk,ball_walk,cdhr,rdhr, billiard);
    NT vol;
    if (CG) {
        vars<NT, RNGType> var2(rnum, n, 10 + n / 10, n_threads, 0.0, e, 0, 0.0, 0, InnerB.second, diam, rng,
                               urdist, urdist1, delta, verbose, rand_only, rounding, NNN, birk, ball_walk, cdhr,
                               rdhr, billiard);
        vars_g<NT, RNGType> var1(n, walk_step, N, win_len, 1, e, InnerB.second, rng, C, frac, ratio, delta, verbose,
                                 rand_only, rounding, NN, birk, ball_walk, cdhr, rdhr);
        vol = volume_gaussian_annealing(P, var1, var2, InnerB);
    } else if (CB) {
        vars_ban <NT> var_ban(lb, ub, p, rmax, alpha, win_len, NN, nu, win2);
        if (!hpoly) {
            vol = vol_cooling_balls(P, var, var_ban, InnerB);
        } else {
            vars_g <NT, RNGType> varg(n, 1, N, 4 * n * n + 500, 1, e, InnerB.second, rng, C, frac, ratio, delta,
                                      verbose, rand_only, false, false, birk, false, true, false);
            vol = vol_cooling_hpoly < HPolytope < Point > > (P, var, var_ban, varg, InnerB);
        }
        if (vol < 0.0) {
            throw Rcpp::exception("Simulated annealing failed! Try to increase the walk length.");
        }
    }else {
        vol = volume(P, var, InnerB);
    }

    return vol * round_val;
}

//' The main function for volume approximation of a convex Polytope (H-polytope, V-polytope or a zonotope)
//'
//' For the volume approximation can be used two algorithms. Either SequenceOfBalls or CoolingGaussian. A H-polytope with \eqn{m} facets is described by a \eqn{m\times d} matrix \eqn{A} and a \eqn{m}-dimensional vector \eqn{b}, s.t.: \eqn{Ax\leq b}. A V-polytope is defined as the convex hull of \eqn{m} \eqn{d}-dimensional points which correspond to the vertices of P. A zonotope is desrcibed by the Minkowski sum of \eqn{m} \eqn{d}-dimensional segments.
//'
//' @param P A convex polytope. It is an object from class (a) Hpolytope or (b) Vpolytope or (c) Zonotope.
//' @param walk_length Optional. The number of the steps for the random walk. The default value is \eqn{\lfloor 10 + d/10\rfloor} for SequenceOfBalls and \eqn{1} otherwise.
//' @param error Optional. Declare the upper bound for the approximation error. The default value is \eqn{1} for SequenceOfBalls and \eqn{0.1} otherwise.
//' @param inner_ball Optional. A \eqn{d+1} vector that contains an inner ball. The first \eqn{d} coordinates corresponds to the center and the last one to the radius of the ball. If it is not given then for H-polytopes the Chebychev ball is computed, for V-polytopes \eqn{d+1} vertices are picked randomly and the Chebychev ball of the defined simplex is computed. For a zonotope that is defined by the Minkowski sum of \eqn{m} segments we compute the maximal \eqn{r} s.t.: \eqn{re_i\in Z} for all \eqn{i=1,\dots ,d}, then the ball centered at the origin with radius \eqn{r/\sqrt{d}} is an inscribed ball.
//' @param algo Optional. A string that declares which algorithm to use: a) \code{'SoB'} for SequenceOfBalls or b) \code{'CG'} for CoolingGaussian or c) \code{'CB'} for cooling bodies.
//' @param random_walk Optional. A string that declares the random walk method: a) \code{'CDHR'} for Coordinate Directions Hit-and-Run, b) \code{'RDHR'} for Random Directions Hit-and-Run, c) \code{'BaW'} for Ball Walk, or \code{'BiW'} for Billiard walk. The default walk is \code{'CDHR'} for H-polytopes and \code{'BiW'} for the other representations.
//' @param rounding Optional. A boolean parameter for rounding. The default value is \code{TRUE} for V-polytopes and \code{FALSE} otherwise.
//' @param parameters Optional. A list for the parameters of the algorithms:
//' \itemize{
//' \item{\code{Window} }{ The length of the sliding window for CG algorithm. The default value is \eqn{500+4dimension^2}.}
//'  \item{\code{C} }{ A constant for the lower bound of \eqn{variance/mean^2} in schedule annealing of CG algorithm. The default value is \eqn{2}.}
//'  \item{\code{M} }{ The number of points we sample in each step of schedule annealing in CG algorithm. The default value is \eqn{500C + dimension^2 / 2}.}
//'  \item{\code{ratio} }{ Parameter of schedule annealing of CG algorithm, larger ratio means larger steps in schedule annealing. The default value is \eqn{1 - 1/dimension}.}
//'  \item{\code{frac} }{ The fraction of the total error to spend in the first gaussian in CG algorithm. The default value is \eqn{0.1}.}
//'  \item{\code{BW_rad} }{ The radius for the ball walk. The default value is \eqn{4r/dimension}, where \eqn{r} is the radius of the inscribed ball of the polytope.}
//'  \item{\code{ub} }{ The lower bound for the ratios in MMC in CB algorithm. The default value is \eqn{0.1}.}
//'  \item{\code{lb} }{ The upper bound for the ratios in MMC in CB algorithm. The default value is \eqn{0.15}.}
//'  \item{\code{N} }{ An integer that controls the number of points \eqn{\nu N} generated in each convex body in annealing schedule of algorithm CB.}
//'  \item{\code{nu} }{ The degrees of freedom for the t-student distribution in t-tests in CB algorithm. The default value is \eqn{10}.}
//'  \item{\code{alpha} }{ The significance level for the t-tests in CB algorithm. The default values is 0.2.}
//'  \item{\code{prob} }{ The probability is used for the empirical confidence interval in ratio estimation of CB algorithm. The default value is \eqn{0.75}.}
//'  \item{\code{hpoly} }{ A boolean parameter to use H-polytopes in MMC of CB algorithm. The default value is \code{FALSE}.}
//'  \item{\code{minmaxW} }{ A boolean parameter to use the sliding window with a minmax values stopping criterion.}
//'  \item{\code{diameter} }{ The diameter of the polytope. It is used to set the maximum length of the trajectory in billiard walk.}
//' }
//'
//' @references \cite{I.Z.Emiris and V. Fisikopoulos,
//' \dQuote{Practical polytope volume approximation,} \emph{ACM Trans. Math. Soft.,} 2014.},
//' @references \cite{A. Chalkis and I.Z.Emiris and V. Fisikopoulos,
//' \dQuote{Practical Volume Estimation by a New Annealing Schedule for Cooling Convex Bodies,} \emph{CoRR, abs/1905.05494,} 2019.},
//' @references \cite{B. Cousins and S. Vempala, \dQuote{A practical volume algorithm,} \emph{Springer-Verlag Berlin Heidelberg and The Mathematical Programming Society,} 2015.}
//'
//'
//' @return The approximation of the volume of a convex polytope.
//' @examples
//' # calling SOB algorithm for a H-polytope (2d unit simplex)
//' P = gen_simplex(2,'H')
//' vol = volume(P)
//'
//' # calling CG algorithm for a V-polytope (3d simplex)
//' P = gen_simplex(2,'V')
//' vol = volume(P, algo = "CG")
//'
//' # calling CG algorithm for a 2-dimensional zonotope defined as the Minkowski sum of 4 segments
//' Z = gen_rand_zonotope(2, 4)
//' vol = volume(Z, random_walk = "RDHR", walk_length = 5)
//' @export
// [[Rcpp::export]]
double volume (Rcpp::Reference P,  Rcpp::Nullable<unsigned int> walk_length = R_NilValue,
                Rcpp::Nullable<double> error = R_NilValue,
                Rcpp::Nullable<Rcpp::NumericVector> inner_ball = R_NilValue,
                Rcpp::Nullable<std::string> algo = R_NilValue,
                Rcpp::Nullable<std::string> random_walk = R_NilValue,
                Rcpp::Nullable<bool> rounding = R_NilValue,
                Rcpp::Nullable<Rcpp::List> parameters = R_NilValue) {

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
    int type = P.field("type");

    bool CG = false, CB = false, cdhr = false, rdhr = false, ball_walk = false, round = false, win2 = false, hpoly = false,
          billiard = false;
    unsigned int win_len = 4*n*n+500, N = 500 * 2 +  n * n / 2, NN = 120 + (n*n)/10, nu = 10;

    NT C = 2.0, ratio = 1.0-1.0/(NT(n)), frac = 0.1, e, delta = -1.0, lb = 0.1, ub = 0.15, p = 0.75, rmax = 0.0,
            alpha = 0.2, diam = -1.0;

    if(!algo.isNotNull()){

        if (type == 2 || type == 3) {
            CB = true;
        } else if (n<=200) {
            CB = true;
        } else {
            CG = true;
        }
        e = (!error.isNotNull()) ? 0.1 : Rcpp::as<NT>(error);
        walkL = (!walk_length.isNotNull()) ? 1 : Rcpp::as<int>(walk_length);

    }else if (Rcpp::as<std::string>(algo).compare(std::string("SOB"))==0){

        walkL = (!walk_length.isNotNull()) ? 10 + n / 10 : Rcpp::as<int>(walk_length);
        e = (!error.isNotNull()) ? 1.0 : Rcpp::as<NT>(error);

    } else if (Rcpp::as<std::string>(algo).compare(std::string("CG"))==0) {

        CG = true;
        e = (!error.isNotNull()) ? 0.1 : Rcpp::as<NT>(error);
        walkL = (!walk_length.isNotNull()) ? 1 : Rcpp::as<int>(walk_length);

    } else if (Rcpp::as<std::string>(algo).compare(std::string("CB")) == 0) {

        CB = true;
        e = (!error.isNotNull()) ? 0.1 : Rcpp::as<NT>(error);
        walkL = (!walk_length.isNotNull()) ? 1 : Rcpp::as<int>(walk_length);

    } else {
        throw Rcpp::exception("Unknown method!");
    }

    if (!random_walk.isNotNull()) {
        if ( type == 1 ){
            cdhr = true;
            if (CB) win_len = 3*n*n+400;
        } else {
            if (CB) {
                billiard = true;
                win_len = 170;
                NN = 125;
            } else {
                rdhr = true;
            }
        }
    }else if (Rcpp::as<std::string>(random_walk).compare(std::string("CDHR")) == 0) {
        cdhr = true;
        if (CB) win_len = 3*n*n+400;
    } else if (Rcpp::as<std::string>(random_walk).compare(std::string("RDHR")) == 0) {
        rdhr = true;
    } else if (Rcpp::as<std::string>(random_walk).compare(std::string("BaW"))==0) {
        ball_walk = true;
    } else if (Rcpp::as<std::string>(random_walk).compare(std::string("BiW"))==0) {
        if (CG) {
            Rcpp::Rcout << "Billiard walk is not supported for CG algorithm. RDHR is used."<<std::endl;
            rdhr = true;
        } else if (SOB) {
            Rcpp::Rcout << "Billiard walk is not supported for SOB algorithm. RDHR is used."<<std::endl;
            rdhr = true;
        } else {
            billiard = true;
            win_len = 170;
            NN = 125;
        }
    }else {
        throw Rcpp::exception("Unknown walk type!");
    }

    if (!rounding.isNotNull() && type == 2){
        round = true;
    } else {
        round = (!rounding.isNotNull()) ? false : Rcpp::as<bool>(rounding);
    }

    if(parameters.isNotNull()) {

        if (Rcpp::as<Rcpp::List>(parameters).containsElementNamed("BW_rad")) {
            delta = Rcpp::as<NT>(Rcpp::as<Rcpp::List>(parameters)["BW_rad"]);
        }
        if (Rcpp::as<Rcpp::List>(parameters).containsElementNamed("C")) {
            C = Rcpp::as<NT>(Rcpp::as<Rcpp::List>(parameters)["C"]);
            N = 500 * ((int) C) + n * n / 2;
        }
        if (Rcpp::as<Rcpp::List>(parameters).containsElementNamed("M")) {
            N = Rcpp::as<int>(Rcpp::as<Rcpp::List>(parameters)["M"]);
        }
        if (Rcpp::as<Rcpp::List>(parameters).containsElementNamed("Window")) {
            win_len = Rcpp::as<int>(Rcpp::as<Rcpp::List>(parameters)["Window"]);
        }
        if (Rcpp::as<Rcpp::List>(parameters).containsElementNamed("frac")) {
            frac = Rcpp::as<NT>(Rcpp::as<Rcpp::List>(parameters)["frac"]);
        }
        if (Rcpp::as<Rcpp::List>(parameters).containsElementNamed("ratio")) {
            ratio = Rcpp::as<NT>(Rcpp::as<Rcpp::List>(parameters)["ratio"]);
        }
        if (Rcpp::as<Rcpp::List>(parameters).containsElementNamed("hpoly")) {
            hpoly = Rcpp::as<bool>(Rcpp::as<Rcpp::List>(parameters)["hpoly"]);
            if ((hpoly && !CB) || (type != 3 && CB && hpoly))
                Rf_warning("flag 'hpoly' can be used to only in MMC of CB algorithm for zonotopes.");
        }
        if (Rcpp::as<Rcpp::List>(parameters).containsElementNamed("lb")) {
            lb = Rcpp::as<NT>(Rcpp::as<Rcpp::List>(parameters)["lb"]);
        }
        if (Rcpp::as<Rcpp::List>(parameters).containsElementNamed("ub")) {
            ub = Rcpp::as<NT>(Rcpp::as<Rcpp::List>(parameters)["ub"]);
        }
        if (Rcpp::as<Rcpp::List>(parameters).containsElementNamed("nu")) {
            nu = Rcpp::as<int>(Rcpp::as<Rcpp::List>(parameters)["nu"]);
        }
        if (Rcpp::as<Rcpp::List>(parameters).containsElementNamed("N")) {
            NN = Rcpp::as<int>(Rcpp::as<Rcpp::List>(parameters)["N"]);
        }
        if (Rcpp::as<Rcpp::List>(parameters).containsElementNamed("minmaxW")) {
            win2 = Rcpp::as<bool>(Rcpp::as<Rcpp::List>(parameters)["minmaxW"]);
        }
        if (Rcpp::as<Rcpp::List>(parameters).containsElementNamed("prob")) {
            p = Rcpp::as<NT>(Rcpp::as<Rcpp::List>(parameters)["prob"]);
        }
        if (Rcpp::as<Rcpp::List>(parameters).containsElementNamed("alpha")) {
            alpha = Rcpp::as<NT>(Rcpp::as<Rcpp::List>(parameters)["alpha"]);
        }
        if (Rcpp::as<Rcpp::List>(parameters).containsElementNamed("diameter")) {
            diam = Rcpp::as<NT>(Rcpp::as<Rcpp::List>(parameters)["diameter"]);
        }
    }

    switch(type) {
        case 1: {
            // Hpolytope
            Hpolytope HP;
            HP.init(n, Rcpp::as<MT>(P.field("A")), Rcpp::as<VT>(P.field("b")));
            return generic_volume<Point, NT>(HP, walkL, e, inner_ball, CG, CB, hpoly, win_len, N, C, ratio, frac, lb, ub, p,
                                             alpha, NN, nu, win2, ball_walk, delta, cdhr, rdhr, billiard, diam, round, type);
        }
        case 2: {
            // Vpolytope
            Vpolytope VP;
            VP.init(n, Rcpp::as<MT>(P.field("V")), VT::Ones(Rcpp::as<MT>(P.field("V")).rows()));
            return generic_volume<Point, NT>(VP, walkL, e, inner_ball, CG, CB, hpoly, win_len, N, C, ratio, frac, lb, ub, p,
                                             alpha, NN, nu, win2, ball_walk, delta, cdhr, rdhr, billiard, diam, round, type);
        }
        case 3: {
            // Zonotope
            zonotope ZP;
            ZP.init(n, Rcpp::as<MT>(P.field("G")), VT::Ones(Rcpp::as<MT>(P.field("G")).rows()));
            return generic_volume<Point, NT>(ZP, walkL, e, inner_ball, CG, CB, hpoly, win_len, N, C, ratio, frac, lb, ub, p,
                                             alpha, NN, nu, win2, ball_walk, delta, cdhr, rdhr, billiard, diam, round, type);
        }
        case 4: {
            // Intersection of two V-polytopes
            Vpolytope VP1;
            Vpolytope VP2;
            InterVP VPcVP;
            VP1.init(n, Rcpp::as<MT>(P.field("V1")), VT::Ones(Rcpp::as<MT>(P.field("V1")).rows()));
            VP2.init(n, Rcpp::as<MT>(P.field("V2")), VT::Ones(Rcpp::as<MT>(P.field("V2")).rows()));
            VPcVP.init(VP1, VP2);
            Rcpp::NumericVector InnerVec(n + 1);
            if (!VPcVP.is_feasible()) throw Rcpp::exception("Empty set!");
            return generic_volume<Point, NT>(VPcVP, walkL, e, inner_ball, CG, CB, hpoly, win_len, N, C, ratio, frac, lb, ub, p,
                                             alpha, NN, nu, win2, ball_walk, delta, cdhr, rdhr, billiard, diam, round, type);
        }
    }

    return 0;
}
