// [[Rcpp::depends(BH)]]

// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

#include <Rcpp.h>
#include <RcppEigen.h>
#include "volume.h"
#include "ball_ann_vol.h"
#include "hzono_vol.h"

//' An internal Rccp function for the over-approximation of a zonotope
//'
//' @param Z A zonotope.
//' @param fit_ratio Optional. A boolean parameter to request the computation of the ratio of fitness.
//' @param walk_length Optional. The number of the steps for the random walk. The default value is \eqn{\lfloor 10 + d/10\rfloor} for SequenceOfBalls and \eqn{1} for CoolingGaussian.
//' @param error Optional. Declare the upper bound for the approximation error. The default value is \eqn{1} for SequenceOfBalls and \eqn{0.1} for CoolingGaussian.
//' @param inner_ball Optional. A \eqn{d+1} vector that contains an inner ball. The first \eqn{d} coordinates corresponds to the center and the last one to the radius of the ball. If it is not given then for H-polytopes the Chebychev ball is computed, for V-polytopes \eqn{d+1} vertices are picked randomly and the Chebychev ball of the defined simplex is computed. For a zonotope that is defined by the Minkowski sum of \eqn{m} segments we compute the maximal \eqn{r} s.t.: \eqn{re_i\in Z} for all \eqn{i=1,\dots ,d}, then the ball centered at the origin with radius \eqn{r/\sqrt{d}} is an inscribed ball.
//' @param random_walk Optional. A string that declares the random walk method: a) \code{'CDHR'} for Coordinate Directions Hit-and-Run, b) \code{'RDHR'} for Random Directions Hit-and-Run or c) \code{'BW'} for Ball Walk. The default walk is \code{'CDHR'}.
//' @param rounding Optional. A boolean parameter for rounding. The default value is \code{FALSE}.
//' @param parameters Optional. A list for the parameters of the volume algorithm
//'
//' @section warning:
//' Do not use this function.
//'
//' @return A List that contains a numerical matrix that describes the PCA approximation as a H-polytope and the ratio of fitness.
// [[Rcpp::export]]
Rcpp::List zono_approx (Rcpp::Reference Z, Rcpp::Nullable<bool> fit_ratio = R_NilValue,
                        Rcpp::Nullable<unsigned int> walk_length = R_NilValue,
                        Rcpp::Nullable<double> error = R_NilValue,
                        Rcpp::Nullable<Rcpp::NumericVector> inner_ball = R_NilValue,
                        Rcpp::Nullable<std::string> random_walk = R_NilValue,
                        Rcpp::Nullable<bool> rounding = R_NilValue,
                        Rcpp::Nullable<Rcpp::List> parameters = R_NilValue) {

    typedef double NT;
    typedef Cartesian <NT> Kernel;
    typedef typename Kernel::Point Point;
    typedef boost::mt19937 RNGType;
    typedef HPolytope <Point> Hpolytope;
    typedef Zonotope <Point> zonotope;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
    typedef Eigen::Matrix <NT, Eigen::Dynamic, Eigen::Dynamic> MT;

    int n = Rcpp::as<int>(Z.field("dimension")), k = Rcpp::as<MT>(Z.field("G")).rows();
    double e = 0.1, delta = -1.0, lb = 0.1, ub = 0.15, p = 0.75, rmax = 0.0, alpha = 0.2, diam = -1.0;
    int win_len = 2 * n * n + 250, NN = 220 + (n * n) / 10, nu =10, walkL = 1;
    bool ball_walk = false, verbose = false, cdhr = false, rdhr = false, billiard = false, round = false, win2 = false, hpoly = false;


    NT ratio = 0.0;

    MT X(n, 2*k);
    X << Rcpp::as<MT>(Z.field("G")).transpose(), -Rcpp::as<MT>(Z.field("G")).transpose();
    Eigen::JacobiSVD<MT> svd(X*X.transpose(), Eigen::ComputeFullU | Eigen::ComputeFullV);
    MT G(k, 2*n);
    G << Rcpp::as<MT>(Z.field("G"))*svd.matrixU(), Rcpp::as<MT>(Z.field("G"))*svd.matrixU();
    VT Gred_ii = G.transpose().cwiseAbs().rowwise().sum();
    MT A(n, 2*n);
    A << -MT::Identity(n,n), MT::Identity(n,n);
    MT Mat(2*n, n+1);

    Mat << Gred_ii, A.transpose()*svd.matrixU().transpose();

    Hpolytope HP;
    HP.init(n, A.transpose()*svd.matrixU().transpose(), Gred_ii);

    if (fit_ratio.isNotNull() && Rcpp::as<bool>(fit_ratio)) {
        NT vol_red = std::abs(svd.matrixU().determinant());
        for (int i = 0; i < n; ++i) {
            vol_red *= 2.0 * Gred_ii(i);
        }

        if (error.isNotNull()) e = Rcpp::as<NT>(error);
        if (walk_length.isNotNull()) walkL = Rcpp::as<int>(walk_length);
        if (rounding.isNotNull()) round = Rcpp::as<bool>(rounding);
        if (!random_walk.isNotNull() || Rcpp::as<std::string>(random_walk).compare(std::string("BilW")) == 0) {
            billiard = true;
        } else if(Rcpp::as<std::string>(random_walk).compare(std::string("RDHR")) == 0) {
            rdhr = true;
        } else if (Rcpp::as<std::string>(random_walk).compare(std::string("CDHR")) == 0) {
            cdhr = true;
        } else if (Rcpp::as<std::string>(random_walk).compare(std::string("BW")) == 0) {
            ball_walk = true;
        }else {
            throw Rcpp::exception("Unknown walk type!");
        }

        if(parameters.isNotNull()) {

            if (Rcpp::as<Rcpp::List>(parameters).containsElementNamed("BW_rad")) {
                delta = Rcpp::as<NT>(Rcpp::as<Rcpp::List>(parameters)["BW_rad"]);
            }
            if (Rcpp::as<Rcpp::List>(parameters).containsElementNamed("Window")) {
                win_len = Rcpp::as<int>(Rcpp::as<Rcpp::List>(parameters)["Window"]);
            }  else if (billiard) {
                win_len = 150;
            }
            if (Rcpp::as<Rcpp::List>(parameters).containsElementNamed("lb")) {
                lb = Rcpp::as<double>(Rcpp::as<Rcpp::List>(parameters)["lb"]);
            }
            if (Rcpp::as<Rcpp::List>(parameters).containsElementNamed("ub")) {
                ub = Rcpp::as<double>(Rcpp::as<Rcpp::List>(parameters)["ub"]);
            }
            if (Rcpp::as<Rcpp::List>(parameters).containsElementNamed("nu")) {
                nu = Rcpp::as<int>(Rcpp::as<Rcpp::List>(parameters)["nu"]);
            }
            if (Rcpp::as<Rcpp::List>(parameters).containsElementNamed("N")) {
                NN = Rcpp::as<int>(Rcpp::as<Rcpp::List>(parameters)["N"]);
            } else if (billiard) {
                NN = 125;
            }
            if (Rcpp::as<Rcpp::List>(parameters).containsElementNamed("minmaxW")) {
                win2 = Rcpp::as<bool>(Rcpp::as<Rcpp::List>(parameters)["minmaxW"]);
            }
            if (Rcpp::as<Rcpp::List>(parameters).containsElementNamed("prob")) {
                p = Rcpp::as<double>(Rcpp::as<Rcpp::List>(parameters)["prob"]);
            }
            if (Rcpp::as<Rcpp::List>(parameters).containsElementNamed("alpha")) {
                alpha = Rcpp::as<double>(Rcpp::as<Rcpp::List>(parameters)["alpha"]);
            }
            if (Rcpp::as<Rcpp::List>(parameters).containsElementNamed("hpoly")) {
                hpoly = Rcpp::as<bool>(Rcpp::as<Rcpp::List>(parameters)["hpoly"]);
            }
            if (Rcpp::as<Rcpp::List>(parameters).containsElementNamed("diameter")) {
                diam = Rcpp::as<double>(Rcpp::as<Rcpp::List>(parameters)["diameter"]);
            }
        }

        zonotope ZP;
        ZP.init(n, Rcpp::as<MT>(Z.field("G")), VT::Ones(Rcpp::as<MT>(Z.field("G")).rows()));

        if (billiard && diam < 0.0) ZP.comp_diam(diam);

        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        // the random engine with this seed
        typedef boost::mt19937    RNGType;
        RNGType rng(seed);
        boost::random::uniform_real_distribution<>(urdist);
        boost::random::uniform_real_distribution<> urdist1(-1,1);

        std::pair<Point,NT> InnerB;

        if(inner_ball.isNotNull()) { //if it is given as an input
            // store internal point hat is given as input
            Rcpp::NumericVector InnerVec = Rcpp::as<Rcpp::NumericVector>(inner_ball);
            std::vector<NT> temp_p;
            for (unsigned int j=0; j<n; j++) temp_p.push_back(InnerVec[j]);

            InnerB.first = Point( n , temp_p.begin() , temp_p.end() );
            // store the radius of the internal ball that is given as input
            InnerB.second = InnerVec[n];
        }else {
            // no internal ball or point is given as input
            InnerB = ZP.ComputeInnerBall();
        }

        vars<NT, RNGType> var(1, n, walkL, 1, 0.0, e, 0, 0.0, 0, InnerB.second, diam, rng,
                               urdist, urdist1, delta, false, false, round, false, false, ball_walk, cdhr,rdhr, billiard);
        vars_ban <NT> var_ban(lb, ub, p, 0.0, alpha, win_len, NN, nu, win2);


        NT vol;
        if (!hpoly) {
            vol = volesti_ball_ann(ZP, var, var_ban, InnerB);
        } else {
            vars_g<NT, RNGType> varg(n, 1, 1000 + n * n / 2, 6*n*n+500, 1, e, InnerB.second, rng, 2.0, 0.1,
                                     1.0 - 1.0 / (NT(n)), delta, false, false, false, false, false, false,
                                     false, true, false);
            vol = vol_hzono<Hpolytope> (ZP, var, var_ban, varg, InnerB);
        }
        ratio = std::pow(vol_red / vol, 1.0/NT(n));
    }

    return Rcpp::List::create(Rcpp::Named("Mat") = Rcpp::wrap(Mat) , Rcpp::Named("fit_ratio") = ratio);

}

