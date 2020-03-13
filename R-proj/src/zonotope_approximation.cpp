// [[Rcpp::depends(BH)]]

// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

#include <Rcpp.h>
#include <RcppEigen.h>
#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include "volume.h"
#include "cooling_balls.h"
#include "cooling_hpoly.h"

//' An internal Rccp function for the over-approximation of a zonotope
//'
//' @param Z A zonotope.
//' @param fit_ratio Optional. A boolean parameter to request the computation of the ratio of fitness.
//' @param algo_parameters Optional. A list that declares the values of the parameters of CB algorithm.
//'
//' @section warning:
//' Do not use this function.
//'
//' @return A List that contains a numerical matrix that describes the PCA approximation as a H-polytope and the ratio of fitness.
// [[Rcpp::export]]
Rcpp::List zono_approx (Rcpp::Reference Z,
                        Rcpp::Nullable<bool> fit_ratio = R_NilValue,
                        Rcpp::Nullable<Rcpp::List> algo_parameters = R_NilValue) {

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
    int win_len = 3 * n * n + 400, NN = 220 + (n * n) / 10, nu = 10, walkL = 1;
    bool ball_walk = false, verbose = false, cdhr = false, rdhr = false, billiard = false, round = false, win2 = false,
         hpoly = false, set_mean_point = false;


    NT ratio = std::numeric_limits<double>::signaling_NaN();

    MT X(n, 2 * k);
    X << Rcpp::as<MT>(Z.field("G")).transpose(), -Rcpp::as<MT>(Z.field("G")).transpose();
    Eigen::JacobiSVD <MT> svd(X * X.transpose(), Eigen::ComputeFullU | Eigen::ComputeFullV);
    MT G(k, 2 * n);
    G << Rcpp::as<MT>(Z.field("G")) * svd.matrixU(), Rcpp::as<MT>(Z.field("G")) * svd.matrixU();
    VT Gred_ii = G.transpose().cwiseAbs().rowwise().sum();
    MT A(n, 2 * n);
    A << -MT::Identity(n, n), MT::Identity(n, n);
    MT Mat(2 * n, n + 1);

    Mat << Gred_ii, A.transpose() * svd.matrixU().transpose();

    Hpolytope HP;
    HP.init(n, A.transpose() * svd.matrixU().transpose(), Gred_ii);

    if (fit_ratio.isNotNull() && Rcpp::as<bool>(fit_ratio)) {
        NT vol_red = std::abs(svd.matrixU().determinant());
        for (int i = 0; i < n; ++i) {
            vol_red *= 2.0 * Gred_ii(i);
        }

        walkL = (!Rcpp::as<Rcpp::List>(algo_parameters).containsElementNamed("walk_length")) ? 1 : Rcpp::as<int>(
                Rcpp::as<Rcpp::List>(algo_parameters)["walk_length"]);
        e = (!Rcpp::as<Rcpp::List>(algo_parameters).containsElementNamed("error")) ? 0.1 : Rcpp::as<NT>(
                Rcpp::as<Rcpp::List>(algo_parameters)["error"]);
        round = (!Rcpp::as<Rcpp::List>(algo_parameters).containsElementNamed("rounding")) ? false : Rcpp::as<bool>(
                Rcpp::as<Rcpp::List>(algo_parameters)["rounding"]);
        //if (rounding.isNotNull()) round = Rcpp::as<bool>(rounding);
        if (!Rcpp::as<Rcpp::List>(algo_parameters).containsElementNamed("random_walk")) {
            billiard = true;
            win_len = 170;
            NN = 125;
        } else if (Rcpp::as<std::string>(Rcpp::as<Rcpp::List>(algo_parameters)["random_walk"]).compare(std::string("BiW")) == 0) {
            billiard = true;
            win_len = 170;
            NN = 125;
        } else if (Rcpp::as<std::string>(Rcpp::as<Rcpp::List>(algo_parameters)["random_walk"]).compare(std::string("RDHR")) == 0) {
            rdhr = true;
        } else if (Rcpp::as<std::string>(Rcpp::as<Rcpp::List>(algo_parameters)["random_walk"]).compare(std::string("CDHR")) == 0) {
            cdhr = true;
        } else if (Rcpp::as<std::string>(Rcpp::as<Rcpp::List>(algo_parameters)["random_walk"]).compare(std::string("Baw")) == 0) {
            ball_walk = true;
        } else {
            throw Rcpp::exception("Unknown walk type!");
        }


        if (Rcpp::as<Rcpp::List>(algo_parameters).containsElementNamed("BW_rad")) {
            delta = Rcpp::as<NT>(Rcpp::as<Rcpp::List>(algo_parameters)["BW_rad"]);
        }
        if (Rcpp::as<Rcpp::List>(algo_parameters).containsElementNamed("len_win")) {
            win_len = Rcpp::as<int>(Rcpp::as<Rcpp::List>(algo_parameters)["len_win"]);
        }
        if (Rcpp::as<Rcpp::List>(algo_parameters).containsElementNamed("lb")) {
            lb = Rcpp::as<NT>(Rcpp::as<Rcpp::List>(algo_parameters)["lb"]);
        }
        if (Rcpp::as<Rcpp::List>(algo_parameters).containsElementNamed("ub")) {
            ub = Rcpp::as<NT>(Rcpp::as<Rcpp::List>(algo_parameters)["ub"]);
        }
        if (Rcpp::as<Rcpp::List>(algo_parameters).containsElementNamed("nu")) {
            nu = Rcpp::as<int>(Rcpp::as<Rcpp::List>(algo_parameters)["nu"]);
        }
        if (Rcpp::as<Rcpp::List>(algo_parameters).containsElementNamed("N")) {
            NN = Rcpp::as<int>(Rcpp::as<Rcpp::List>(algo_parameters)["N"]);
        }
        if (Rcpp::as<Rcpp::List>(algo_parameters).containsElementNamed("minmaxW")) {
            win2 = Rcpp::as<bool>(Rcpp::as<Rcpp::List>(algo_parameters)["minmaxW"]);
        }
        if (Rcpp::as<Rcpp::List>(algo_parameters).containsElementNamed("prob")) {
            p = Rcpp::as<NT>(Rcpp::as<Rcpp::List>(algo_parameters)["prob"]);
        }
        if (Rcpp::as<Rcpp::List>(algo_parameters).containsElementNamed("alpha")) {
            alpha = Rcpp::as<NT>(Rcpp::as<Rcpp::List>(algo_parameters)["alpha"]);
        }
        if (Rcpp::as<Rcpp::List>(algo_parameters).containsElementNamed("hpoly")) {
            hpoly = Rcpp::as<bool>(Rcpp::as<Rcpp::List>(algo_parameters)["hpoly"]);
        }
        if (Rcpp::as<Rcpp::List>(algo_parameters).containsElementNamed("L")) {
            diam = Rcpp::as<NT>(Rcpp::as<Rcpp::List>(algo_parameters)["L"]);
        }

        zonotope ZP;
        ZP.init(n, Rcpp::as<MT>(Z.field("G")), VT::Ones(Rcpp::as<MT>(Z.field("G")).rows()));

        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        // the random engine with this seed
        typedef boost::mt19937 RNGType;
        RNGType rng(seed);
        boost::random::uniform_real_distribution<>(urdist);
        boost::random::uniform_real_distribution<> urdist1(-1, 1);


        std::pair <Point, NT> InnerB;
        InnerB.second = -1.0;
        if (Rcpp::as<Rcpp::List>(algo_parameters).containsElementNamed("inner_ball")) {
            if (Rcpp::as<VT>(Rcpp::as<Rcpp::List>(algo_parameters)["inner_ball"]).size() != n + 1) {
                throw Rcpp::exception("Inscribed ball has to lie in the same dimension as the polytope P");
            } else {
                set_mean_point = true;
                VT temp = Rcpp::as<VT>(Rcpp::as<Rcpp::List>(algo_parameters)["inner_ball"]);
                InnerB.first = Point(n, std::vector<NT>(&temp[0], temp.data() + temp.cols() * temp.rows() - 1));
                InnerB.second = temp(n);
                if (InnerB.second <= 0.0)
                    throw Rcpp::exception("The radius of the given inscribed ball has to be a positive number.");
            }
        } else {
            InnerB = ZP.ComputeInnerBall();
        }

        if (billiard && diam < 0.0) {
            ZP.comp_diam(diam, 0.0);
        } else if (ball_walk && delta < 0.0) {
            delta = 4.0 * InnerB.second / NT(n);
        }

        vars <NT, RNGType> var(1, n, walkL, 1, 0.0, e, 0, 0.0, 0, InnerB.second, diam, rng,
                               urdist, urdist1, delta, false, false, round, false, false, ball_walk, cdhr, rdhr,
                               billiard);
        vars_ban <NT> var_ban(lb, ub, p, 0.0, alpha, win_len, NN, nu, win2);

        NT vol;
        if (!hpoly) {
            vol = vol_cooling_balls(ZP, var, var_ban, InnerB);
        } else {
            vars_g <NT, RNGType> varg(n, 1, 1000 + n * n / 2, 6 * n * n + 500, 1, e, InnerB.second, rng, 2.0, 0.1,
                                      1.0 - 1.0 / (NT(n)), delta, false, false, false, false, false, false, true,
                                      false);
            vol = vol_cooling_hpoly<Hpolytope>(ZP, var, var_ban, varg, InnerB);
        }
        ratio = std::pow(vol_red / vol, 1.0 / NT(n));
    }

    return Rcpp::List::create(Rcpp::Named("Mat") = Rcpp::wrap(Mat), Rcpp::Named("fit_ratio") = ratio);

}
