// [[Rcpp::depends(BH)]]

// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

#include <Rcpp.h>
#include <RcppEigen.h>
#include "volume.h"
#include "ball_ann_vol.h"

// [[Rcpp::export]]
Rcpp::List zono_approx (Rcpp::Reference P, Rcpp::Nullable<bool> fit_ratio = R_NilValue,
                        Rcpp::Nullable<Rcpp::List> Parameters = R_NilValue) {

    typedef double NT;
    typedef Cartesian <NT> Kernel;
    typedef typename Kernel::Point Point;
    typedef boost::mt19937 RNGType;
    typedef HPolytope <Point> Hpolytope;
    typedef Zonotope <Point> zonotope;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
    typedef Eigen::Matrix <NT, Eigen::Dynamic, Eigen::Dynamic> MT;

    int n = Rcpp::as<int>(P.field("dimension")), k = Rcpp::as<MT>(P.field("G")).rows();
    double e = 0.1, delta = -1.0, lb = 0.1, ub = 0.15, p = 0.75, rmax = 0.0, alpha = 0.2;
    int win_len = 2 * n * n + 250, NN = 220 + (n * n) / 10, nu =10, walk_step = 1;
    bool ball_walk = false, verbose = false, cdhr = false, rdhr = true, rounding = false, win2 = false;


    NT ratio = 0.0;

    MT X(n, 2*k);
    X << Rcpp::as<MT>(P.field("G")).transpose(), -Rcpp::as<MT>(P.field("G")).transpose();
    Eigen::JacobiSVD<MT> svd(X*X.transpose(), Eigen::ComputeFullU | Eigen::ComputeFullV);
    MT G(k, 2*n);
    G << Rcpp::as<MT>(P.field("G"))*svd.matrixU(), Rcpp::as<MT>(P.field("G"))*svd.matrixU();
    VT Gred_ii = G.transpose().cwiseAbs().rowwise().sum();
    //std::cout<<Gred_ii<<"\n"<<std::endl;
    MT A(n, 2*n);
    //std::cout<<MT::Identity(n,n)<<"\n"<<std::endl;
    A << -MT::Identity(n,n), MT::Identity(n,n);
    MT Mat(2*n, n+1);

    Mat << Gred_ii, A.transpose()*svd.matrixU().transpose();

   // std::cout<<A<<svd.matrixU().transpose()<<"\n"<<std::endl;

    Hpolytope HP;
    HP.init(n, A.transpose()*svd.matrixU().transpose(), Gred_ii);

    if (fit_ratio.isNotNull() && Rcpp::as<bool>(fit_ratio)) {
        NT vol_red = std::abs(svd.matrixU().determinant());
        for (int i = 0; i < n; ++i) {
            vol_red *= 2.0 * Gred_ii(i);
        }


        zonotope ZP;
        ZP.init(n, Rcpp::as<MT>(P.field("G")), VT::Ones(Rcpp::as<MT>(P.field("G")).rows()));

        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        // the random engine with this seed
        typedef boost::mt19937    RNGType;
        RNGType rng(seed);
        boost::random::uniform_real_distribution<>(urdist);
        boost::random::uniform_real_distribution<> urdist1(-1,1);

        std::pair<Point,NT> InnerBall = ZP.ComputeInnerBall();

        vars<NT, RNGType> var(1, n, walk_step, 1, 0.0, e, 0, 0.0, 0, InnerBall.second, rng,
                               urdist, urdist1, delta, false, false, rounding, false, false, ball_walk, cdhr,rdhr);
        vars_ban <NT> var_ban(lb, ub, p, 0.0, alpha, win_len, NN, nu, win2);

        std::list<Point> randPoints;
        Point q(n);
        rand_point_generator(ZP, q, 1000, 1, randPoints, var);
        int count = 0;
        for (typename std::list<Point>::iterator pit = randPoints.begin();  pit!=randPoints.end(); ++pit) {
            if (HP.is_in(*pit)==-1) count++;
        }
        std::cout<<"points in HP = "<<count<<std::endl;

        NT vol = volesti_ball_ann(ZP, var, var_ban, InnerBall);

        //std::cout<<"vol_pca = "<<vol_red<<" volZ = "<<vol<<std::endl;

        ratio = std::pow(vol_red / vol, 1.0/NT(n));
        std::cout<<"ratio = "<<ratio<<std::endl;
    }

    return Rcpp::List::create(Rcpp::Named("Mat") = Rcpp::wrap(Mat) , Rcpp::Named("fit_ratio") = ratio);

}

