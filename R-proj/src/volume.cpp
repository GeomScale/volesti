// [[Rcpp::depends(BH)]]

// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 and 2019 program.

#include <Rcpp.h>
#include <RcppEigen.h>
#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include "random_walks/random_walks.hpp"
#include "volume/volume_cooling_gaussians.hpp"
#include "volume/volume_sequence_of_balls.hpp"
#include "volume/volume_cooling_gaussians.hpp"
#include "volume/volume_cooling_balls.hpp"
#include "volume/volume_cooling_hpoly.hpp"


template <typename Polytope,  typename RNGType,  typename NT>
double generic_volume(Polytope& P, RNGType &rng, unsigned int walk_length, NT e,
                      bool CG, bool CB, unsigned int win_len,
                      bool rounding, bool cdhr, bool rdhr, bool ball_walk,
                      bool billiard, int type)
{
    typedef typename Polytope::MT MT;
    typedef typename Polytope::VT VT;
    typedef typename Polytope::PointType Point;

    typedef HPolytope<Point> Hpolytope;

    NT round_val = 1.0;
    unsigned int n = P.dimension();

    if (rounding) {
        std::pair<Point, NT> InnerBall = P.ComputeInnerBall();

        if (type == 1) {
            round_val = round_polytope<CDHRWalk, MT, VT>(P, InnerBall, 10 + 10 * n, rng).second;
        } else {
            round_val = round_polytope<BilliardWalk, MT, VT>(P, InnerBall, 2, rng).second;
        }
    }

    NT vol;
    if (CG) {
        if (cdhr) {
            vol = volume_cooling_gaussians<GaussianCDHRWalk>(P, rng, e, walk_length);
        } else if (rdhr) {
            vol = volume_cooling_gaussians<GaussianRDHRWalk>(P, rng, e, walk_length);
        } else {
            vol = volume_cooling_gaussians<GaussianBallWalk>(P, rng, e, walk_length);
        }
    } else if (CB) {
        if (cdhr) {
            vol = volume_cooling_balls<CDHRWalk>(P, rng, e, walk_length, win_len);
        } else if (rdhr) {
            vol = volume_cooling_balls<RDHRWalk>(P, rng, e, walk_length, win_len);
        } else if (ball_walk) {
            vol = volume_cooling_balls<BallWalk>(P, rng, e, walk_length, win_len);
        } else {
            vol = volume_cooling_balls<BilliardWalk>(P, rng, e, walk_length, win_len);
        }
    } else {
        if (cdhr) {
            vol = volume_sequence_of_balls<CDHRWalk>(P, rng, e, walk_length);
        } else if (rdhr) {
            vol = volume_sequence_of_balls<RDHRWalk>(P, rng, e, walk_length);
        } else if (ball_walk) {
            vol = volume_sequence_of_balls<BallWalk>(P, rng, e, walk_length);
        } else {
            vol = volume_sequence_of_balls<BilliardWalk>(P, rng, e, walk_length);
        }
    }
    vol *= round_val;
    return vol;
}

//' The main function for volume approximation of a convex Polytope (H-polytope, V-polytope or a zonotope)
//'
//' For the volume approximation can be used two algorithms. Either SequenceOfBalls or CoolingGaussian. A H-polytope with \eqn{m} facets is described by a \eqn{m\times d} matrix \eqn{A} and a \eqn{m}-dimensional vector \eqn{b}, s.t.: \eqn{Ax\leq b}. A V-polytope is defined as the convex hull of \eqn{m} \eqn{d}-dimensional points which correspond to the vertices of P. A zonotope is desrcibed by the Minkowski sum of \eqn{m} \eqn{d}-dimensional segments.
//'
//' @param P A convex polytope. It is an object from class (a) Hpolytope or (b) Vpolytope or (c) Zonotope.
//' @param settings Optional. A list that declares which algorithm, random walk and values of parameters to use, as follows:
//' \itemize{
//' \item{\code{algorithm} }{ A string to set the algorithm to use: a) \code{'SoB'} for SequenceOfBalls or b) \code{'CG'} for CoolingGaussian or c) \code{'CB'} for cooling bodies. The defalut algorithm for H-polytopes is \code{'CB'} when \eqn{d\leq 200} and \code{'CG'} when \eqn{d>200}. For the other representations the default algorithm is \code{'CB'}.}
//' \item{\code{error} }{ A numeric value to set the upper bound for the approximation error. The default value is \eqn{1} for \code{'SOB'} and \eqn{0.1} otherwise.}
//' \item{\code{random_walk} }{ A string that declares the random walk method: a) \code{'CDHR'} for Coordinate Directions Hit-and-Run, b) \code{'RDHR'} for Random Directions Hit-and-Run, c) \code{'BaW'} for Ball Walk, or \code{'BiW'} for Billiard walk. The default walk is \code{'CDHR'} for H-polytopes and \code{'BiW'} for the other representations.}
//' \item{\code{walk_length} }{ An integer to set the number of the steps for the random walk. The default value is \eqn{\lfloor 10 + d/10\rfloor} for \code{'SOB'} and \eqn{1} otherwise.}
//' \item{\code{win_len} }{ The length of the sliding window for CG algorithm. The default value is \eqn{500+4dimension^2}.}
//' \item{\code{hpoly} }{ A boolean parameter to use H-polytopes in MMC of CB algorithm. The default value is \code{FALSE}.}
//' }
//' @param rounding Optional. A boolean parameter for rounding. The default value is \code{TRUE} for V-polytopes and \code{FALSE} otherwise.
//' @param seed Optional. A fixed seed for the number generator.
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
//' vol = volume(P, settings = list("algorithm" = "CG"))
//'
//' # calling CG algorithm for a 2-dimensional zonotope defined as the Minkowski sum of 4 segments
//' Z = gen_rand_zonotope(2, 4)
//' vol = volume(Z, settings = list("random_walk" = "RDHR", "walk_length" = 5))
//' @export
// [[Rcpp::export]]
double volume (Rcpp::Reference P,
               Rcpp::Nullable<Rcpp::List> settings = R_NilValue,
               Rcpp::Nullable<bool> rounding = R_NilValue,
               Rcpp::Nullable<double> seed = R_NilValue) {

    typedef double NT;
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT> RNGType;
    typedef HPolytope <Point> Hpolytope;
    typedef VPolytope<Point> Vpolytope;
    typedef Zonotope <Point> zonotope;
    typedef IntersectionOfVpoly<Vpolytope, RNGType> InterVP;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    unsigned int n = P.field("dimension"), walkL, type = P.field("type");

    RNGType rng(n);
    if (seed.isNotNull()) {
        unsigned seed2 = Rcpp::as<double>(seed);
        rng.set_seed(seed2);
    }

    bool CG = false, CB = false, cdhr = false, rdhr = false, ball_walk = false, round = false,
             hpoly = false, billiard = false;
    unsigned int win_len = 4*n*n+500;
    
    NT e;

    if (!Rcpp::as<Rcpp::List>(settings).containsElementNamed("algorithm")) {
        if (type == 2 || type == 3) {
            CB = true;
        } else if (n <= 200) {
            CB = true;
        } else {
            CG = true;
        }
        walkL = (!Rcpp::as<Rcpp::List>(settings).containsElementNamed("walk_length")) ? 1 : Rcpp::as<int>(
                Rcpp::as<Rcpp::List>(settings)["walk_length"]);
        e = (!Rcpp::as<Rcpp::List>(settings).containsElementNamed("error")) ? 0.1 : Rcpp::as<NT>(
                Rcpp::as<Rcpp::List>(settings)["error"]);
    } else if (Rcpp::as<std::string>(Rcpp::as<Rcpp::List>(settings)["algorithm"]).compare(std::string("SOB")) == 0) {

        walkL = (!Rcpp::as<Rcpp::List>(settings).containsElementNamed("walk_length")) ? 10 + n / 10 : Rcpp::as<int>(
                Rcpp::as<Rcpp::List>(settings)["walk_length"]);
        e = (!Rcpp::as<Rcpp::List>(settings).containsElementNamed("error")) ? 1.0 : Rcpp::as<NT>(
                Rcpp::as<Rcpp::List>(settings)["error"]);

    } else if (Rcpp::as<std::string>(Rcpp::as<Rcpp::List>(settings)["algorithm"]).compare(std::string("CG")) == 0) {

        CG = true;
        walkL = (!Rcpp::as<Rcpp::List>(settings).containsElementNamed("walk_length")) ? 1 : Rcpp::as<int>(
                Rcpp::as<Rcpp::List>(settings)["walk_length"]);
        e = (!Rcpp::as<Rcpp::List>(settings).containsElementNamed("error")) ? 0.1 : Rcpp::as<NT>(
                Rcpp::as<Rcpp::List>(settings)["error"]);

    } else if (Rcpp::as<std::string>(Rcpp::as<Rcpp::List>(settings)["algorithm"]).compare(std::string("CB")) == 0) {

        CB = true;
        walkL = (!Rcpp::as<Rcpp::List>(settings).containsElementNamed("walk_length")) ? 1 : Rcpp::as<int>(
                Rcpp::as<Rcpp::List>(settings)["walk_length"]);
        e = (!Rcpp::as<Rcpp::List>(settings).containsElementNamed("error")) ? 0.1 : Rcpp::as<NT>(
                Rcpp::as<Rcpp::List>(settings)["error"]);

    } else {
        throw Rcpp::exception("Unknown method!");
    }


    if (!Rcpp::as<Rcpp::List>(settings).containsElementNamed("random_walk")) {
        if ( type == 1 ){
            cdhr = true;
            if (CB) win_len = 3*n*n+400;
        } else {
            if (CB) {
                billiard = true;
                win_len = 200;
            } else {
                rdhr = true;
            }
        }
    }else if (Rcpp::as<std::string>(Rcpp::as<Rcpp::List>(settings)["random_walk"]).compare(std::string("CDHR")) == 0) {
        cdhr = true;
        if (CB) win_len = 3*n*n+400;
    } else if (Rcpp::as<std::string>(Rcpp::as<Rcpp::List>(settings)["random_walk"]).compare(std::string("RDHR")) == 0) {
        rdhr = true;
        if (CB) win_len = 3*n*n+400;
    } else if (Rcpp::as<std::string>(Rcpp::as<Rcpp::List>(settings)["random_walk"]).compare(std::string("BaW")) == 0) {
        ball_walk = true;
    } else if (Rcpp::as<std::string>(Rcpp::as<Rcpp::List>(settings)["random_walk"]).compare(std::string("BiW")) == 0) {
        if (CG) {
            if (type !=1){
                Rcpp::Rcout << "Billiard walk is not supported for CG algorithm. RDHR is used."<<std::endl;
                rdhr = true;
            } else {
                Rcpp::Rcout << "Billiard walk is not supported for CG algorithm. CDHR is used."<<std::endl;
                cdhr = true;
            }
        } else {
            billiard = true;
            win_len = 200;
        }
    }else {
        throw Rcpp::exception("Unknown walk type!");
    }

    if (e <= 0.0) {
        throw Rcpp::exception("The error parameter has to be a positive number!");
    }

    if (!rounding.isNotNull() && type == 2){
        round = true;
    } else {
        round = (!rounding.isNotNull()) ? false : Rcpp::as<bool>(rounding);
    }

    if (Rcpp::as<Rcpp::List>(settings).containsElementNamed("win_len")) {
        win_len = Rcpp::as<int>(Rcpp::as<Rcpp::List>(settings)["win_len"]);
        if (!CB && !CG) Rf_warning("flag 'win_len' can be used only for CG or CB algorithms.");
    }

    switch(type) {
        case 1: {
            // Hpolytope
            Hpolytope HP;
            HP.init(n, Rcpp::as<MT>(P.field("A")), Rcpp::as<VT>(P.field("b")));
            return generic_volume(HP, rng, walkL, e, CG, CB, win_len, round,
                                             cdhr, rdhr, ball_walk, billiard, type);
        }
        case 2: {
            // Vpolytope
            Vpolytope VP;
            VP.init(n, Rcpp::as<MT>(P.field("V")), VT::Ones(Rcpp::as<MT>(P.field("V")).rows()));
            return generic_volume(VP, rng, walkL, e, CG, CB, win_len, round,
                                             cdhr, rdhr, ball_walk, billiard, type);
        }
        case 3: {
            // Zonotope
            zonotope ZP;
            ZP.init(n, Rcpp::as<MT>(P.field("G")), VT::Ones(Rcpp::as<MT>(P.field("G")).rows()));
            if (Rcpp::as<Rcpp::List>(settings).containsElementNamed("hpoly")) {
                hpoly = Rcpp::as<bool>(Rcpp::as<Rcpp::List>(settings)["hpoly"]);
                if (hpoly && !CB)
                    Rf_warning("flag 'hpoly' can be used to only in MMC of CB algorithm for zonotopes.");
            } else if (ZP.num_of_generators() / ZP.dimension() < 5) {
                hpoly = true;
            } else {
                hpoly = false;
            }
            if (hpoly && CB) {
                if (cdhr) {
                    return volume_cooling_hpoly<CDHRWalk, Hpolytope>(ZP, rng, e, walkL, win_len);
                } else if (rdhr) {
                    return volume_cooling_hpoly<RDHRWalk, Hpolytope>(ZP, rng, e, walkL, win_len);
                } else if (ball_walk) {
                    return volume_cooling_hpoly<BallWalk, Hpolytope>(ZP, rng, e, walkL, win_len);
                } else {
                    return volume_cooling_hpoly<BilliardWalk, Hpolytope>(ZP, rng, e, walkL, win_len);
                }
            }
            return generic_volume(ZP, rng, walkL, e, CG, CB, win_len, round,
                                             cdhr, rdhr, ball_walk, billiard, type);
        }
        case 4: {
            // Intersection of two V-polytopes
            Vpolytope VP1;
            Vpolytope VP2;
            InterVP VPcVP;
            VP1.init(n, Rcpp::as<MT>(P.field("V1")), VT::Ones(Rcpp::as<MT>(P.field("V1")).rows()));
            VP2.init(n, Rcpp::as<MT>(P.field("V2")), VT::Ones(Rcpp::as<MT>(P.field("V2")).rows()));
            VPcVP.init(VP1, VP2);
            if (!VPcVP.is_feasible()) throw Rcpp::exception("Empty set!");
            return generic_volume(VPcVP, rng, walkL, e, CG, CB, win_len, round,
                                             cdhr, rdhr, ball_walk, billiard, type);
        }
    }

    return 0;
}
