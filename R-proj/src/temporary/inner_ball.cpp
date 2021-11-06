// Copyright (c) 2012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

#include <Rcpp.h>
#include <RcppEigen.h>
#include <chrono>
#include "cartesian_geom/cartesian_kernel.h"
#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include "volume/volume_sequence_of_balls.hpp"
#include "preprocess/max_inscribed_ball.hpp"

//' Compute an inscribed ball of a convex polytope
//'
//' For a H-polytope described by a \eqn{m\times d} matrix \eqn{A} and a \eqn{m}-dimensional vector \eqn{b}, s.t.: \eqn{P=\{x\ |\  Ax\leq b\} }, this function computes the largest inscribed ball (Chebychev ball) by solving the corresponding linear program.
//' For both zonotopes and V-polytopes the function computes the minimum \eqn{r} s.t.: \eqn{ r e_i \in P} for all \eqn{i=1, \dots ,d}. Then the ball centered at the origin with radius \eqn{r/ \sqrt{d}} is an inscribed ball.
//'
//' @param P A convex polytope. It is an object from class (a) Hpolytope or (b) Vpolytope or (c) Zonotope or (d) VpolytopeIntersection.
//' @param lpsolve Optional. A boolean variable to compute the Chebychev ball of an H-polytope using the lpsolve library.
//'
//' @return A \eqn{(d+1)}-dimensional vector that describes the inscribed ball. The first \eqn{d} coordinates corresponds to the center of the ball and the last one to the radius.
//'
//' @examples
//' # compute the Chebychev ball of the 2d unit simplex
//' P = gen_cube(10,'H')
//' ball_vec = inner_ball(P)
//'
//' # compute an inscribed ball of the 3-dimensional unit cube in V-representation
//' P = gen_cube(3, 'V')
//' ball_vec = inner_ball(P, lpsolve = TRUE)
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector inner_ball(Rcpp::Reference P, 
                               Rcpp::Nullable<bool> lpsolve = R_NilValue) {

    typedef double NT;
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT> RNGType;
    typedef HPolytope<Point> Hpolytope;
    typedef VPolytope<Point> Vpolytope;
    typedef Zonotope<Point> zonotope;
    typedef IntersectionOfVpoly<Vpolytope, RNGType> InterVP;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    unsigned int n = P.field("dimension"), type = P.field("type");

    std::pair <Point, NT> InnerBall;
    bool lp_solve = (lpsolve.isNotNull()) ? Rcpp::as<bool>(lpsolve) : false;

    switch (type) {
        case 1: {
            // Hpolytope
            Hpolytope HP(n, Rcpp::as<MT>(P.field("A")), Rcpp::as<VT>(P.field("b")));
            if (lp_solve) {
                InnerBall = ComputeChebychevBall<NT, Point>(HP.get_mat(), HP.get_vec());
            } else {
                InnerBall = HP.ComputeInnerBall();
            }
            if (InnerBall.second < 0.0) throw Rcpp::exception("Unable to compute a feasible point.");
            break;
        }
        case 2: {
            // Vpolytope
            Vpolytope VP(n, Rcpp::as<MT>(P.field("V")), VT::Ones(Rcpp::as<MT>(P.field("V")).rows()));
            InnerBall = VP.ComputeInnerBall();
            if (InnerBall.second < 0.0) throw Rcpp::exception("Unable to compute a feasible point.");
            break;
        }
        case 3: {
            // Zonotope
            zonotope ZP(n, Rcpp::as<MT>(P.field("G")), VT::Ones(Rcpp::as<MT>(P.field("G")).rows()));
            InnerBall = ZP.ComputeInnerBall();
            if (InnerBall.second < 0.0) throw Rcpp::exception("Unable to compute a feasible point.");
            break;
        }
        case 4: {
            // Intersection of two V-polytopes
            Vpolytope VP1(n, Rcpp::as<MT>(P.field("V1")), VT::Ones(Rcpp::as<MT>(P.field("V1")).rows()));
            Vpolytope VP2(n, Rcpp::as<MT>(P.field("V2")), VT::Ones(Rcpp::as<MT>(P.field("V2")).rows()));
            InterVP VPcVP(VP1, VP2);
            if (!VPcVP.is_feasible()) throw Rcpp::exception("Empty set!");
            InnerBall = VPcVP.ComputeInnerBall();
            if (InnerBall.second < 0.0) throw Rcpp::exception("Unable to compute a feasible point.");
            break;
        }
    }

    Rcpp::NumericVector vec(n + 1);
    for (unsigned int k = 0; k < n; ++k){
        vec[k] = InnerBall.first[k];
    }

    vec[n] = InnerBall.second;
    return vec;
}
