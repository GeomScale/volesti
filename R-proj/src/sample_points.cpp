// [[Rcpp::depends(BH)]]

// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

#include <Rcpp.h>
#include <RcppEigen.h>
#include <chrono>
#include "cartesian_geom/cartesian_kernel.h"
#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include "vars.h"
#include "hpolytope.h"
#include "vpolytope.h"
#include "zpolytope.h"
#include "samplers.h"
#include "gaussian_samplers.h"
#include "sample_only.h"
#include "simplex_samplers.h"
#include "vpolyintersectvpoly.h"


//' Sample points from a convex Polytope (H-polytope, V-polytope or a zonotope) or use direct methods for uniform sampling from the unit or the canonical or an arbitrary \eqn{d}-dimensional simplex and the boundary or the interior of a \eqn{d}-dimensional hypersphere
//'
//' Sample N points with uniform or multidimensional spherical gaussian -centered in an internal point- target distribution.
//' The \eqn{d}-dimensional unit simplex is the set of points \eqn{\vec{x}\in \R^d}, s.t.: \eqn{\sum_i x_i\leq 1}, \eqn{x_i\geq 0}. The \eqn{d}-dimensional canonical simplex is the set of points \eqn{\vec{x}\in \R^d}, s.t.: \eqn{\sum_i x_i = 1}, \eqn{x_i\geq 0}.
//'
//' @param P A convex polytope. It is an object from class (a) Hpolytope or (b) Vpolytope or (c) Zonotope.
//' @param N The number of points that the function is going to sample from the convex polytope. The default value is \eqn{100}.
//' @param distribution Optional. A string that declares the target distribution: a) \code{'uniform'} for the uniform distribution or b) \code{'gaussian'} for the multidimensional spherical distribution. The default target distribution is uniform.
//' @param random_walk Optional. A string that declares the random walk method: a) \code{'CDHR'} for Coordinate Directions Hit-and-Run, b) \code{'RDHR'} for Random Directions Hit-and-Run, c) \code{'BaW'} for Ball Walk or d) \code{'BiW'} for Billiard walk. The default walk is \code{'BiW'} for the uniform distribution or \code{'CDHR'} for the Normal distribution.
//' @param walk_length Optional. The number of the steps for the random walk. The default value is \eqn{1} for \code{'BiW'} and \eqn{\lfloor 10 + d/10\rfloor} otherwise, where \eqn{d} is the dimension that the polytope lies.
//' @param exact A boolean parameter. It should be used for the uniform sampling from the boundary or the interior of a hypersphere centered at the origin or from the unit or the canonical or an arbitrary simplex. The arbitrary simplex has to be given as a V-polytope. For the rest well known convex bodies the dimension has to be declared and the type of body as well as the radius of the hypersphere.
//' @param body A string that declares the type of the body for the exact sampling: a) \code{'unit simplex'} for the unit simplex, b) \code{'canonical simplex'} for the canonical simplex, c) \code{'hypersphere'} for the boundary of a hypersphere centered at the origin, d) \code{'ball'} for the interior of a hypersphere centered at the origin.
//' @param parameters Optional. A list for the parameters of the methods:
//' \itemize{
//' \item{\code{variance} }{ The variance of the multidimensional spherical gaussian. The default value is 1.}
//' \item{\code{dimension} }{ An integer that declares the dimension when exact sampling is enabled for a simplex or a hypersphere.}
//' \item{\code{radius} }{ The radius of the \eqn{d}-dimensional hypersphere. The default value is \eqn{1}.}
//' \item{\code{BW_rad} }{ The radius for the ball walk.}
//' \item{\code{diameter} }{ The diameter of the polytope. It is used to set the maximum length of the trajectory in billiard walk.}
//' }
//' @param InnerPoint A \eqn{d}-dimensional numerical vector that defines a point in the interior of polytope P.
//'
//' @references \cite{R.Y. Rubinstein and B. Melamed,
//' \dQuote{Modern simulation and modeling} \emph{ Wiley Series in Probability and Statistics,} 1998.}
//' @references \cite{A Smith, Noah and W Tromble, Roy,
//' \dQuote{Sampling Uniformly from the Unit Simplex,} \emph{ Center for Language and Speech Processing Johns Hopkins University,} 2004.}
//' @references \cite{Art B. Owen,
//' \dQuote{Monte Carlo theory, methods and examples,} \emph{ Art Owen,} 2009.}
//'
//' @return A \eqn{d\times N} matrix that contains, column-wise, the sampled points from the convex polytope P.
//' @examples
//' # uniform distribution from the 3d unit cube in V-representation using ball walk
//' P = gen_cube(3, 'V')
//' points = sample_points(P, random_walk = "BaW", walk_length = 5)
//'
//' # gaussian distribution from the 2d unit simplex in H-representation with variance = 2
//' A = matrix(c(-1,0,0,-1,1,1), ncol=2, nrow=3, byrow=TRUE)
//' b = c(0,0,1)
//' P = Hpolytope$new(A,b)
//' points = sample_points(P, distribution = "gaussian", parameters = list("variance" = 2))
//'
//' # uniform points from the boundary of a 10-dimensional hypersphere
//' points = sample_points(exact = TRUE, body = "hypersphere", parameters = list("dimension" = 10))
//'
//' # 10000 uniform points from a 2-d arbitrary simplex
//' V = matrix(c(2,3,-1,7,0,0),ncol = 2, nrow = 3, byrow = TRUE)
//' P = Vpolytope$new(V)
//' points = sample_points(P, N = 10000, exact = TRUE)
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix sample_points(Rcpp::Nullable<Rcpp::Reference> P = R_NilValue,
                                  Rcpp::Nullable<unsigned int> N = R_NilValue,
                                  Rcpp::Nullable<std::string> distribution = R_NilValue,
                                  Rcpp::Nullable<std::string> random_walk = R_NilValue,
                                  Rcpp::Nullable<unsigned int> walk_length = R_NilValue,
                                  Rcpp::Nullable<bool> exact = R_NilValue,
                                  Rcpp::Nullable<std::string> body = R_NilValue,
                                  Rcpp::Nullable<Rcpp::List> parameters = R_NilValue,
                                  Rcpp::Nullable<Rcpp::NumericVector> InnerPoint = R_NilValue){

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

    Hpolytope HP;
    Vpolytope VP;
    zonotope ZP;
    InterVP VPcVP;

    int type, dim, numpoints;
    NT radius = 1.0, delta = -1.0, diam = -1.0;
    bool set_mean_point = false, cdhr = false, rdhr = false, ball_walk = false, gaussian = false, billiard = false;
    std::list<Point> randPoints;
    std::pair<Point, NT> InnerBall;

    numpoints = (!N.isNotNull()) ? 100 : Rcpp::as<unsigned int>(N);

    if (exact.isNotNull()) {
        if (P.isNotNull()) {
            type = Rcpp::as<Rcpp::Reference>(P).field("type");
            dim = Rcpp::as<Rcpp::Reference>(P).field("dimension");
            if (Rcpp::as<bool>(exact) && type==2) {
                if (Rcpp::as<MT>(Rcpp::as<Rcpp::Reference>(P).field("V")).rows() == Rcpp::as<MT>(Rcpp::as<Rcpp::Reference>(P).field("V")).cols()+1) {
                    Vpolytope VP;
                    VP.init(dim,
                            Rcpp::as<MT>(Rcpp::as<Rcpp::Reference>(P).field("V")),
                            VT::Ones(Rcpp::as<MT>(Rcpp::as<Rcpp::Reference>(P).field("V")).rows()));
                    Sam_arb_simplex(VP, numpoints, randPoints);
                } else {
                    throw Rcpp::exception("Not a simplex!");
                }
            } else if (Rcpp::as<bool>(exact) && type!=2) {
                throw Rcpp::exception("Not a simplex in V-representation!");
            }
        } else {
            if (!body.isNotNull()) {

                throw Rcpp::exception("Wrong input!");

            } else if (!Rcpp::as<Rcpp::List>(parameters).containsElementNamed("dimension")) {

                throw Rcpp::exception("Wrong input!");

            }
            dim = Rcpp::as<NT>(Rcpp::as<Rcpp::List>(parameters)["dimension"]);
            if (Rcpp::as<Rcpp::List>(parameters).containsElementNamed("radius")) {

                radius = Rcpp::as<NT>(Rcpp::as<Rcpp::List>(parameters)["radius"]);

            }
            if (Rcpp::as<std::string>(body).compare(std::string("hypersphere"))==0) {

                for (unsigned int k = 0; k < numpoints; ++k) {
                    randPoints.push_back(get_point_on_Dsphere<RNGType , Point >(dim, radius));
                }

            } else if (Rcpp::as<std::string>(body).compare(std::string("ball"))==0) {

                for (unsigned int k = 0; k < numpoints; ++k) {
                    randPoints.push_back(get_point_in_Dsphere<RNGType , Point >(dim, radius));
                }

            } else if (Rcpp::as<std::string>(body).compare(std::string("unit simplex"))==0) {

                Sam_Unit<NT, RNGType >(dim, numpoints, randPoints);

            } else if (Rcpp::as<std::string>(body).compare(std::string("canonical simplex"))==0) {

                Sam_Canon_Unit<NT, RNGType >(dim, numpoints, randPoints);

            } else {

                throw Rcpp::exception("Wrong input!");

            }
        }
    } else if (P.isNotNull()) {

        type = Rcpp::as<Rcpp::Reference>(P).field("type");
        dim = Rcpp::as<Rcpp::Reference>(P).field("dimension");
        unsigned int walkL = 10 + dim / 10;

        if (!distribution.isNotNull() || Rcpp::as<std::string>(distribution).compare(std::string("uniform")) == 0) {
            walkL = 3;
            billiard = true;
        } else if (Rcpp::as<std::string>(distribution).compare(std::string("gaussian")) == 0) {
            gaussian = true;
            cdhr = true;
        } else {
            throw Rcpp::exception("Wrong distribution!");
        }

        Point MeanPoint;
        if (InnerPoint.isNotNull()) {
            if (Rcpp::as<Rcpp::NumericVector>(InnerPoint).size() != dim) {
                Rcpp::warning("Internal Point has to lie in the same dimension as the polytope P");
            } else {
                set_mean_point = true;
                MeanPoint = Point(dim, Rcpp::as < std::vector < NT > > (InnerPoint).begin(),
                                  Rcpp::as < std::vector < NT > > (InnerPoint).end());
            }
        }
        if (walk_length.isNotNull()) walkL = Rcpp::as<unsigned int>(walk_length);

        NT a = 0.5;

        if (Rcpp::as<Rcpp::List>(parameters).containsElementNamed("variance"))
            a = 1.0 / (2.0 * Rcpp::as<NT>(Rcpp::as<Rcpp::List>(parameters)["variance"]));

        if (!random_walk.isNotNull()) {
            if (gaussian) {
                cdhr = true;
            } else {
                billiard = true;
            }
        } else if (Rcpp::as<std::string>(random_walk).compare(std::string("CDHR")) == 0) {
            cdhr = true;
        } else if (Rcpp::as<std::string>(random_walk).compare(std::string("RDHR")) == 0) {
            rdhr = true;
        } else if (Rcpp::as<std::string>(random_walk).compare(std::string("BaW")) == 0) {
            if (Rcpp::as<Rcpp::List>(parameters).containsElementNamed("BW_rad"))
                delta = Rcpp::as<NT>(Rcpp::as<Rcpp::List>(parameters)["BW_rad"]);
            ball_walk = true;
        } else if (Rcpp::as<std::string>(random_walk).compare(std::string("BiW")) == 0) {
            if (gaussian) throw Rcpp::exception("Billiard walk can be used only for uniform sampling!");
            billiard = true;
            if (Rcpp::as<Rcpp::List>(parameters).containsElementNamed("diameter"))
                diam = Rcpp::as<NT>(Rcpp::as<Rcpp::List>(parameters)["diameter"]);
        } else {
            throw Rcpp::exception("Unknown walk type!");
        }

        bool rand_only=false,
                NN=false,
                birk=false,
                verbose=false;

        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        // the random engine with this seed
        RNGType rng(seed);
        boost::random::uniform_real_distribution<>(urdist);
        boost::random::uniform_real_distribution<> urdist1(-1,1);

        switch(type) {
            case 1: {
                // Hpolytope
                HP.init(dim, Rcpp::as<MT>(Rcpp::as<Rcpp::Reference>(P).field("A")),
                        Rcpp::as<VT>(Rcpp::as<Rcpp::Reference>(P).field("b")));

                if (!set_mean_point || ball_walk || billiard) {
                    InnerBall = HP.ComputeInnerBall();
                    if (!set_mean_point) MeanPoint = InnerBall.first;
                }
                if (billiard && diam < 0.0) HP.comp_diam(diam, InnerBall.second);
                break;
            }
            case 2: {
                // Vpolytope
                VP.init(dim, Rcpp::as<MT>(Rcpp::as<Rcpp::Reference>(P).field("V")),
                        VT::Ones(Rcpp::as<MT>(Rcpp::as<Rcpp::Reference>(P).field("V")).rows()));

                if (!set_mean_point || ball_walk) {
                    InnerBall = VP.ComputeInnerBall();
                    if (!set_mean_point) MeanPoint = InnerBall.first;
                }
                if (billiard && diam < 0.0) VP.comp_diam(diam, 0.0);
                break;
            }
            case 3: {
                // Zonotope
                ZP.init(dim, Rcpp::as<MT>(Rcpp::as<Rcpp::Reference>(P).field("G")),
                        VT::Ones(Rcpp::as<MT>(Rcpp::as<Rcpp::Reference>(P).field("G")).rows()));

                if (!set_mean_point || ball_walk) {
                    InnerBall = ZP.ComputeInnerBall();
                    if (!set_mean_point) MeanPoint = InnerBall.first;
                }
                if (billiard && diam < 0.0) ZP.comp_diam(diam, 0.0);
                break;
            }
            case 4: {
                // Intersection of two V-polytopes
                Vpolytope VP1;
                Vpolytope VP2;
                VP1.init(dim, Rcpp::as<MT>(Rcpp::as<Rcpp::Reference>(P).field("V1")),
                         VT::Ones(Rcpp::as<MT>(Rcpp::as<Rcpp::Reference>(P).field("V1")).rows()));
                VP2.init(dim, Rcpp::as<MT>(Rcpp::as<Rcpp::Reference>(P).field("V2")),
                         VT::Ones(Rcpp::as<MT>(Rcpp::as<Rcpp::Reference>(P).field("V2")).rows()));
                VPcVP.init(VP1, VP2);

                if (!VPcVP.is_feasible()) throw Rcpp::exception("Empty set!");
                InnerBall = VPcVP.ComputeInnerBall();
                if (!set_mean_point) MeanPoint = InnerBall.first;
                if (billiard && diam < 0.0) {
                    VPcVP.comp_diam(diam, InnerBall.second);
                }
                break;
            }
        }

        if (ball_walk && delta < 0.0) {
            if (gaussian) {
                delta = 4.0 * InnerBall.second / std::sqrt(std::max(NT(1.0), a) * NT(dim));
            } else {
                delta = 4.0 * InnerBall.second / std::sqrt(NT(dim));
            }
        }

        vars<NT, RNGType> var1(1,dim,walkL,1,0.0,0.0,0,0.0,0,InnerBall.second,diam,rng,urdist,urdist1,
                               delta,verbose,rand_only,false,NN,birk,ball_walk,cdhr,rdhr, billiard);
        vars_g<NT, RNGType> var2(dim, walkL, 0, 0, 1, 0, InnerBall.second, rng, 0, 0, 0, delta, verbose,
                                 rand_only, false, NN, birk, ball_walk, cdhr, rdhr);

        switch (type) {
            case 1: {
                sampling_only<Point>(randPoints, HP, walkL, numpoints, gaussian,
                                     a, MeanPoint, var1, var2);
                break;
            }
            case 2: {
                sampling_only<Point>(randPoints, VP, walkL, numpoints, gaussian,
                                     a, MeanPoint, var1, var2);
                break;
            }
            case 3: {
                sampling_only<Point>(randPoints, ZP, walkL, numpoints, gaussian,
                                     a, MeanPoint, var1, var2);
                break;
            }
            case 4: {
                sampling_only<Point>(randPoints, VPcVP, walkL, numpoints, gaussian,
                                     a, MeanPoint, var1, var2);
                break;
            }
        }

    } else {

        throw Rcpp::exception("Wrong input!");

    }

    MT RetMat(dim, numpoints);
    unsigned int jj = 0;

    for (typename std::list<Point>::iterator rpit = randPoints.begin(); rpit!=randPoints.end(); rpit++, jj++)
        RetMat.col(jj) = rpit->getCoefficients();
    return Rcpp::wrap(RetMat);

}
