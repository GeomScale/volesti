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
#include "polytopes.h"
#include "samplers.h"
#include "gaussian_samplers.h"
#include "sample_only.h"
#include "simplex_samplers.h"
#include "vpolyintersectvpoly.h"


//' Sample many points from a convex Polytope (H-polytope, V-polytope or a zonotope) or use direct methods for uniform sampling from unit simplex and hypersphere
//'
//' Sample N points from a H or a V-polytope or a zonotope with uniform or spherical gaussian -centered in an internal point- target distribution.
//' The \eqn{d}-dimensional unit simplex is the set of points \eqn{\vec{x}\in \R^d}, s.t.: \eqn{\sum_i x_i\leq 1}, \eqn{x_i\geq 0}. The \eqn{d}-dimensional canonical simplex is the set of points \eqn{\vec{x}\in \R^d}, s.t.: \eqn{\sum_i x_i = 1}, \eqn{x_i\geq 0}.
//'
//' @param P A convex polytope. It is an object from class (a) Hpolytope or (b) Vpolytope or (c) Zonotope.
//' @param N The number of points that the function is going to sample from the convex polytope. The default value is \eqn{100}.
//' @param distribution Optional. A string that declares the target distribution: a) 'uniform' for uniform distribution or b) 'gaussian' for spherical multidimensional distribution. The default target distribution is uniform.
//' @param WalkType Optional. A string that declares the random walk method: a) 'CDHR' for Coordinate Directions Hit-and-Run, b) 'RDHR' for Random Directions Hit-and-Run or c) 'BW' for Ball Walk. The default walk is  'CDHR'.
//' @param walk_length Optional. The number of the steps for the random walk. The default value is \eqn{\lfloor 10 + d/10\rfloor}.
//' @param exact A boolean parameter. It should be used for uniform sampling from the boundary or the interior of a hypersphere centered at the origin or from a unit or an arbitrary simplex. The arbitrary simplex has to be given as a V-polytope. For the rest well known convex bodies it has to be declared the dimension and the type of body (simplex, sphere, ball) as well as the radius of the hypersphere.
//' @param body A string that declares the type of the body for the exact sampling: a) 'unit simplex' for the unit simplex, b) 'canonical simplex' for the canonical simplex, c) 'hypersphere' for the boundary of a hypersphere centered at the origin, d) 'ball' for the interior of a hypersphere centered at the origin.
//' @param Parameters A list for the parameters of the methods:
//' \itemize{
//' \item{variance }{The variance of the spherical multidimensional gaussian. The default value is 1.}
//' \item{dimension }{An integer that declares the dimension when exact sampling is enabled for a simplex or a hypersphere.}
//' \item{radius }{The radius of the \eqn{d}-dimensional hypersphere. Default value is \eqn{1}.}
//' \item{BW_rad }{The radius for the ball walk.}
//' }
//' @param InnerPoint A \eqn{d}-dimensional numerical vector that defines a point in the interior of polytope P.
//'
//' @references \cite{R.Y. Rubinstein and B. Melamed,
//' \dQuote{Modern simulation and modeling} \emph{ Wiley Series in Probability and Statistics,} 1998.}
//' @references \cite{A Smith, Noah and W Tromble, Roy,
//' \dQuote{Sampling Uniformly from the Unit Simplex,} \emph{ Center for Language and Speech Processing Johns Hopkins University,} 2004.}
//' @references \cite{Art B. Owen,
//' \dQuote{Monte Carlo theory, methods and examples,} \emph{ Copyright Art Owen,} 2009-2013.}
//'
//' @return A \eqn{d\times N} matrix that contains, column-wise, the sampled points from the convex polytope.
//' @examples
//' # uniform distribution from a 3d cube in V-representation using ball walk
//' P = GenCube(3, 'V')
//' points = sample_points(P, WalkType = "BW", walk_length = 5)
//'
//' # gaussian distribution from a 2d unit simplex in H-representation with variance = 2
//' A = matrix(c(-1,0,0,-1,1,1), ncol=2, nrow=3, byrow=TRUE)
//' b = c(0,0,1)
//' P = Hpolytope$new(A,b)
//' points = sample_points(P, distribution = "gaussian", Parameters = list("variance" = 2))
//'
//' # uniform points from the boundary of a 10-dimensional hypersphere
//' points = sample_points(exact = TRUE, body = "hypersphere", Parameters = list("dimension" = 10))
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
                                  Rcpp::Nullable<std::string> WalkType = R_NilValue,
                                  Rcpp::Nullable<unsigned int> walk_length = R_NilValue,
                                  Rcpp::Nullable<bool> exact = R_NilValue,
                                  Rcpp::Nullable<std::string> body = R_NilValue,
                                  Rcpp::Nullable<Rcpp::List> Parameters = R_NilValue,
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
    NT radius = 1.0, delta = -1.0;
    bool set_mean_point = false, coordinate = true, ball_walk = false, gaussian = false;
    std::list<Point> randPoints;
    std::pair<Point, NT> InnerBall;

    if (!N.isNotNull()) {
        numpoints = 100;
    } else {
        numpoints = Rcpp::as<unsigned int>(N);
    }

    if (exact.isNotNull()) {
        if (P.isNotNull()) {
            type = Rcpp::as<Rcpp::Reference>(P).field("t");
            if (Rcpp::as<bool>(exact) && type==2) {
                if (Rcpp::as<MT>(Rcpp::as<Rcpp::Reference>(P).field("V")).rows() == Rcpp::as<MT>(Rcpp::as<Rcpp::Reference>(P).field("V")).cols()+1) {
                    Vpolytope VP;
                    VP.init(Rcpp::as<MT>(Rcpp::as<Rcpp::Reference>(P).field("V")).cols(),
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

            } else if (!Rcpp::as<Rcpp::List>(Parameters).containsElementNamed("dimension")) {

                throw Rcpp::exception("Wrong input!");

            }
            dim = Rcpp::as<NT>(Rcpp::as<Rcpp::List>(Parameters)["dimension"]);
            if (Rcpp::as<Rcpp::List>(Parameters).containsElementNamed("radius")) {

                radius = Rcpp::as<NT>(Rcpp::as<Rcpp::List>(Parameters)["radius"]);

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
        unsigned int walkL = 10+dim/10;
        Point MeanPoint;
        if (InnerPoint.isNotNull()) {
            if (Rcpp::as<Rcpp::NumericVector>(InnerPoint).size()!=dim) {
                Rcpp::warning("Internal Point has to lie in the same dimension as the polytope P");
            } else {
                set_mean_point = true;
                MeanPoint = Point( dim , Rcpp::as<std::vector<NT> >(InnerPoint).begin(),
                        Rcpp::as<std::vector<NT> >(InnerPoint).end() );
            }
        }
        if(walk_length.isNotNull()) walkL = Rcpp::as<unsigned int>(walk_length);

        NT a = 0.5;

        if (Rcpp::as<Rcpp::List>(Parameters).containsElementNamed("variance"))
            a = 1.0 / (2.0 * Rcpp::as<NT>(Rcpp::as<Rcpp::List>(Parameters)["variance"]));

        if(!WalkType.isNotNull() || Rcpp::as<std::string>(WalkType).compare(std::string("CDHR"))==0){
            coordinate = true;
            ball_walk = false;
        } else if (Rcpp::as<std::string>(WalkType).compare(std::string("RDHR"))==0) {
            coordinate = false;
            ball_walk = false;
        } else if (Rcpp::as<std::string>(WalkType).compare(std::string("BW"))==0) {
            if (Rcpp::as<Rcpp::List>(Parameters).containsElementNamed("BW_rad")) {
                delta = Rcpp::as<NT>(Rcpp::as<Rcpp::List>(Parameters)["BW_rad"]);
            }
            coordinate = false;
            ball_walk = true;
        } else {
            throw Rcpp::exception("Unknown walk type!");
        }

        if (distribution.isNotNull()) {
            if (Rcpp::as<std::string>(distribution).compare(std::string("gaussian"))==0) {
                gaussian = true;
            } else if(Rcpp::as<std::string>(distribution).compare(std::string("uniform"))!=0) {
                throw Rcpp::exception("Wrong distribution!");
            }
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
        vars<NT, RNGType> var1(1,dim,walkL,1,0.0,0.0,0,0.0,0,InnerBall.second,rng,urdist,urdist1,
                               delta,verbose,rand_only,false,NN,birk,ball_walk,coordinate);
        vars_g<NT, RNGType> var2(dim, walkL, 0, 0, 1, 0, InnerBall.second, rng, 0, 0, 0, delta, false, verbose,
                                 rand_only, false, NN, birk, ball_walk, coordinate);

        switch(type) {
            case 1: {
                // Hpolytope
                HP.init(dim, Rcpp::as<MT>(Rcpp::as<Rcpp::Reference>(P).field("A")),
                        Rcpp::as<VT>(Rcpp::as<Rcpp::Reference>(P).field("b")));

                if (!set_mean_point || ball_walk) {
                    InnerBall = HP.ComputeInnerBall();
                    if (!set_mean_point) MeanPoint = InnerBall.first;
                }
            }
            case 2: {
                // Vpolytope
                VP.init(dim, Rcpp::as<MT>(Rcpp::as<Rcpp::Reference>(P).field("V")),
                        VT::Ones(Rcpp::as<MT>(Rcpp::as<Rcpp::Reference>(P).field("V")).rows()));

                if (!set_mean_point || ball_walk) {
                    InnerBall = VP.ComputeInnerBall();
                    if (!set_mean_point) MeanPoint = InnerBall.first;
                }
            }
            case 3: {
                // Zonotope
                ZP.init(dim, Rcpp::as<MT>(Rcpp::as<Rcpp::Reference>(P).field("G")),
                        VT::Ones(Rcpp::as<MT>(Rcpp::as<Rcpp::Reference>(P).field("G")).rows()));

                if (!set_mean_point || ball_walk) {
                    InnerBall = ZP.ComputeInnerBall();
                    if (!set_mean_point) MeanPoint = InnerBall.first;
                }
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

                if (!set_mean_point || ball_walk) {
                    InnerBall = VP.ComputeInnerBall();
                    if (!set_mean_point) MeanPoint = InnerBall.first;
                }
            }
        }

        if (ball_walk) {
            if (gaussian) {
                delta = 4.0 * InnerBall.second / std::sqrt(std::max(NT(1.0), a) * NT(dim));
            } else {
                delta = 4.0 * InnerBall.second / std::sqrt(NT(dim));
            }
        }

        switch (type) {
            case 1: {
                sampling_only<Point>(randPoints, HP, walkL, numpoints, gaussian,
                                     a, MeanPoint, var1, var2);
            }
            case 2: {
                sampling_only<Point>(randPoints, VP, walkL, numpoints, gaussian,
                                     a, MeanPoint, var1, var2);
            }
            case 3: {
                sampling_only<Point>(randPoints, ZP, walkL, numpoints, gaussian,
                                     a, MeanPoint, var1, var2);
            }
            case 4: {
                sampling_only<Point>(randPoints, VPcVP, walkL, numpoints, gaussian,
                                     a, MeanPoint, var1, var2);
            }
        }

    } else {

        throw Rcpp::exception("Wrong input!");

    }

    Rcpp::NumericMatrix PointSet(dim,numpoints);
    typename std::list<Point>::iterator rpit=randPoints.begin();
    typename std::vector<NT>::iterator qit;
    unsigned int j = 0, i;
    for ( ; rpit!=randPoints.end(); rpit++, j++) {
        qit = (*rpit).iter_begin(); i=0;
        for ( ; qit!=(*rpit).iter_end(); qit++, i++){
            PointSet(i,j)=*qit;
        }
    }
    return PointSet;

}
