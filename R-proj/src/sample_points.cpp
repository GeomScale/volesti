// [[Rcpp::depends(BH)]]

// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 and 2019 program.

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
#include "vpolyintersectvpoly.h"


//' Sample uniformly or normally distributed points from a convex Polytope (H-polytope, V-polytope or a zonotope).
//'
//' Sample n points with uniform or multidimensional spherical gaussian -with a mode at any point- target distribution.
//'
//' @param P A convex polytope. It is an object from class (a) Hpolytope or (b) Vpolytope or (c) Zonotope or (d) an intersection of two V-polytopes.
//' @param n The number of points that the function is going to sample from the convex polytope.
//' @param random_walk Optional. A list that declares the random walk and some related parameters as follows:
//' \itemize{
//' \item{\code{walk} }{ A string to declare the random walk: i) \code{'CDHR'} for Coordinate Directions Hit-and-Run, ii) \code{'RDHR'} for Random Directions Hit-and-Run, iii) \code{'BaW'} for Ball Walk, iv) \code{'BiW'} for Billiard walk, v) \code{'BCDHR'} boundary sampling by keeping the extreme points of CDHR or vi) \code{'BRDHR'} boundary sampling by keeping the extreme points of RDHR. The default walk is \code{'BiW'} for the uniform distribution or \code{'CDHR'} for the Gaussian distribution.}
//' \item{\code{walk_length} }{ The number of the steps for the random walk. The default value is \eqn{5} for \code{'BiW'} and \eqn{\lfloor 10 + d/10\rfloor} otherwise.}
//' \item{\code{nburns} }{ The number of points to burn before start sampling.}
//' \item{\code{BaW_rad} }{ The radius for the ball walk.}
//' \item{\code{L} }{ The maximum length of the billiard trajectory.}
//' }
//' @param distribution Optional. A list that declares the target density and some related parameters as follows:
//' \itemize{
//' \item{\code{density}}{A string: (a) \code{'uniform'} for the uniform distribution or b) \code{'gaussian'} for the multidimensional spherical distribution. The default target distribution is uniform.}
//' \item{\code{variance} }{ The variance of the multidimensional spherical gaussian. The default value is 1.}
//' \item{\code{StartingPoint} }{ A \eqn{d}-dimensional numerical vector that declares a starting point in the interior of the polytope for the random walk. The default choice is the center of the Chebychev ball.}
//'  \item{\code{mode} }{ A \eqn{d}-dimensional numerical vector that declares the mode of the Gaussian distribution. The default choice is the center of the Chebychev ball.}
//' }
//'
//' @return A \eqn{d\times n} matrix that contains, column-wise, the sampled points from the convex polytope P.
//' @examples
//' # uniform distribution from the 3d unit cube in V-representation using ball walk
//' P = gen_cube(3, 'V')
//' points = sample_points(P, n = 100, random_walk = list("walk" = "BaW", "walk_length" = 5))
//'
//' # gaussian distribution from the 2d unit simplex in H-representation with variance = 2
//' A = matrix(c(-1,0,0,-1,1,1), ncol=2, nrow=3, byrow=TRUE)
//' b = c(0,0,1)
//' P = Hpolytope$new(A,b)
//' points = sample_points(P, n = 100, distribution = list("density" = "gaussian", "variance" = 2))
//'
//' # uniform points from the boundary of a 2-dimensional random H-polytope
//' P = gen_rand_hpoly(2,20)
//' points = sample_points(P, n = 5000, random_walk = list("walk" = "BRDHR"))
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix sample_points(Rcpp::Nullable<Rcpp::Reference> P = R_NilValue,
                                  Rcpp::Nullable<unsigned int> n = R_NilValue,
                                  Rcpp::Nullable<Rcpp::List> random_walk = R_NilValue,
                                  Rcpp::Nullable<Rcpp::List> distribution = R_NilValue){

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

    int type, dim, numpoints, nburns = 0;
    NT radius = 1.0, delta = -1.0, diam = -1.0;
    bool set_mode = false, cdhr = false, rdhr = false, ball_walk = false, gaussian = false,
          billiard = false, boundary = false, set_starting_point = false;
    std::list<Point> randPoints;
    std::pair<Point, NT> InnerBall;
    Point mode(dim);

    numpoints = (!n.isNotNull()) ? 100 : Rcpp::as<unsigned int>(n);

    if (!n.isNotNull()) {
      throw Rcpp::exception("The number of samples is not declared!");
    } else {
        numpoints = Rcpp::as<unsigned int>(n);
        if (numpoints <= 0) throw Rcpp::exception("The number of samples has to be a positice integer!");
    }

    if (!P.isNotNull()) {
        throw Rcpp::exception("No polytope is given as input!");
    }

    type = Rcpp::as<Rcpp::Reference>(P).field("type");
    dim = Rcpp::as<Rcpp::Reference>(P).field("dimension");
    int walkL = 10 + dim / 10;

    if (!distribution.isNotNull() || !Rcpp::as<Rcpp::List>(distribution).containsElementNamed("density")) {
        billiard = true;
    } else if (
            Rcpp::as<std::string>(Rcpp::as<Rcpp::List>(distribution)["density"]).compare(std::string("uniform")) == 0) {
            billiard = true;
    } else if (
            Rcpp::as<std::string>(Rcpp::as<Rcpp::List>(distribution)["density"]).compare(std::string("gaussian")) == 0) {
            gaussian = true;
    } else {
        throw Rcpp::exception("Wrong distribution!");
    }

    Point StartingPoint;
    if (Rcpp::as<Rcpp::List>(distribution).containsElementNamed("starting_point")) {
        if (Rcpp::as<VT>(Rcpp::as<Rcpp::List>(distribution)["starting_point"]).size() != dim) {
            throw Rcpp::exception("Starting Point has to lie in the same dimension as the polytope P");
        } else {
            set_starting_point = true;
            VT temp = Rcpp::as<VT>(Rcpp::as<Rcpp::List>(distribution)["starting_point"]);
            StartingPoint = Point(dim, std::vector<NT>(&temp[0], temp.data() + temp.cols() * temp.rows()));
        }
    }

    if (Rcpp::as<Rcpp::List>(distribution).containsElementNamed("mode")) {
        if (!gaussian) throw Rcpp::exception("Mode is given only for Gaussian sampling!");
        if (Rcpp::as<VT>(Rcpp::as<Rcpp::List>(distribution)["mode"]).size() != dim) {
            throw Rcpp::exception("Mode has to be a point in the same dimension as the polytope P");
        } else {
            set_mode = true;
            VT temp = Rcpp::as<VT>(Rcpp::as<Rcpp::List>(distribution)["mode"]);
            mode = Point(dim, std::vector<NT>(&temp[0], temp.data() + temp.cols() * temp.rows()));
        }
    }

    NT a = 0.5;
    if (Rcpp::as<Rcpp::List>(distribution).containsElementNamed("variance")) {
        a = 1.0 / (2.0 * Rcpp::as<NT>(Rcpp::as<Rcpp::List>(distribution)["variance"]));
        if (!gaussian) {
            Rcpp::warning("The variance can be set only for Gaussian sampling!");
        } else if (a <= 0.0) {
            throw Rcpp::exception("The variance has to be positive!");
        }
    }

    if (!random_walk.isNotNull() || !Rcpp::as<Rcpp::List>(random_walk).containsElementNamed("walk")) {
        if (gaussian) {
            cdhr = true;
        } else {
            billiard = true;
            walkL = 5;
        }
    } else if (Rcpp::as<std::string>(Rcpp::as<Rcpp::List>(random_walk)["walk"]).compare(std::string("CDHR")) == 0) {
        cdhr = true;
    } else if (Rcpp::as<std::string>(Rcpp::as<Rcpp::List>(random_walk)["walk"]).compare(std::string("RDHR")) == 0) {
        rdhr = true;
    } else if (Rcpp::as<std::string>(Rcpp::as<Rcpp::List>(random_walk)["walk"]).compare(std::string("BaW")) == 0) {
        if (Rcpp::as<Rcpp::List>(random_walk).containsElementNamed("BaW_rad")) {
            delta = Rcpp::as<NT>(Rcpp::as<Rcpp::List>(random_walk)["BaW_rad"]);
        }
        ball_walk = true;
    } else if (Rcpp::as<std::string>(Rcpp::as<Rcpp::List>(random_walk)["walk"]).compare(std::string("BiW")) == 0) {
        if (gaussian) throw Rcpp::exception("Billiard walk can be used only for uniform sampling!");
        billiard = true;
        walkL = 5;
        if (Rcpp::as<Rcpp::List>(random_walk).containsElementNamed("L")) {
            diam = Rcpp::as<NT>(Rcpp::as<Rcpp::List>(random_walk)["L"]);
        }
    } else if (Rcpp::as<std::string>(Rcpp::as<Rcpp::List>(random_walk)["walk"]).compare(std::string("BRDHR")) == 0) {
        if (gaussian) throw Rcpp::exception("Gaussian sampling from the boundary is not supported!");
        rdhr = true;
        boundary = true;
    }else if (Rcpp::as<std::string>(Rcpp::as<Rcpp::List>(random_walk)["walk"]).compare(std::string("BCDHR")) == 0) {
        if (gaussian) throw Rcpp::exception("Gaussian sampling from the boundary is not supported!");
        cdhr = true;
        boundary = true;
    }else {
        throw Rcpp::exception("Unknown walk type!");
    }

    if (Rcpp::as<Rcpp::List>(random_walk).containsElementNamed("walk_length")) {
        walkL = Rcpp::as<int>(Rcpp::as<Rcpp::List>(random_walk)["walk_length"]);
        if (walkL <= 0) {
            throw Rcpp::exception("The walk length has to be a positive integer!");
        }
    }

    if (Rcpp::as<Rcpp::List>(random_walk).containsElementNamed("nburns")) {
        nburns = Rcpp::as<int>(Rcpp::as<Rcpp::List>(random_walk)["nburns"]);
        if (nburns < 0) {
            throw Rcpp::exception("The number of points to burn before sampling has to be a positive integer!");
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

    switch(type) {
        case 1: {
            // Hpolytope
            HP.init(dim, Rcpp::as<MT>(Rcpp::as<Rcpp::Reference>(P).field("A")),
                    Rcpp::as<VT>(Rcpp::as<Rcpp::Reference>(P).field("b")));

            if (!set_starting_point || (!set_mode && gaussian) || ball_walk || billiard) {
                InnerBall = HP.ComputeInnerBall();
                if (!set_starting_point) StartingPoint = InnerBall.first;
                if (!set_mode && gaussian) mode = InnerBall.first;
            }
            if (HP.is_in(StartingPoint) == 0) {
                throw Rcpp::exception("The given point is not in the interior of the polytope!");
            }
            if (billiard && diam < 0.0) HP.comp_diam(diam, InnerBall.second);
            HP.normalize();
            if (gaussian) {
                StartingPoint = StartingPoint - mode;
                HP.shift(Eigen::Map<VT>(&mode.get_coeffs()[0], mode.dimension()));
            }
            break;
        }
        case 2: {
            // Vpolytope
            VP.init(dim, Rcpp::as<MT>(Rcpp::as<Rcpp::Reference>(P).field("V")),
                    VT::Ones(Rcpp::as<MT>(Rcpp::as<Rcpp::Reference>(P).field("V")).rows()));

            if (!set_starting_point || (!set_mode && gaussian) || ball_walk) {
                InnerBall = VP.ComputeInnerBall();
                if (!set_starting_point) StartingPoint = InnerBall.first;
                if (!set_mode && gaussian) mode = InnerBall.first;
            }
            if (VP.is_in(StartingPoint) == 0)
                throw Rcpp::exception("The given point is not in the interior of the polytope!");
            if (billiard && diam < 0.0) VP.comp_diam(diam);
            if (gaussian) {
                StartingPoint = StartingPoint - mode;
                VP.shift(Eigen::Map<VT>(&mode.get_coeffs()[0], mode.dimension()));
            }
            break;
        }
        case 3: {
            // Zonotope
            ZP.init(dim, Rcpp::as<MT>(Rcpp::as<Rcpp::Reference>(P).field("G")),
                    VT::Ones(Rcpp::as<MT>(Rcpp::as<Rcpp::Reference>(P).field("G")).rows()));

            if (!set_starting_point || (!set_mode && gaussian) || ball_walk) {
                InnerBall = ZP.ComputeInnerBall();
                if (!set_starting_point) StartingPoint = InnerBall.first;
                if (!set_mode && gaussian) mode = InnerBall.first;
            }
            if (ZP.is_in(StartingPoint) == 0)
                throw Rcpp::exception("The given point is not in the interior of the polytope!");
            if (billiard && diam < 0.0) ZP.comp_diam(diam);
            if (gaussian) {
                StartingPoint = StartingPoint - mode;
                ZP.shift(Eigen::Map<VT>(&mode.get_coeffs()[0], mode.dimension()));
            }
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
            if (!set_starting_point) StartingPoint = InnerBall.first;
            if (!set_mode && gaussian) mode = InnerBall.first;
            if (VPcVP.is_in(StartingPoint) == 0)
                throw Rcpp::exception("The given point is not in the interior of the polytope!");
            if (billiard && diam < 0.0) {
                VPcVP.comp_diam(diam, InnerBall.second);
            }
            if (gaussian) {
                StartingPoint = StartingPoint - mode;
                VPcVP.shift(Eigen::Map<VT>(&mode.get_coeffs()[0], mode.dimension()));
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
                                 a, boundary, StartingPoint, nburns, var1, var2);
            break;
        }
        case 2: {
            sampling_only<Point>(randPoints, VP, walkL, numpoints, gaussian,
                                 a, boundary, StartingPoint, nburns, var1, var2);
            break;
        }
        case 3: {
            sampling_only<Point>(randPoints, ZP, walkL, numpoints, gaussian,
                                 a, boundary, StartingPoint, nburns, var1, var2);
            break;
        }
        case 4: {
            sampling_only<Point>(randPoints, VPcVP, walkL, numpoints, gaussian,
                                 a, boundary, StartingPoint, nburns, var1, var2);
            break;
        }
    }

    if (numpoints % 2 == 1 && boundary) numpoints--;
    MT RetMat(dim, numpoints);
    unsigned int jj = 0;

    for (typename std::list<Point>::iterator rpit = randPoints.begin(); rpit!=randPoints.end(); rpit++, jj++) {
        if (gaussian) {
            RetMat.col(jj) = Eigen::Map<VT>(&(*rpit).get_coeffs()[0], (*rpit).dimension()) +
                             Eigen::Map<VT>(&mode.get_coeffs()[0], mode.dimension());
        } else {
            RetMat.col(jj) = Eigen::Map<VT>(&(*rpit).get_coeffs()[0], (*rpit).dimension());
        }
    }
    return Rcpp::wrap(RetMat);

}
